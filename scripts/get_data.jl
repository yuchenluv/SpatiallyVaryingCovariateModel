"""
Download and save all required input data, 
including 1d annual maximum precipitation from GHCN, CO2 concentrations from Law Dome and Mauna Loa
"""


using CategoricalArrays: categorical
using CSV
using DataFrames
using DataFramesMeta
using Dates
using DrWatson
using GHCNData
using HTTP
using ProgressMeter
using StatsBase


"""
Get the GHCN inventories, filtered to only include stations within the bounds of a square bounding box
given by values of `params`
"""
function produce_station_inventory(params)
    station_inventory = @chain GHCNData.load_inventories() begin
        @rsubset :ELEMENT == "PRCP"
        @rsubset params.latmin <= :LATITUDE <= params.latmax
        @rsubset params.lonmin <= :LONGITUDE <= params.lonmax
        @transform :NYEARS = :LASTYEAR - :FIRSTYEAR
        @rsubset :NYEARS >= 30
        @select :ID :LATITUDE :LONGITUDE
        rename(:ID => :stnid, :LATITUDE => :lat, :LONGITUDE => :lon)
    end
    return Dict("station_inventory" => station_inventory)
end

function get_station_inventory(params)
    needed_params = params[[:latmin, :lonmin, :latmax, :lonmax]]
    return produce_or_load(
        datadir("processed", "raw_1d", "station_inventory"), needed_params, produce_station_inventory
    )[1]["station_inventory"]
end


"""
Get the GHCN data for a single station as a DataFrame
See https://www1.ncdc.noaa.gov/pub/data/ghcn/daily/readme.txt for a description of the invalid QFLAGS
"""
function stnid_to_df(id)
    invalid_qflags = ["D", "G", "I", "K", "L", "M", "N", "O", "R", "S", "T", "W", "X", "Z"]
    df = @chain id begin
        GHCNData.load_data_file
        GHCNData.convert_to_time_series
        rename(:DATE => :date, :PRCP_QFLAG => :flag)
        dropmissing(:PRCP_VALUE) # toss missing prcp data
        @rsubset !in(:flag, invalid_qflags) # toss any data points that have been flagged
        @transform :flag = categorical(replace(:flag, " " => missing)) # categorical is easier to store
        @transform :prcp = :PRCP_VALUE * 0.1#u"mm" # add units to avoid confusion
        @transform :stnid = id
        @select :stnid :date :prcp :flag # keep only what we need
    end
    return df
end


"""
Produce the daily precip time series for all stations
"""
function produce_daily_precip(params)
    station_inventory = get_station_inventory(params)
    stnids = station_inventory.stnid
    daily_precip = @showprogress "Reading daily precipitation files" map(stnid_to_df, stnids)
    return Dict("daily_precip" => daily_precip)
end

function get_daily_precip(params)
    needed_params = params[[:latmin, :lonmin, :latmax, :lonmax]]
    return produce_or_load(
        datadir("processed", "raw_1d", "daily_precip"), needed_params, produce_daily_precip
    )[1]["daily_precip"]
end


"""
Convert daily time series to annual maximum series
"""
function daily_to_annmax(df_daily, min_obs_per_year, min_years)
    stnid = first(df_daily[!, :stnid])
    df2 = @chain df_daily begin
        @rtransform :month = Dates.month(:date)
        @rtransform :year = Dates.year(:date)
        dropmissing(:prcp)
        groupby(:year)
        @transform :N = length(:date)
        @rsubset :N >= min_obs_per_year
    end
    nrow(df2) == 0 && return Nothing
    annual = @chain df2 begin
        groupby(:year)
        @combine begin
            :prcp = maximum(:prcp)
        end
        rename(:prcp => stnid)
    end
    nrow(annual) <= min_years && return Nothing
    return annual
end

function produce_annmax_precip(params)
    daily_precip = get_daily_precip(params)
    annual_dfs = [
        daily_to_annmax(df, params[:min_obs_per_year], params[:min_years]) for
        df in daily_precip
    ]
    annual_dfs = [df for df in annual_dfs if (df != Nothing)]
    annual = outerjoin(annual_dfs...; on=:year)
    annual = sort(annual, :year)
    return Dict("annual" => annual)
end

function get_annmax_precip(params)
    needed_params = params[[
        :latmin, :lonmin, :latmax, :lonmax, :min_obs_per_year, :min_years
    ]]
    return produce_or_load(
        datadir("processed", "raw_1d", "AMS"), needed_params, produce_annmax_precip;
    )[1]["annual"]
end


"""
Combining Law Dome and Mauna Loa CO2 data
Largest difference between these two sets of data is about 1.71 ppm
"""
function produce_combined_CO2(params)

    # 1958 ~ 2023 Mauna Loa
    url_ML = "https://gml.noaa.gov/webdata/ccgg/trends/co2/co2_mm_mlo.csv"
    monthly_ML = DataFrame(CSV.File(download(url_ML); header=41))
    annual_ML = @chain monthly_ML begin
        groupby(:year)
        @combine :CO2 = mean(:average)
        @rtransform :log_CO2_ML = log(:CO2)
        @select :year :log_CO2_ML
    end

    # 154 ~ 1996 Law Dome
    url_LD = "https://www.ncei.noaa.gov/pub/data/paleo/icecore/antarctica/law/law2018splines.txt"
    annual_LD = DataFrame(CSV.File(HTTP.get(url_LD).body; header=118, delim="\t", ignorerepeated=true))[:, 1:2]
    rename!(annual_LD, [:year, :CO2])
    annual_LD = @chain annual_LD begin
        @rtransform :log_CO2_LD = log(:CO2)
        @select :year :log_CO2_LD
        @rsubset :year >= 1849
    end

    CO2_merged = outerjoin(annual_ML, annual_LD, on=:year, makeunique=true)
    CO2_merged = sort(CO2_merged, :year)
    CO2_merged.log_CO2_ML .= coalesce.(CO2_merged.log_CO2_ML, CO2_merged.log_CO2_LD)
    CO2_merged.log_CO2_LD .= coalesce.(CO2_merged.log_CO2_ML, CO2_merged.log_CO2_LD)
    CO2_merged.log_CO2 = (CO2_merged.log_CO2_ML + CO2_merged.log_CO2_LD) / 2
    select!(CO2_merged, Not([:log_CO2_ML, :log_CO2_LD]))

    return Dict("CO2" => CO2_merged)
end

function get_CO2(params)
    needed_params = Dict()
    return produce_or_load(
        datadir("processed", "predictors"), needed_params, produce_combined_CO2; prefix="CO2"
    )[1]["CO2"]
end
