"""
Common functions for use
"""

"""Functions for maps"""

# Function for single scatter plots
function map_points(point_df, var_mean, label_name, range, colorscheme, title_name, res)
    # point_df: stations with coordinates
    # var_mean: values to be mapped
    # res: resolution of the plot
    states = download(
        "https://raw.githubusercontent.com/PublicaMundi/MappingAPI/master/data/geojson/us-states.json",
    )
    states_geo = GeoJSON.read(read(states, String))

    Δ = 0.25
    lonlims = (minimum(point_df[!, :lon]) - Δ, maximum(point_df[!, :lon]) + Δ)
    latlims = (minimum(point_df[!, :lat]) - Δ, maximum(point_df[!, :lat]) + Δ)

    f = Figure(resolution=res, fontsize=22)
    ax = GeoAxis(
        f[1, 1];
        dest="+proj=eqc",
        title=title_name,
        lonlims=lonlims,
        latlims=latlims
    )

    poly!(ax, states_geo; color=:white, strokewidth=1)

    h1 = CairoMakie.scatter!(ax, point_df[!, :lon], point_df[!, :lat]; color=var_mean, colormap=colorscheme, colorrange=range, markersize=20, strokewidth=0, strokecolor=:grey17)

    Colorbar(f[1, 2], h1, label=label_name)
    poly!(ax, states_geo, edgecolor=:gray, strokewidth=1, color=:transparent)
    return f
end


# Function for subplots of mapped points
function map_points_subplots(points_df, all_data, colorschemes, rows, cols, row_names, res, row_hs, title_names, ranges; diff_coord=false, bar_all=true, diff_colscheme=false, bar_name=false)
    # bar_all: if would need color bar for all subplots
    # diff_colscheme: if different color schemes for the color bars
    # point_df: stations with coordinates
    # all_data: a vector with all the plotting vectors
    # diff_coord: if plotting for the same points
    states = download(
        "https://raw.githubusercontent.com/PublicaMundi/MappingAPI/master/data/geojson/us-states.json",
    )
    states_geo = GeoJSON.read(read(states, String))

    Δ = 0.25

    f = Figure(resolution=res, fontsize=22, layout=GridLayout(row_heights=row_hs))

    point_df = points_df
    colorscheme = colorschemes

    for i in 1:length(all_data)
        if diff_colscheme == true
            colorscheme = colorschemes[i]
        end
        if diff_coord == true
            point_df = points_df[i]
        end
        lonlims = (minimum(point_df[!, :lon]) - Δ, maximum(point_df[!, :lon]) + Δ)
        latlims = (minimum(point_df[!, :lat]) - Δ, maximum(point_df[!, :lat]) + Δ)
        ax = GeoAxis(
            f[rows[i], cols[i]];
            dest="+proj=eqc",
            title=title_names[i],
            lonlims=lonlims,
            latlims=latlims,
            ylabel=row_names[i]
        )
        poly!(ax, states_geo; color=:white, strokewidth=1)
        h1 = CairoMakie.scatter!(ax, point_df[!, :lon], point_df[!, :lat]; color=all_data[i], colormap=colorscheme, colorrange=ranges[i], markersize=15, strokewidth=0, strokecolor=:grey17)

        if bar_name == false
            b_name = ""
        else
            b_name = bar_name[i]
        end

        if bar_all == true
            Colorbar(f[rows[i], cols[i]+1], h1, label=b_name, height=Relative(0.9))
        else
            if cols[i] == maximum(cols)
                Colorbar(f[rows[i], cols[i]+1], h1, label=b_name, height=Relative(0.9))
            end
        end

        poly!(ax, states_geo, edgecolor=:gray, strokewidth=1, color=:transparent)
    end
    return f
end


"""Functions for return level estimates"""

# nonstationary return level given GEV parameters
function rl_estimate_N0(p, x, μ0, μ_beta, logσ0, logσ_beta, ξ)
    μ = μ0 + x * μ_beta
    σ = exp(logσ0 + x * logσ_beta)
    return quantile(GeneralizedExtremeValue(μ, σ, ξ), p)
end

# estimate for the Spatially Varying Covariate Model (mean and std for MCMC posteriors)
function rl_estimate_full(point_df, μ_beta, logσ_beta, μ0_posterior, logσ0_posterior, ξ_posterior, x, p)
    mean_stations = zeros(size(point_df)[1])
    std_stations = zeros(size(point_df)[1])
    for k in 1:size(point_df)[1]
        μ_beta_k = μ_beta[:, k]
        logσ_beta_k = logσ_beta[:, k]
        μ0_k = μ0_posterior[:, k]
        logσ0_k = logσ0_posterior[:, k]
        rl_k = rl_estimate1.(p, x, μ0_k, μ_beta_k, logσ0_k, logσ_beta_k, ξ_posterior)
        mean_stations[k] = mean(rl_k)
        std_stations[k] = std(rl_k)
    end
    return mean_stations, std_stations
end

# estimate for the Nonpooled Nonstationary Model (mean and std for MCMC posteriors)
function rl_estimate_N(point_df, μ_beta, logσ_beta, μ0_posterior, logσ0_posterior, ξ_posterior, x, p)
    mean_stations = zeros(size(point_df)[1])
    std_stations = zeros(size(point_df)[1])
    for k in 1:size(point_df)[1]
        μ_beta_k = μ_beta[:, k]
        logσ_beta_k = logσ_beta[:, k]
        μ0_k = μ0_posterior[:, k]
        logσ0_k = logσ0_posterior[:, k]
        ξ_k = ξ_posterior[:, k]
        rl_k = rl_estimate1.(p, x, μ0_k, μ_beta_k, logσ0_k, logσ_beta_k, ξ_k)
        mean_stations[k] = mean(rl_k)
        std_stations[k] = std(rl_k)
    end
    return mean_stations, std_stations
end

# estimate for the Pooled Stationary Model (mean and std for MCMC posteriors)
function rl_etimate_S(point_df, mu, logs, xi, p)
    mean_stations = zeros(size(point_df)[1])
    std_stations = zeros(size(point_df)[1])
    for k in 1:size(point_df)[1]
        mu_k = mu[:, k]
        logs_k = logs[:, k]
        rl_k = quantile.(GeneralizedExtremeValue.(mu_k, exp.(logs_k), xi), p)
        mean_stations[k] = mean(rl_k)
        std_stations[k] = std(rl_k)
    end
    return mean_stations, std_stations
end