"""
Common functions for plotting
"""

"""Functions for maps"""

# function for plotting single axes that can be used for
# main references
# https://docs.makie.org/v0.21/tutorials/layout-tutorial
# https://docs.makie.org/dev/reference/blocks/gridlayout

function point_axis(ax, title_name, lonlims, latlims, ylabel_name, width, height, states_geo, point_df, var, colorscheme, range, with_xticks, with_yticks, xticks_size, yticks_size)
    ax1a = GeoAxis(
            ax;
            dest="+proj=eqc",
            title=title_name,
            lonlims=lonlims,
            latlims=latlims, 
            ylabel = ylabel_name,
            yticks = [29, 30, 31],
            width = width,
            height = height,
            xticklabelsize = xticks_size,
            yticklabelsize = yticks_size,
            ylabelfont=:bold,
    )
    poly!(ax1a, states_geo; color=:white, strokewidth=1)
    h1 = CairoMakie.scatter!(ax1a, point_df[!, :lon], point_df[!, :lat]; color=var, colormap=colorscheme, colorrange=range, markersize=20, strokewidth=0, strokecolor=:grey17)
    if with_xticks == false
        hidexdecorations!(ax1a, grid = false)
    end
    if with_yticks == false
        hideydecorations!(ax1a, grid = false)
    end
    return h1
end

# when need to plot two variables
function point_axis_2(ax, title_name, lonlims, latlims, ylabel_name, width, height, states_geo, point_df, var, colorscheme, range, with_xticks, with_yticks, xticks_size, yticks_size, var2, point_df2)
    ax1a = GeoAxis(
            ax;
            dest="+proj=eqc",
            title=title_name,
            lonlims=lonlims,
            latlims=latlims, 
            ylabel = ylabel_name,
            yticks = [29, 30, 31],
            width = width,
            height = height,
            xticklabelsize = xticks_size,
            yticklabelsize = yticks_size,
            ylabelfont=:bold,
    )
    poly!(ax1a, states_geo; color=:white, strokewidth=1)
    h1 = CairoMakie.scatter!(ax1a, point_df[!, :lon], point_df[!, :lat]; color=var, colormap=colorscheme, colorrange=range, markersize=20, strokewidth=0)
    CairoMakie.scatter!(ax1a, point_df2[!, :lon], point_df2[!, :lat]; color=var2, colormap=colorscheme, colorrange=range, markersize=20, strokewidth=2, strokecolor=:red)
    if with_xticks == false
        hidexdecorations!(ax1a, grid = false)
    end
    if with_yticks == false
        hideydecorations!(ax1a, grid = false)
    end
    return h1
end


"""plot time series of return level estimates"""

function plot_time_series(years, city_title, dists_ts; credible_interval = 0.95, q = 0.99, lengend_loc = :topleft, x_label = "Year", y_label = "Return Level [inches/day]", Atlas14 = -99)

    p = Plots.plot(;
        xlabel = x_label,
        ylabel = y_label,
        legend = lengend_loc,
        title = city_title,
    )

    qup = 1 - (1 - credible_interval) / 2
    qlow = (1 - credible_interval) / 2

    ub = [quantile([quantile(d, q) for d in dists_ts[y]], qup) for y = 1:length(years)]
    lb = [quantile([quantile(d, q) for d in dists_ts[y]], qlow) for y = 1:length(years)]
    
    # historical time
    Plots.plot!(
        p,
        years[1:(length(years)-1)],
        ub[1:(length(years)-1)];
        fillbetween = lb[1:(length(years)-1)],
        fillcolor = :blue,
        fillalpha = 0.4,
        linecolor = false,
        label = "95 Credible Interval",
    )

    # future projection
    Plots.plot!(
        p,
        years[(length(years)-1):length(years)],
        ub[(length(years)-1):length(years)];
        fillbetween = lb[(length(years)-1):length(years)],
        fillcolor = :blue,
        fillalpha = 0.2,
        linecolor = false,
        label = false,
    )

    median = [quantile([quantile(d, q) for d in dists_ts[y]], 0.5) for y = 1:length(years)]
    
    # historical time
    Plots.plot!(
        p,
        years[1:(length(years)-1)],
        median[1:(length(years)-1)];
        color = :blue,
        label = "Median",
        linewidth = 2,
    )

    # future projection
    Plots.plot!(
        p,
        years[(length(years)-1):length(years)],
        median[(length(years)-1):length(years)];
        color = :blue,
        label = "RCP6",
        linewidth = 2,
        linestyle = :dash
    )

    if Atlas14 >= 0
        hline!(
            [Atlas14],
            linestyle = :dash,
            linewidth = 2,
            color = :red,
            label = "Atlas 14",
        )
    end

    return p
end

"""plot different return level with uncertainty bounds"""

# from James' codes
# https://github.com/jdossgollin/2022-elevation-robustness/blob/e5d6adb3dc3cffd9cd80ee91064b6decae483f7e/scripts/plotutils.jl#L50
function plot_return_period(
    Atlas14_rls,
    plot_title,
    y_range,
    gevs::Vector{<:Distributions.GeneralizedExtremeValue};
    color_scheme = ColorSchemes.algae,
    type = "Posterior",
    lengend_loc = :bottomright,
    x_label = "Return Period [years]",
    y_label = "Return Level [inches/day]",
    include_Atlas14 = true,
)
    rts = range(1.25, 500; length = 250) # return periods
    aeps = 1 .- 1 ./ rts # annual exceedance probability
    xticks = [2, 5, 10, 25, 50, 100, 250, 500]
    ranges = [0.95, 0.80, 0.5]

    p = Plots.plot(;
        xlabel = x_label,
        ylabel = y_label,
        xscale = :log,
        legend = lengend_loc,
        xticks = (xticks, string.(xticks)),
        ylim = y_range,
        title = plot_title,
    )

    for range in ranges
        qup = 1 - (1 - range) / 2
        qlow = (1 - range) / 2
        ub = [quantile([quantile(d, xi) for d in gevs], qup) for xi in aeps]
        lb = [quantile([quantile(d, xi) for d in gevs], qlow) for xi in aeps]
        range_pct = Int(range * 100)
        fillcolor = ColorSchemes.get(color_scheme, range - 0.5)
        Plots.plot!(
            p,
            rts,
            ub;
            fillbetween = lb,
            fillcolor = fillcolor,
            fillalpha = 1,
            linecolor = false,
            label = "$(range_pct)% Credible Interval",
        )
    end

    median = [quantile([quantile(d, xi) for d in gevs], 0.50) for xi in aeps]
    Plots.plot!(
        p,
        rts,
        median;
        color = ColorSchemes.get(color_scheme, 1.0),
        label = "$(type) Median",
        linewidth = 2,
    )

    if include_Atlas14 == true
        Plots.plot!(
            p,
            [1, 2, 5, 10, 25, 50, 100, 200, 500],
            Atlas14_rls;
            color = :orangered2,
            label = "Altas 14",
            linewidth = 3.5,
        )
    end

    return p
end

"""plot gridded estimates"""

function map_grids(
    point_df,
    lons_transformed,
    lats_transformed,
    var_transformed,
    label_name,
    range,
    colorscheme,
    title_name,
    res,
)

    states = download(
        "https://raw.githubusercontent.com/PublicaMundi/MappingAPI/master/data/geojson/us-states.json",
    )
    states_geo = GeoJSON.read(read(states, String))

    Δ = 0.25
    lonlims = (minimum(point_df[!, :lon]) - Δ, maximum(point_df[!, :lon]) + Δ)
    latlims = (minimum(point_df[!, :lat]) - Δ, maximum(point_df[!, :lat]) + Δ)

    f = Figure(resolution = res)
    ax = GeoAxis(
        f[1, 1];
        dest = "+proj=eqc",
        title = title_name,
        lonlims = lonlims,
        latlims = latlims,
    )

    h = CairoMakie.heatmap!(
        ax,
        lons_transformed,
        lats_transformed,
        var_transformed,
        colormap = colorscheme,
        colorrange = range,
    )

    Colorbar(f[1, 2], h, width = 10, label = label_name, height = Relative(0.85))
    poly!(ax, states_geo, edgecolor = :gray, strokewidth = 1, color = :transparent)

    return f
end

function map_grids_subplots(
    point_df,
    lons_transformed,
    lats_transformed,
    var_transformed_all,
    colorschemes,
    rows,
    cols,
    row_names,
    res,
    row_hs,
    title_names,
    ranges;
    diff_coord = false,
    bar_all = true,
    diff_colbar = false,
    bar_name = false,
)

    states = download(
        "https://raw.githubusercontent.com/PublicaMundi/MappingAPI/master/data/geojson/us-states.json",
    )
    states_geo = GeoJSON.read(read(states, String))

    f = Figure(resolution = res, fontsize = 20, layout = GridLayout(row_heights = row_hs))

    Δ = 0.25
    lonlims = (minimum(point_df[!, :lon]) - Δ, maximum(point_df[!, :lon]) + Δ)
    latlims = (minimum(point_df[!, :lat]) - Δ, maximum(point_df[!, :lat]) + Δ)

    for i = 1:length(var_transformed_all)
        if diff_colbar == true
            colorscheme = colorschemes[i]
        else
            colorscheme = colorschemes
        end

        ga = GeoAxis(
            f[rows[i], cols[i]];
            dest = "+proj=eqc",
            title = title_names[i],
            lonlims = lonlims,
            latlims = latlims,
            ylabel = row_names[i],
        )

        h = CairoMakie.heatmap!(
            ga,
            lons_transformed,
            lats_transformed,
            var_transformed_all[i],
            colormap = colorscheme,
            colorrange = ranges[i],
            markersize = 20,
        )

        if bar_name == false
            b_name = false
        else
            b_name = bar_name[i]
        end

        if bar_all == true
            Colorbar(f[rows[i], cols[i]+1], h, label = b_name, height = Relative(0.9))
        else
            if cols[i] == maximum(cols)
                Colorbar(f[rows[i], cols[i]+1], h, label = b_name, height = Relative(0.9))
            end
        end

        poly!(ga, states_geo, edgecolor = :gray, strokewidth = 1, color = :transparent)
    end
    return f
end