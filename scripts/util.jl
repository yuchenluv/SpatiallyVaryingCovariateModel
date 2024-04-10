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

"""functions to read MCMC csv files saved from r"""

# pooled stationary model
function read_MCMC_results_S(file_path, n_stations)
    stationary_pooled_posterior = DataFrame(CSV.File(datadir(file_path)))
    select!(stationary_pooled_posterior, Not(:Column1))

    # mu_w = stationary_pooled_posterior[:, 1]
    # logs_w = stationary_pooled_posterior[:, 2]
    # rho = stationary_pooled_posterior[:, 3]
    # alpha = stationary_pooled_posterior[:, 4]
    # xi = stationary_pooled_posterior[:, 5]
    # mu = stationary_pooled_posterior[:, 6:5+n_stations]
    # logs = stationary_pooled_posterior[:, 6+n_stations:5+n_stations*2]
    mu_w = stationary_pooled_posterior[:, 1]
    logs_w = stationary_pooled_posterior[:, 2]
    rho = stationary_pooled_posterior[:, 3]
    alpha = stationary_pooled_posterior[:, 4]
    mu = stationary_pooled_posterior[:, 5:4+n_stations]
    logs = stationary_pooled_posterior[:, 5+n_stations:4+n_stations*2]
    xi = stationary_pooled_posterior[:, 5+n_stations*2]

    return mu_w, logs_w, rho, alpha, mu, logs, xi
end

# nonpooled nonstationary model
function read_MCMC_results_N(file_path, n_stations)
    nonstationary_nonpooled_posterior = DataFrame(CSV.File(datadir(file_path)))
    select!(nonstationary_nonpooled_posterior, Not(:Column1))

    mu0 = nonstationary_nonpooled_posterior[:, 1:n_stations]
    logs0 = nonstationary_nonpooled_posterior[:, n_stations+1:n_stations*2]
    xi = nonstationary_nonpooled_posterior[:, n_stations*2+1:n_stations*3]
    mu_beta = nonstationary_nonpooled_posterior[:, n_stations*3+1:n_stations*4]
    logs_beta = nonstationary_nonpooled_posterior[:, n_stations*4+1:n_stations*5]
    return mu0, logs0, xi, mu_beta, logs_beta
end

# spatially varying covariate model
# "rho", "alpha", "mu_w", "logs_w", "mu0_w", "logs0_w", "xi"
# same rho for all spatially varying parameters
function read_MCMC_results_full_shareRho1(file_path, n_stations)
    nonstationary_multiGP_posterior = DataFrame(CSV.File(datadir(file_path)))
    select!(nonstationary_multiGP_posterior, Not(:Column1))
    rho = nonstationary_multiGP_posterior[:, 1]
    alpha = nonstationary_multiGP_posterior[:, 2]
    mu_w = nonstationary_multiGP_posterior[:, 3]
    logs_w = nonstationary_multiGP_posterior[:, 4]
    μ0_w = nonstationary_multiGP_posterior[:, 5]
    logσ0_w = nonstationary_multiGP_posterior[:, 6]
    xi = nonstationary_multiGP_posterior[:, 7]

    μ_beta = nonstationary_multiGP_posterior[:, 8:7+n_stations]
    logσ_beta = nonstationary_multiGP_posterior[:, 8+n_stations:7+n_stations*2]
    μ0 = nonstationary_multiGP_posterior[:, 8+n_stations*2:7+n_stations*3]
    logσ0 = nonstationary_multiGP_posterior[:, 8+n_stations*3:7+n_stations*4]
    return rho, alpha, mu_w, logs_w, μ0_w, logσ0_w, xi, μ_beta, logσ_beta, μ0, logσ0
end

# same rho for mu0 and mu_beta; logs0 and logs_beta
function read_MCMC_results_full_shareRho2(file_path, n_stations)
    nonstationary_multiGP_posterior = DataFrame(CSV.File(datadir(file_path)))
    select!(nonstationary_multiGP_posterior, Not(:Column1))
    mu_rho = nonstationary_multiGP_posterior[:, 1]
    mu_alpha = nonstationary_multiGP_posterior[:, 2]
    logs_rho = nonstationary_multiGP_posterior[:, 3]
    logs_alpha = nonstationary_multiGP_posterior[:, 4]
    mu_w = nonstationary_multiGP_posterior[:, 5]
    mu0_w = nonstationary_multiGP_posterior[:, 6]
    logs_w = nonstationary_multiGP_posterior[:, 7]
    logs0_w = nonstationary_multiGP_posterior[:, 8]
    xi = nonstationary_multiGP_posterior[:, 9]

    μ_beta = nonstationary_multiGP_posterior[:, 10:9+n_stations]
    logσ_beta = nonstationary_multiGP_posterior[:, 10+n_stations:9+n_stations*2]
    μ0 = nonstationary_multiGP_posterior[:, 10+n_stations*2:9+n_stations*3]
    logσ0 = nonstationary_multiGP_posterior[:, 10+n_stations*3:9+n_stations*4]
    return mu_rho, mu_alpha, logs_rho, logs_alpha, mu_w, mu0_w, logs_w, logs0_w, xi, μ_beta, logσ_beta, μ0, logσ0
end

# different rho for all spatially varying covariates
function read_MCMC_results_full_diffRho(file_path, n_stations)
    nonstationary_multiGP_posterior = DataFrame(CSV.File(datadir(file_path)))
    select!(nonstationary_multiGP_posterior, Not(:Column1))
    mu_rho = nonstationary_multiGP_posterior[:, 1]
    mu_alpha = nonstationary_multiGP_posterior[:, 2]
    logs_rho = nonstationary_multiGP_posterior[:, 3]
    logs_alpha = nonstationary_multiGP_posterior[:, 4]
    mu0_rho = nonstationary_multiGP_posterior[:, 5]
    mu0_alpha = nonstationary_multiGP_posterior[:, 6]
    logs0_rho = nonstationary_multiGP_posterior[:, 7]
    logs0_alpha = nonstationary_multiGP_posterior[:, 8]
    xi = nonstationary_multiGP_posterior[:, 9]

    μ_beta = nonstationary_multiGP_posterior[:, 10:9+n_stations]
    logσ_beta = nonstationary_multiGP_posterior[:, 10+n_stations:9+n_stations*2]
    μ0 = nonstationary_multiGP_posterior[:, 10+n_stations*2:9+n_stations*3]
    logσ0 = nonstationary_multiGP_posterior[:, 10+n_stations*3:9+n_stations*4]
    return mu_rho, mu_alpha, logs_rho, logs_alpha, mu0_rho, mu0_alpha, logs0_rho, logs0_alpha, xi, μ_beta, logσ_beta, μ0, logσ0
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
        rl_k = rl_estimate_N0.(p, x, μ0_k, μ_beta_k, logσ0_k, logσ_beta_k, ξ_posterior)
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
        rl_k = rl_estimate_N0.(p, x, μ0_k, μ_beta_k, logσ0_k, logσ_beta_k, ξ_k)
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

"""Functions for quantile plots"""

# get the simulated distribution for the full model
function get_gev(x, μ0, μ_beta, logσ0, logσ_beta, ξ)
    μ = μ0 + x * μ_beta
    σ = exp(logσ0 + x * logσ_beta)
    return GeneralizedExtremeValue(μ, σ, ξ)
end

# quantile of one observation given the simulated distributions
function sim_quantile(prcp_s, μ_beta_s, logσ_beta_s, μ0_s, logσ0_s, ξ_s, x)
    sim_dist = get_gev.(x, μ0_s, μ_beta_s, logσ0_s, logσ_beta_s, ξ_s)
    sim_q = cdf.(sim_dist, prcp_s)
    return sim_q
end

function sim_quantiles(prcp, μ_beta, logσ_beta, μ0, logσ0, ξ, logCO2, total_records)
    i = 0
    sim_qs = zeros(10000 * total_records)
    for y in 1:size(prcp)[1]
        x = logCO2[y]
        for s in 1:size(prcp)[2]
            if ismissing(prcp[y, s]) || prcp[y, s] < 0
                nothing
            else
                i = i + 1
                sim_q_s = sim_quantile(prcp[y, s] / 25.4, μ_beta[:, s], logσ_beta[:, s], μ0[:, s], logσ0[:, s], ξ, x)
                sim_qs[10000*(i-1)+1:10000*i] = sim_q_s
            end
        end
    end
    return sim_qs
end

"""GP interpolation with MVN"""

using GaussianProcesses

function GP_param(rho, alpha, w, var_sim, x_old, lons_new, lats_new)
    gp = GPE(x_old, var_sim, MeanZero(), Mat12Iso(log(rho), log(sqrt(alpha^2 * w))), 1e-6)
    mean_new, cov_new = predict_y(gp, hcat(lons_new, lats_new)')
    return mean_new
end

# interprete results from one simulation
function GP_dist(rho, alpha, mu_w, logs_w, mu0_w, logs0_w, xi, mu_beta, logs_beta, mu0, logs0, x_old, lons_new, lats_new, x)
    mu0_new = GP_param(rho, alpha, mu0_w, mu0, x_old, lons_new, lats_new)
    mu_beta_new = GP_param(rho, alpha, mu_w, mu_beta, x_old, lons_new, lats_new)
    logs0_new = GP_param(rho, alpha, logs0_w, logs0, x_old, lons_new, lats_new)
    logs_beta_new = GP_param(rho, alpha, logs_w, logs_beta, x_old, lons_new, lats_new)
    return get_gev.(x, mu0_new, mu_beta_new, logs0_new, logs_beta_new, xi)
end

