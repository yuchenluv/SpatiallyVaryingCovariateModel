"""
Common functions for calculations
"""


"""functions to read MCMC csv files saved from r"""

# pooled stationary model
function read_MCMC_results_S(file_path, n_stations)
    stationary_pooled_posterior = DataFrame(CSV.File(datadir(file_path)))
    select!(stationary_pooled_posterior, Not(:Column1))

    mu_rho = stationary_pooled_posterior[:, 1]
    mu_alpha = stationary_pooled_posterior[:, 2]
    logs_rho = stationary_pooled_posterior[:, 3]
    logs_alpha = stationary_pooled_posterior[:, 4]
    xi = stationary_pooled_posterior[:, 5]
    mu = stationary_pooled_posterior[:, 6:5+n_stations]
    logs = stationary_pooled_posterior[:, 6+n_stations:5+n_stations*2]

    return mu_rho, mu_alpha, logs_rho, logs_alpha, xi, mu, logs
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


"""Functions to get the simulated distributions for points"""

# Pooled Stationary Model
function get_dist_MCMC_S(station_i, mu, logs, xi)
    mu_k = mu[:, station_i]
    logs_k = logs[:, station_i]
    return GeneralizedExtremeValue.(mu_k, exp.(logs_k), xi)
end

# Nonpooled Nonstationary Model
function get_dist_MCMC_N(station_i, μ_beta, logσ_beta, μ0_posterior, logσ0_posterior, ξ_posterior, x)
    μ_beta_k = μ_beta[:, station_i]
    logσ_beta_k = logσ_beta[:, station_i]
    μ0_k = μ0_posterior[:, station_i]
    logσ0_k = logσ0_posterior[:, station_i]
    ξ_k = ξ_posterior[:, station_i]
    μ_k = μ0_k .+ x .* μ_beta_k
    σ_k = exp.(logσ0_k .+ x .* logσ_beta_k)
    return GeneralizedExtremeValue.(μ_k, σ_k, ξ_k)
end

# Spatially Varying Covariate Model
function get_dist_MCMC_full(station_i, μ_beta, logσ_beta, μ0_posterior, logσ0_posterior, ξ_posterior, x)
    μ_beta_k = μ_beta[:, station_i]
    logσ_beta_k = logσ_beta[:, station_i]
    μ0_k = μ0_posterior[:, station_i]
    logσ0_k = logσ0_posterior[:, station_i]
    μ_k = μ0_k .+ x .* μ_beta_k
    σ_k = exp.(logσ0_k .+ x .* logσ_beta_k)
    return GeneralizedExtremeValue.(μ_k, σ_k, ξ_posterior)
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

function sim_quantiles(prcp, μ_beta, logσ_beta, μ0, logσ0, ξ, logCO2, total_records, n_sims)
    i = 0
    sim_qs = zeros(n_sims * total_records)
    for y in 1:size(prcp)[1]
        x = logCO2[y]
        for s in 1:size(prcp)[2]
            if ismissing(prcp[y, s]) || prcp[y, s] < 0
                nothing
            else
                i = i + 1
                sim_q_s = sim_quantile(prcp[y, s] / 25.4, μ_beta[:, s], logσ_beta[:, s], μ0[:, s], logσ0[:, s], ξ, x)
                sim_qs[n_sims*(i-1)+1:n_sims*i] = sim_q_s
            end
        end
    end
    return sim_qs
end


"""GP interpolation with MVN"""

using GaussianProcesses

function GP_param(rho, alpha, var_sim, x_old, lons_new, lats_new)
    # alpha here is actually estimated as std
    gp = GPE(x_old, var_sim, MeanZero(), Mat12Iso(log(rho), log(alpha)), 1e-6)
    mean_new, cov_new = predict_y(gp, hcat(lats_new, lons_new)')
    return mean_new
end

# interprete results from one simulation
function GP_dist(mu_rho, mu_alpha, logs_rho, logs_alpha, mu0_rho, mu0_alpha, logs0_rho, logs0_alpha, xi, mu_beta, logs_beta, mu0, logs0, x_old, lons_new, lats_new, x)
    mu0_new = GP_param(mu0_rho, mu0_alpha, mu0, x_old, lons_new, lats_new)
    mu_beta_new = GP_param(mu_rho, mu_alpha, mu_beta, x_old, lons_new, lats_new)
    logs0_new = GP_param(logs0_rho, logs0_alpha, logs0, x_old, lons_new, lats_new)
    logs_beta_new = GP_param(logs_rho, logs_alpha, logs_beta, x_old, lons_new, lats_new)
    return get_gev.(x, mu0_new, mu_beta_new, logs0_new, logs_beta_new, xi)
end


"""web scrape Atlas 14 estimates"""

# Function to extract and parse the quantiles from the fetched data
function parse_quantiles_from_url(url)
    # Fetching data from the URL
    response = HTTP.get(url, require_ssl_verification=false)
    data = String(response.body)

    # Extracting the quantiles data
    start_idx = findfirst("quantiles = ", data)
    end_idx = findnext(";", start_idx) - 1
    quantiles_str = data[start_idx:end_idx]

    # Cleaning up the quantiles string for JSON parsing
    quantiles_str = replace(quantiles_str, "quantiles = " => "")
    quantiles_str = replace(quantiles_str, "'" => "\"")

    # Parsing the JSON string to a Julia array
    quantiles_array = JSON.parse(quantiles_str)

    # Converting the string values to floats
    quantiles_floats = map(x -> parse(Float64, x), quantiles_array)

    return quantiles_floats
end

function get_Atlas14(lat, lon, row, col)
    # row is linked with durations
    url = "https://hdsc.nws.noaa.gov/cgi-bin/hdsc/new/cgi_readH5.py?lat=" * lat * "&lon=" * lon * "&type=pf&data=depth&units=english&series=pds"
    a = DataFrame(CSV.File(download(url)))[1, 1]
    a = a[15:length(a)-2]
    a = replace(a, "'" => "")
    b = split(a, "], [")
    c1 = [parse(Float64, s) for s in split(b[row], ", ")][col]
    return c1
end

function get_Atlas14_IDF(lat, lon)
    url = "https://hdsc.nws.noaa.gov/cgi-bin/hdsc/new/cgi_readH5.py?lat=" * lat * "&lon=" * lon * "&type=pf&data=depth&units=english&series=pds"
    a = DataFrame(CSV.File(download(url)))[1, 1]
    a = a[15:length(a)-2]
    a = replace(a, "'" => "")
    b = split(a, "], [")
    rls = [[parse(Float64, s) for s in split(b[10], ", ")][i] for i in [1, 2, 3, 4, 5, 6, 7, 8, 9]]
    return rls
end

