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
# different rho (kernel length parameter) for all spatially varying covariates
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


"""functions to calculate correlation coefficients"""

function kendall_cor(idx, dataset, x_dataset)
    y = dataset[completecases(DataFrame(id=dataset[:, idx])), :]
    X = x_dataset[1:(length(x_dataset)-1)][completecases(DataFrame(id=dataset[:, idx])), :]
    # y = ustrip.(u"mm", y[:, idx])
    y = y[:, idx]
    y = y ./ 25.4
    cor = StatsBase.corkendall(X, y)
    return cor
end

function cor_df(dataset_obs, x_dataset)
    cor = [kendall_cor(idx, dataset_obs, x_dataset)[1] for idx in names(dataset_obs)]
    cor = DataFrame(stnid=names(dataset_obs), name=cor)
    return cor
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


"""Functions for quantile plots for the full model"""

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


"""Validation metrics"""

# LogS

# for stationary model
function logs_S(mu, logs, xi, obs_d)
    total_logs = 0.0
    for y in 1:size(obs_d)[1]
        for s in 1:size(obs_d)[2]
            if !ismissing(obs_d[y, s])
                # posterior probabilities for one observation given all MCMC posterior distributions
                pdf_yi = pdf.(GeneralizedExtremeValue.(mu[:, s], exp.(logs[:, s]), xi), obs_d[y, s])
                # negative log of the average among all MCMC posterior distributions
                logpdf_yi = -log(mean(pdf_yi))
                # sum up the negative log
                total_logs = total_logs + logpdf_yi
            end
        end
    end
    # average_logs = total_logs / (sum(!ismissing(x) for x in eachcol(obs_d) for x in x) * length(xi))
    return total_logs
end

# For nonstationary framework (Nonpooled Nonstationary & full model)
function logs_N(mu_beta, logs_beta, mu0, logs0, xi, CO2, obs_d, model_N)
    total_logs = 0.0
    for y in 1:size(obs_d)[1]
        for s in 1:size(obs_d)[2]
            if !ismissing(obs_d[y, s])
                mu = mu0[:, s] .+ CO2[y] .* mu_beta[:, s]
                sigma = exp.(logs0[:, s] .+ CO2[y] .* logs_beta[:, s])
                if model_N == "nonpooled"
                    xi_s = xi[:, s]
                else 
                    xi_s = xi
                end

                # posterior probabilities for one observation given all MCMC posterior distributions
                pdf_yi = pdf.(GeneralizedExtremeValue.(mu, sigma, xi_s), obs_d[y, s])
                # negative log of the average among all MCMC posterior distributions
                logpdf_yi = -log(mean(pdf_yi))
                # sum up the negative log
                total_logs = total_logs + logpdf_yi
            end
        end
    end
    return total_logs
end


# Quantile score

# mean for all simulations of one year and one station
function calculate_quantile_score0(obs, est_p, p)
    qs = zeros(Float64, length(est_p))
    for i in 1:length(est_p)
        if obs - est_p[i] > 0
            qs[i] = p * (obs - est_p[i])
        else
            qs[i] = (p - 1) * (obs - est_p[i])
        end
    end
    return mean(qs)
end

# for the stationary model
function calculate_quantile_score_S(obs, p, mu_s, logs_s, xi_s)
    # p is non-exceedance probability
    qs = Array{Union{Float64,Missing}}(zeros(size(obs)[1], size(obs)[2]))
    for t in 1:size(obs)[1]
        for s in 1:size(obs)[2]
            if ismissing(obs[t, s]) || obs[t, s] < 0
                qs[t, s] = missing
            else
                rl_p_yi = quantile.(GeneralizedExtremeValue.(mu_s[:, s], exp.(logs_s[:, s]), xi_s), p)
                qs[t, s] = calculate_quantile_score0(obs[t, s], rl_p_yi, p)
            end
        end
    end
    return sum(skipmissing(qs))
end

# for the nonstationary model
function calculate_quantile_score_N(obs, p, mu_beta, logs_beta, mu0, logs0, xi, CO2, model_N)
    # p is non-exceedance probability
    qs = Array{Union{Float64,Missing}}(zeros(size(obs)[1], size(obs)[2]))
    for t in 1:size(obs)[1]
        for s in 1:size(obs)[2]
            if ismissing(obs[t, s]) || obs[t, s] < 0
                qs[t, s] = missing
            else
                mu = mu0[:, s] .+ CO2[t] .* mu_beta[:, s]
                sigma = exp.(logs0[:, s] .+ CO2[t] .* logs_beta[:, s])
                if model_N == "nonpooled"
                    xi_s = xi[:, s]
                else 
                    xi_s = xi
                end

                rl_p_yi = quantile.(GeneralizedExtremeValue.(mu, sigma, xi_s), p)
                qs[t, s] = calculate_quantile_score0(obs[t, s], rl_p_yi, p)
            end
        end
    end
    return sum(skipmissing(qs))
end

# CRPS

# to calculate for one observation
using QuadGK  # A Julia package for numerical integration

function heaviside(y, observation)
    return y >= observation ? 1.0 : 0.0
end

function crps(observation, forecast_cdf)
    integrand(y) = (forecast_cdf(y) - heaviside(y, observation))^2
    result, _ = quadgk(integrand, obs_min, obs_max)
    return result
end

# for the stationary model
function cdf_S(μ, σ, ξ)
    function cdf_gev(y)
        return cdf(GeneralizedExtremeValue(μ, σ, ξ), y)
    end
    return cdf_gev
end

function crps_S(mu, logs, xi, obs_d)
    total_crps = 0.0
    for y in 1:size(obs_d)[1]
        print(y)
        print("  ")
        for s in 1:size(obs_d)[2]
            if !ismissing(obs_d[y, s])
                GEV_cdf_s = cdf_S.(mu[:, s], exp.(logs[:, s]), xi)
                crps_value = crps.(obs_d[y, s], GEV_cdf_s)
                total_crps += mean(crps_value)
            end
        end
    end
    return total_crps
end

# for nonstationary models (Full Model and Nonpooled Nonstationary Model)

function cdf_N(μ_beta, logσ_beta, μ0, logσ0, ξ, x)
    μ = μ0 + x * μ_beta
    σ = exp(logσ0 + x * logσ_beta)
    function cdf_gev(y)
        return cdf(GeneralizedExtremeValue(μ, σ, ξ), y)
    end
    return cdf_gev
end

function crps_N(mu_beta, logs_beta, mu0, logs0, xi, CO2, obs_d, model_N)
    total_crps = 0.0
    for y in 1:size(obs_d)[1]
        print(y)
        print("  ")
        for s in 1:size(obs_d)[2]
            if !ismissing(obs_d[y, s])
                if model_N == "nonpooled"
                    xi_s = xi[:, s]
                else 
                    xi_s = xi
                end

                GEV_cdf_y = cdf_N.(mu_beta[:, s], logs_beta[:, s], mu0[:, s], logs0[:, s], xi_s, CO2[y])
                crps_value = crps.(obs_d[y, s], GEV_cdf_y)
                total_crps += mean(crps_value)
            end
        end
    end
    return total_crps
end