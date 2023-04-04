#' Required packages
pacman::p_load(data.table, Rcpp, RcppArmadillo, inline, deSolve, tmvtnorm, readxl, magrittr, binom, qs, units)

#' Create folder to store output
OUTPUT_FOLDER = "~/workspace/pcvm_toy/output/"
if(!dir.exists(OUTPUT_FOLDER)) dir.create(OUTPUT_FOLDER)

#' Load PCVm
PCVM_FOLDER = "~/workspace/pcvm/" #folder where pcvm is stored
PCVM_COMPILE = TRUE #should model be recompiled?
source(sprintf("%s/index.R", PCVM_FOLDER)) #build the model

#' Set up project specific functions and parameters
source("./functions.R")
source("./data_load_all.R")
source("./model_setup.R")

#' Parameters to fit in MCMC
priors = rbindlist(list(
  data.table(variable = "beta_VT_1", min = 0, max = 1, mean = 0.001, sd = 0.01,
             prior = function(x) dunif(x, min=0, max=1, log=TRUE)),
  data.table(variable = "beta_VT_2", min = 0, max = 1, mean = 0.001, sd = 0.01,
             prior = function(x) dunif(x, min=0, max=1, log=TRUE)),
  data.table(variable = "beta_VT_3", min = 0, max = 1, mean = 0.001, sd = 0.01,
             prior = function(x) dunif(x, min=0, max=1, log=TRUE)),
  data.table(variable = "beta_VT_NVT", min = 0, max = 2, mean = 1, sd = 0.1,
             prior = function(x) dunif(x, min=0, max=2, log=TRUE))))

#' Function that will be used to replace the model parameters with new values
#' - will be used when fitting parameters
updateParameters = function(params, model_params){
  #' Assign betaVT to the correct age groups
  betaVT = set_units(c(0, 5, 15), "year") %>% setAgeBreaks() %>%
    .[, value := c(params[variable == "beta_VT_1", value],
                   params[variable == "beta_VT_2", value],
                   params[variable == "beta_VT_3", value])] %>%
    combineAgeBreaks(target = age_groups_model, additional = .) %>% .[, value]
  
  betaNVT = betaVT/params[variable == "beta_VT_NVT", value]
  
  new_matrix_VT = sweep(contact_matrix_model, 2, betaVT, "*") %>% t
  new_matrix_NVT = sweep(contact_matrix_model, 2, betaNVT, "*") %>% t
  model_params[["params_unvac"]][["trial_arms"]][["1+1"]][["parameters"]][["betaVT"]] =
    new_matrix_VT
  model_params[["params_unvac"]][["trial_arms"]][["2+1"]][["parameters"]][["betaVT"]] =
    new_matrix_VT
  model_params[["params_unvac"]][["trial_arms"]][["3+0"]][["parameters"]][["betaVT"]] =
    new_matrix_VT
  
  model_params[["params_unvac"]][["trial_arms"]][["1+1"]][["parameters"]][["betaNVT"]] =
    new_matrix_NVT
  model_params[["params_unvac"]][["trial_arms"]][["2+1"]][["parameters"]][["betaNVT"]] =
    new_matrix_NVT
  model_params[["params_unvac"]][["trial_arms"]][["3+0"]][["parameters"]][["betaNVT"]] =
    new_matrix_NVT
  
  model_params[["params_vac"]][["trial_arms"]][["1+1"]][["parameters"]][["betaVT"]] =
    new_matrix_VT
  model_params[["params_vac"]][["trial_arms"]][["2+1"]][["parameters"]][["betaVT"]] =
    new_matrix_VT
  model_params[["params_vac"]][["trial_arms"]][["3+0"]][["parameters"]][["betaVT"]] =
    new_matrix_VT
  
  model_params[["params_vac"]][["trial_arms"]][["1+1"]][["parameters"]][["betaNVT"]] =
    new_matrix_NVT
  model_params[["params_vac"]][["trial_arms"]][["2+1"]][["parameters"]][["betaNVT"]] =
    new_matrix_NVT
  model_params[["params_vac"]][["trial_arms"]][["3+0"]][["parameters"]][["betaNVT"]] =
    new_matrix_NVT
  
  return(model_params)
}

#' Use this to test a single iteration of the model
model_params_current = priors[, c("variable", "mean")] %>%
  (function(x){z=x[, mean]; names(z)=x[, variable]; return(z)})
model_params = updateParameters(
  model_params_current %>% t %>% as.data.table %>% melt(id.vars=character(0)),
  model_params)
calculateLL(model_params)

#' Run MCMC algorithm
out = adaptiveMCMC(calculateLL, model_params, updateParameters, priors, mcmc_steps = 20000,
                   mcmc_adapt_size_start = 100, mcmc_adapt_shape_start = 500, mcmc_adapt_shape_stop = 1500,
                   output_file = sprintf("%s/mcmc_toymodel_fit.txt", OUTPUT_FOLDER))

#' Run diagnostics
x = coda::mcmc(out)
plot(x)
coda::autocorr.plot(x)