pacman::p_load(data.table, Rcpp, RcppArmadillo, inline, deSolve, tmvtnorm, readxl, magrittr, binom, qs, units, ggplot2)

#' Create folder to store output
OUTPUT_FOLDER = "C:/workspace/pcvm_toy/output/"
OUTPUT_SUBFOLDER = "pcv_campaign_simulations"
if(!dir.exists(OUTPUT_FOLDER)) dir.create(OUTPUT_FOLDER)
if(!dir.exists(paste0(OUTPUT_FOLDER, "/", OUTPUT_SUBFOLDER))) dir.create(paste0(OUTPUT_FOLDER, OUTPUT_SUBFOLDER))

#' Load PCVm
PCVM_FOLDER = "C:/workspace/pcvm/"
PCVM_COMPILE = FALSE
source(sprintf("%s/index.R", PCVM_FOLDER))

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

#' Read the MCMC output file
output_file = "mcmc_toymodel_fit"
mcmc_settings = fread(sprintf("%s/mcmc_out_settings.csv", OUTPUT_FOLDER))
            
posterior = fread(sprintf("%s/mcmc_toymodel_fit", OUTPUT_FOLDER))
colnames(posterior)[1] = "iteration"
posterior = posterior[!is.na(loglikelihood)]
posterior[, stage := "init"]
            
if(!is.na(mcmc_settings[param == "mcmc_adapt_shape_start_at", value]))
  posterior[iteration >= mcmc_settings[param == "mcmc_adapt_shape_start_at", value], stage := "adapt_shape"]
if(!is.na(mcmc_settings[param == "mcmc_adapt_shape_stop_at", value]))
  posterior[iteration >= mcmc_settings[param == "mcmc_adapt_shape_stop_at", value], stage := "no_adapt"]
posterior[, stage := factor(stage, c("init", "adapt_size", "adapt_shape", "no_adapt"))]
posterior = melt(posterior, id.vars=c("iteration", "stage")) %>% setorder(iteration, variable) %>%
  .[as.numeric(stage) == max(as.numeric(stage))]
            
out_files = list.files(sprintf("%s/model_output/", OUTPUT_FOLDER))
out_iterations = gsub(".*?([0-9]+).*", "\\1", out_files) %>% as.numeric %>% sort
posterior_iterations = posterior$iteration %>% unique %>% sort %>% tail(10000)

i=1
for(i in 1:500){
  mcmc_it = sample(posterior_iterations, 1)
  
  #ccr_sampled = age_groups_model %>%
  #  sampleCaseCarrierRatio(case_carrier_data, case_carrier_data_agegroups)
  
  model_params = updateParameters(
    posterior[iteration == mcmc_it & variable != "loglikelihood"],
    model_params)
                
  res_prevacc = out_iterations[out_iterations <= mcmc_it] %>%
    rev %>% .[1] %>% (function(x) sprintf("%s/model_output/out_%s.qs", OUTPUT_FOLDER, x)) %>% qread
              
  postvac_out_carrier_cases = list()
  postvac_out_prevalence = list()
  
  state_prop_equil = res_prevacc %>% eqStatesVaccinate(model_params$params_vac)
                
  message("start model")
  
  res_postvacc = lsoda(y=state_prop_equil[, state], times=seq(0, 365*5, 1), func = "derivs", parms = model_params$params_vac,
                         dllname = "pcvm", nout = state_prop_equil[, .N], initfunc = "initmod",
                         outnames = paste0("incidence_", seq_len(state_prop_equil[, .N])))
    
  message("process model")
  res_postvacc = as.data.table(res_postvacc)
  res_postvacc_incidence = res_postvacc %>% reshapeModelOutput(model_params$params_vac)

  carrier_prevalence = res_postvacc_incidence %>%
    .[variable == "prevalence"] %>%
    merge(age_groups_model %>% copy %>% .[, age := 1:.N] %>% .[, c("age", "name")], by="age") %>% .[, -"age"] %>%
    .[, .(value = sum(value)), by=c("cluster", "compartment", "name", "time")] %>%
    dcast(cluster+name+time~compartment, value.var="value") %>% 
    .[, VT := VT+B] %>% .[, -"B"]
    
  carrier_prevalence[, vaccine_strategy := strategy]
  carrier_prevalence[, mcmciter := mcmc_it]
  carrier_prevalence[, pvac_iter := i]
    
  population_size_model = model_params$params_vac$trial_arms %>%
    lapplyNamed(function(x) data.table(cluster = names(x),
                                       N = x[[1]]$parameters$N) %>% .[, age := 1:.N]) %>%
    rbindlist() %>% .[, c("cluster", "age", "N")]
    
  carrier_cases = res_postvacc_incidence %>%
    .[variable == "incidence"] %>%
    merge(population_size_model, by=c("cluster", "age")) %>%
    merge(age_groups_model %>% copy %>% .[, age := 1:.N] %>% .[, c("age", "name")], by="age") %>%
    .[, -c("age")] %>% 
    .[, total := value*N] %>%
    .[, .(incidence = sum(total)), by=c("cluster", "dose", "compartment", "name", "time")] %>% 
    dcast(cluster+dose+name+time~compartment, value.var="incidence") %>% 
    .[, .(vt=S_VT+NVT_B, nvt=S_NVT+VT_B), by=c("cluster", "dose", "name", "time")]

    carrier_cases %<>% merge(ccr_sampled %>%
                           merge(age_groups_model %>% combineAgeBreaks(vaccine_efficacy_disease %>% copy %>% rbind(vaccine_efficacy_disease) %>%
                                                                         .[, c("dose", "value") := .(c("unvaccinated", "dose_optimal"), c(0, 0.75))],
                                                                       method = "mean", value.var = "value", additional_group = "dose") %>% .[, -c("from", "to")], by="name") %>%
                            .[, c("name", "NVT", "VT", "value", "dose")], by=c("name", "dose")) %>%
      .[, VT := VT * (1 - value)] %>% .[, -"value"] %>%
      .[, c("ipd_vt", "ipd_nvt") := .(vt*VT, nvt*NVT)] %>% .[, -c("NVT", "VT")] %>%
      .[, .(vt = sum(vt), nvt = sum(nvt), ipd_vt = sum(ipd_vt), ipd_nvt = sum(ipd_nvt)), by=c("cluster", "name", "time")]
    
    carrier_cases[, vaccine_strategy := strategy]
    carrier_cases[, mcmciter := mcmc_it]
    carrier_cases[, pvac_iter := i]
                    
    postvac_out_carrier_cases[[length(postvac_out_carrier_cases) + 1]] = carrier_cases
    postvac_out_prevalence[[length(postvac_out_prevalence) + 1]] = carrier_prevalence

    qs::qsave(postvac_out_prevalence, sprintf("%s/%s/iter_%s_prevalence.qs", OUTPUT_FOLDER, OUTPUT_SUBFOLDER, i))
    qs::qsave(postvac_out_carrier_cases, sprintf("%s/%s/iter_%s_carrier_cases.qs", OUTPUT_FOLDER, OUTPUT_SUBFOLDER, i))
}              
