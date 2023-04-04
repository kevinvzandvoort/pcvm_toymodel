#' Age groups that will be used in the model
#' - good to use granular age groups in the early ages, to allow different vaccine strategies
age_groups_model = c(set_units(c(0:12), "month"), set_units(c(1:5, 6, 10, 15, 20, 30, 40, 50), "years")) %>%
  setAgeBreaks()

#' A - Generic constant values
MODEL_TIMESTEP = 1 #days. Nb model uses ODEs so this doesn't matter much, but all rates will be specified per MODEL_TIMESTEP
MODEL_YEARS_PREVACC = 4 #how long to run the model before vaccine introduction (to reach equilibrium)?
MODEL_YEARS_POSTVACC = 5 #how long to run the model after vaccine introduction?
NULL_RR = c(0) %>% setAgeBreaks() %>% .[, value := 1] %>%
  combineAgeBreaks(target = age_groups_model, additional = .) %>% .[, value] #Risk-ratio that will do nothing (TODO: refactor)
no_coverage = age_groups_model %>% getVaccineCoverage(0, 0) #used when coverage is 0 for all age groups

#' B - Infection specific parameters
#' clearance rates per day (divided by 7 as initial values provided per week)
clearance_rates = c(0, 2, 5) %>% setAgeBreaks() %>% .[, c("st", "value") := .("VT", c(0.062, 0.12, 0.34)/7)] %>%
  rbind(c(0, 2, 5) %>% setAgeBreaks() %>% .[, c("st", "value") := .("NVT", c(0.086, 0.15, 0.34)/7)])

#' The competition parameter reflects the extent to which carrying one or more pneumococcal serotypes of one group (carrying
#' VT or NVT) protects against acquisition of serotypes in the other group (NVT or VT). I.e. if competition is 20%, those
#' who carry VT serotypes experience 0.8 times the rate of infection with NVTs compared to someone who is susceptible
pneumo_competition = 0.9

#' C - Vaccine specific parameters
vaccine_coverage = 0.85
vaccine_efficacy_transmission = c(0) %>% setAgeBreaks() %>% .[, value := 0.5]
vaccine_efficacy_disease = c(0) %>% setAgeBreaks() %>% .[, value := 0.5] #combined efficacy of 75%
vaccine_waning_rates = c(0) %>% setAgeBreaks() %>% .[, value := 1/(8*365)] #on average, 8 years of vaccine protection

#' D - Population specific parameters
population_size_model = age_groups_model %>% combineAgeBreaks(hargeisa_population_data, method = "sum")
contact_matrix_model = adjustContactMatrixAgeGroups(age_groups_model, ethiopia_contact_matrix, hargeisa_contact_matrix_agegroups, population_size_model)

#' Specify populations to model
#' - Make sure to transpose the contact matrices so that columns denote contactees and rows contactors
#'   as this is required for the matrix multiplication when calculating the FOIs in the model
#'   i.e row i will show average number of contacts made by those of age i in all age groups, which will
#'   be multiplied with a column vector with prevalence in each contacted age groups and summed to give
#'   the average number of effective contacts with infectious individuals made by someone aged i (or j in
#'   contact matrix before transposing)
#' - TODO: need to refactor start and stop of acq adjustments. Set for negative and very large timesteps now
#'   to make sure they are always active
#' - This is where the majority of the model specification happens

model_populations = list(
  #Create one list for every population that you want to model (it could also just be a single population).
  "1+1" = list(
    "parameters" = list(
      #adjust_acq_vt and adjust_acq_nvt are age-specific risk-ratios that temporarily adjust the susceptibility towards
      # transmission, between timesteps adjust_acq_start and adjust_acq_stop. These are used in the Vietnam model to implement
      # the impact of NPIs during COVID-19. To disable, apply the NULL_RR vector and have them always be active
      adjust_acq_vt = NULL_RR, adjust_acq_nvt = NULL_RR,
      adjust_acq_start = -1, adjust_acq_stop = 1e6,
      #this is the contact matrix scaled by age-specific beta parameters for VTs (for lambda_VT calculation). Will be
      # overwritten when fitting the model
      betaVT = contact_matrix_model %>% t,
      #the same for NVTs (lambda_NVT)
      betaNVT = contact_matrix_model %>% t,
      #total population size in each age group, at the start of the model simulation. Nb the model simulates proportions in
      # each age group, NOT total numbers, but this allows calculation of IPD cases and appropriate travel rates
      N = population_size_model[, value]),
    #these are the different vaccinated arms in the population, one for each dose of vaccine protection
    #in addition to these, every population will ALWAYS have an unvaccinated arm
    #Nb arms is a bit poorly named - vaccine strata would be better
    "arms" = list(
        "1+0" = list(
          #coverage level, specified at the age people are vaccinated from the previous vaccine arm (here, the unvaccinated arm)
          "coverage" = c(0, 0, 0.68, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0),
          #catch-up coverage by age group, for a catch-up dose given at the start of the post-vaccination schedule 
          "catchup_coverage" = no_coverage,
          #waning rate by age, at which people wane to the previous vaccine arm (here, the unvaccinated arm)
          "waning" = age_groups_model %>% combineAgeBreaks(vaccine_waning_rates) %>% .[, value],
          #vaccine efficacy against transmission, by age, for people in this vaccinated arm
          "efficacy" = age_groups_model %>% combineAgeBreaks(vaccine_efficacy_transmission) %>% .[, value * 0.25]),
        "1+1" = list(
          "coverage" = c(0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0),
          "catchup_coverage" = age_groups_model %>% getVaccineCoverage(c(set_units(6, "weeks"), set_units(5, "years")), vaccine_coverage),
          "waning" = age_groups_model %>% combineAgeBreaks(vaccine_waning_rates) %>% .[, value],
          "efficacy" = age_groups_model %>% combineAgeBreaks(vaccine_efficacy_transmission) %>% .[, value]))),
  "2+1" = list(
    "parameters" = list(
      adjust_acq_vt = NULL_RR, adjust_acq_nvt = NULL_RR,
      adjust_acq_start = -1, adjust_acq_stop = 1e6,
      betaVT = contact_matrix_model %>% t,
      betaNVT = contact_matrix_model %>% t,
      N = population_size_model[, value]),
    "arms" = list(
      "1+0" = list(
        "coverage" = c(0, 0, 0.68, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0),
        "catchup_coverage" = no_coverage,
        "waning" = age_groups_model %>% combineAgeBreaks(vaccine_waning_rates) %>% .[, value],
        "efficacy" = age_groups_model %>% combineAgeBreaks(vaccine_efficacy_transmission) %>% .[, value * 0.25]),
      "2+0" = list(
        "coverage" = c(0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0),
        "catchup_coverage" = no_coverage,
        "waning" = age_groups_model %>% combineAgeBreaks(vaccine_waning_rates) %>% .[, value],
        "efficacy" = age_groups_model %>% combineAgeBreaks(vaccine_efficacy_transmission) %>% .[, value * 0.75]),
      "2+1" = list(
        "coverage" = c(0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0),
        "catchup_coverage" = age_groups_model %>% getVaccineCoverage(c(set_units(6, "weeks"), set_units(5, "years")), vaccine_coverage),
        "waning" = age_groups_model %>% combineAgeBreaks(vaccine_waning_rates) %>% .[, value],
        "efficacy" = age_groups_model %>% combineAgeBreaks(vaccine_efficacy_transmission) %>% .[, value]))),
  "3+0" = list(
    "parameters" = list(
      adjust_acq_vt = NULL_RR, adjust_acq_nvt = NULL_RR,
      adjust_acq_start = -1, adjust_acq_stop = 1e6,
      betaVT = contact_matrix_model %>% t,
      betaNVT = contact_matrix_model %>% t,
      N = population_size_model[, value]),
    "arms" = list(
      "1+0" = list(
        "coverage" = c(0, 0, 0.68, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0),
        "catchup_coverage" = no_coverage,
        "waning" = age_groups_model %>% combineAgeBreaks(vaccine_waning_rates) %>% .[, value],
        "efficacy" = age_groups_model %>% combineAgeBreaks(vaccine_efficacy_transmission) %>% .[, value * 0.25]),
      "2+0" = list(
        "coverage" = c(0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0),
        "catchup_coverage" = no_coverage,
        "waning" = age_groups_model %>% combineAgeBreaks(vaccine_waning_rates) %>% .[, value],
        "efficacy" = age_groups_model %>% combineAgeBreaks(vaccine_efficacy_transmission) %>% .[, value * 0.75]),
      "3+0" = list(
        "coverage" = c(0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0),
        "catchup_coverage" = age_groups_model %>% getVaccineCoverage(c(set_units(6, "weeks"), set_units(5, "years")), vaccine_coverage),
        "waning" = age_groups_model %>% combineAgeBreaks(vaccine_waning_rates) %>% .[, value],
        "efficacy" = age_groups_model %>% combineAgeBreaks(vaccine_efficacy_transmission) %>% .[, value]))))

#' E - Generic parameters for all populations used in the model
global_params = list(
  #total number of agegroups in the model
  n_agrp = age_groups_model[, .N],
  #competition parameter
  comp = pneumo_competition,
  #rate at which VTs and NVTs are cleared
  clearVT = age_groups_model %>% combineAgeBreaks(clearance_rates[st == "VT"]) %>% .[, value],
  clearNVT = age_groups_model %>% combineAgeBreaks(clearance_rates[st == "NVT"]) %>% .[, value],
  #rate at which individuals age between age groups
  ageout = age_groups_model %>% .[, .(duration = (to - from) %>% set_units("days"))] %>% .[, 1/as.numeric(duration)],
  #model populations
  trial_arms = model_populations,
  #travel matrix for contact between populations. travel[i, j] is proportion of contacts of someone in population j, that are made with those in
  # population i. All columns should sum to 1. Ideally only specified to one half of the matrix, where the other half is given the value -1 (will be automatically calculated)
  # to maintain reciprocity of contacts
  travel = diag(rep(1, 3)),
  #migration matrix for contact between populations. migration[i, j] is the rate at which someone from population j migrates to population i.
  migration = diag(rep(0, 3)))

#' F - Setup model parameters passed to the model
p_unvac = global_params
p_unvac$trial_arms %<>% lapply(function(x){x$arms = list(); return(x)}) #there are no vaccinated strata in the before vaccination
model_params = list(
  #parameters to use before vaccination
  params_unvac = p_unvac,
  #parameters to use after vaccination
  params_vac = global_params,
  #timesteps to output data before vaccination 
  times_out_prevacc = (set_units(MODEL_YEARS_PREVACC, "years") %>% set_units(days) %>% as.numeric() / MODEL_TIMESTEP) %>% ceiling() %>% (function(x) seq(0, x)),
  #timesteps to output data after vaccination
  times_postvacc_eval = (set_units(MODEL_YEARS_POSTVACC, "years") %>% set_units(days) %>% as.numeric() / MODEL_TIMESTEP) %>% ceiling() %>% (function(x) seq(0, x)),
  #initial stage of each cluster before vaccination
  state_prop_cluster = c(
    rep(0.8, age_groups_model[, .N]), #S
    rep(0.1, age_groups_model[, .N]), #VT
    rep(0.1, age_groups_model[, .N]), #NVT
    rep(0, age_groups_model[, .N]) #B
  ) %>% rep(length(model_populations)))

#' Adjust for timestep
model_params = model_params %>% adjustForTimeStep()
contact_matrix_model = contact_matrix_model %>% adjustForTimeStep()
