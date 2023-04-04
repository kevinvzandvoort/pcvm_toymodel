prepareDataForPosterior = function(res_prevacc, data, stratify_cluster = FALSE){
  bycols = c("run", "compartment", "name.y")
  if(stratify_cluster) bycols = c(bycols, "cluster")
  
  age_groups_model %>% matchingAgeBreaks(data) %>% .[, age := 1:.N] %>% merge(res_prevacc, by="age") %>%
    .[, .(value = sum(value), N = sum(N), prev = sum(value)/sum(N)), by=bycols] %>%
    (function(x){
      if(stratify_cluster) x %>% dcast(run + cluster + name.y ~ compartment, value.var = "prev")
      else x %>% dcast(run + name.y ~ compartment, value.var = "prev")
    }) %>%
    #Count B in both VT and NVT
    #.[, VT := VT + B*1] %>% #count B in VT, assume it is the dominant serotype
    #.[, NVT := NVT + B*0] %>% #count B in NVT, assume it is the dominant serotype
    #.[, -"B"] %>%
    melt(measure.vars = c("S", "VT", "NVT", "B"), variable.name = "compartment", value.name = "modelled")
}

prepareDataForLL = function(res_prevacc, data){
  age_groups_model %>% matchingAgeBreaks(data) %>% .[, age := 1:.N] %>% merge(res_prevacc, by="age") %>%
    .[, .(value = sum(value), N = sum(N), prev = sum(value)/sum(N)), by=c("compartment", "name.y")] %>%
    dcast(name.y ~ compartment, value.var = "prev") %>%
    #Count B in both VT and NVT
    .[, VT := VT + B*1] %>% #count B in VT, assume it is the dominant serotype
    .[, NVT := NVT + B*0] %>% #count B in NVT, assume it is the dominant serotype
    .[, -"B"] %>%
    melt(measure.vars = c("S", "VT", "NVT"), variable.name = "compartment", value.name = "modelled")
}

#' Calculate the log likelihood based on the pre-vaccination model output (only the last timestep)
calculateLLprevacc = function(res_prevacc, data, return_data = FALSE){
  res_prevacc %>% prepareDataForLL(data) %>%
    merge(data %>% .[, -c("from", "to")] %>% melt(id.vars = "name", variable.name = "compartment", value.name = "observed"),
          by.x = c("name.y", "compartment"), by.y = c("name", "compartment"), all = TRUE) %>%
    .[, age := "name.y"] %>% .[, -"name.y"] %>% multinomial_log_ll() %>% sum(na.rm = TRUE)
}

#' Function that will be used in MCMC algorithm
calculateLL = function(model_params){
  
  #' Run model pre-vaccination
  res_prevacc = lsoda(
    y=model_params$state_prop_cluster, times=model_params$times_out_prevacc, func = "derivs", parms = model_params$params_unvac,
    dllname = "pcvm", nout = 1, initfunc = "initmod", outnames = "output") %>% tryCatchWE
  
  #' Only continue if no errors or warning, otherwise reject sample by returning -Inf
  if(res_prevacc$status == 0){
    
    #' Get states at steady state
    res_prevacc = res_prevacc$value %>% as.data.table %>% .[, -"output"] %>% .[time == max(time)] %>% 
      reshapeModelOutput(model_params$params_unvac)
    save_data = res_prevacc
    
    #' Get population size by age in all populations to calculate total prevalence for weighted mean
    N = model_params$params_unvac$trial_arms %>% seq_along %>%
      lapply(function(x, clusters) data.table(N = clusters[[x]]$parameters$N, cluster = names(clusters)[x]) %>% .[, age := 1:.N],
             model_params$params_unvac$trial_arms) %>% rbindlist
    res_prevacc = res_prevacc %>% merge(N, by = c("cluster", "age")) %>% .[, value := value * N]
    
    #' Calculate log likelihood only for the idp population
    log_ll_total = res_prevacc %>% calculateLLprevacc(data_unvacll)
  } else {
    log_ll_total = -Inf
  }
  
  if(is.infinite(log_ll_total)) save_data = NULL
  return(list(log_ll_total = log_ll_total, save_data = save_data))
}

lapplyNamed = function(x, fun){
  if(is.null(names(x))) stop("Element has no names")
  if(length(formals(fun)) != 1) stop("Function does not only have a single argument")
  if(is.list(x)){
    lapply(names(x), function(x, vals = x){ x = setNames(list(vals[[x]]), x); fun(x) }, x)
  } else {
    lapply(names(x), function(x, vals = x){ x = vals[x]; fun(x) }, x) 
  }
}

#' Calculate the multinomial log likelihood  
multinomial_log_ll = function(data){
  if(nrow(data[modelled <= 0]) > 0){
    warning("At least one modelled value was <0, rejecting the proposal.")
    return(-Inf)
  } else
    data[, LL := observed * log(modelled)] %>%
    .[, .(LL = sum(LL)), by = c("age")] %>%
    .[, LL] %>%
    return
}