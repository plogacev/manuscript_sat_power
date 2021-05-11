
suppressPackageStartupMessages({
  library(plyr)
  library(dplyr)
  source("./simulations_functions.R")
  library(mrsat)
  library(doMC)
})

.S3method("logLik", "SATcurve", "logLik.SATcurve")
.S3method("summary", "SATcurve", "summary.SATcurve")

cat("using this many cores:", parallel::detectCores(logical=FALSE))
doMC::registerDoMC(parallel::detectCores(logical=FALSE))
#cat("using this many cores: 4")
#doMC::registerDoMC(4)

run_parallel = T
base_path <- "../simulation_output"
base_path_responses <- file.path(base_path, "1_responses")



fit_twocondition_models_by_participant <- function(fname_data)
{
  fname <- fname_data %>% gsub("/1_responses/", "/2_model_fits/", .)

  if (file.exists(fname) || !create_lockfile(fname)) {
    return (NULL)
  }

  cat("\n", basename(fname), "\n")
  
  responses_lst <- readRDS(fname_data)

  estimates <-
    plyr::ldply(seq_along(responses_lst), function(i)
    {
      cur_responses <- responses_lst[[i]]
      
      simulation_id <- unique(cur_responses$simulation_id)
      stopifnot(length(simulation_id) == 1)
      
      n_participants <- length(unique(cur_responses$participant_id))
      stopifnot( n_participants == max(cur_responses$participant_id) )
      
      pars_lst <- attributes(cur_responses)$pars
      pars <- plyr::ldply(pars_lst, function(par) par %>% .[-1:-2] %>% unlist() )
      pars %<>% dplyr::select(par_delta_asymptote.mean, par_delta_invrate.mean, par_delta_intercepts.mean, 
                              #par_delta_asymptote.sd, par_delta_invrate.sd, par_delta_intercepts.sd,
                              pop_asymptote1 = mean_asymptotes1, pop_asymptote2 = mean_asymptotes2, 
                              pop_invrate1 = mean_invrates1, pop_invrates2 = mean_invrates2,
                              pop_intercept1 = mean_intercepts1, pop_intercept2 = mean_intercepts2,
                              subj_intercept1 = intercepts1, subj_intercept2 = intercepts2,
                              subj_rate1 = rates1, subj_rates2 = rates2,
                              subj_asymptote1 = asymptotes1, subj_asymptote2 = asymptotes2
      )
      stopifnot( nrow(pars) == n_participants )
      
      #t <- system.time({
      estimates_i <- plyr::ddply(cur_responses, .(participant_id), function(cur_responses_ij) {
          j <- cur_responses_ij$participant_id[1]
          estimates_ij <- fit_models_two_conditions(cur_responses_ij)
          par_indices <- rep(j, nrow(estimates_ij))
          estimates_ij %<>% bind_cols( slice(pars, par_indices), .)
          estimates_ij
      })
      #})
      
      estimates_i %<>% cbind(simulation_id = rep(simulation_id, n_participants), .)
      
      estimates_i
    }, .parallel = run_parallel, .progress = "text")
  
  print(c("Finished", date()))
  
  # For some reason the number of simulations returned is not always the same as the number of simulations in the input,
  # and it doesn't seem to be systematic. Simulations go missing at random, and when code is re-run, usually either no simulations 
  # are missing, or different ones. There may be an interaction with the parallelization.
  # Make sure that results for all simulations are present.
  n_sims_match <- length(unique(estimates$simulation_id)) == length(responses_lst)
  if (n_sims_match) {
    saveRDS(estimates, file = fname)
    remove_lockfile(fname)
  } else {
    remove_lockfile(fname)
    stop("Simulations missing.")    
  }
}


# while(TRUE)
# {
    fnames_data <- dir(base_path_responses, pattern = "*.rds$", full.names = T, recursive = T)
    
    fnc <- fnames_data %>% 
            gsub(base_path_responses, "", .) %>% 
            gsub("sim_", "", .) %>% 
            gsub(".rds$", "", .) %>% 
            stringr::str_split_fixed("[-|_|/]", n=8) %>% 
            .[,-1]
    # tail(fnc)

    df_fnames <- data.frame(fnames = fnames_data,
                            response_type = fnc[,1],
                            n_subjects = fnc[,2] %>% as.double(),
                            n_responses_bysubj = fnc[,3] %>% as.double(),
                            #deltas_dist = fnc[,4],
                            delta_asym = fnc[,5] %>% as.double(),
                            delta_invrate = fnc[,6] %>% as.double(),
                            delta_incp = fnc[,7] %>% as.double()
                            )
    
    for (fname_data in sample(df_fnames$fnames)) {
      fit_twocondition_models_by_participant(fname_data = fname_data)
    }
# }
