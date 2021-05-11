
library(plyr)
library(dplyr)
library(magrittr)
library(ggplot2)
library(ggpubr)

run_parallel = T


fname2sim_type <- function(fnames)
{
  df <- data.frame(fname = fnames, 
                   base_fname = basename(fnames))
  df$sim_id <- df$base_fname %>% gsub(".rds", "", .) %>% gsub("sim-", "", .)

  fnc <- stringr::str_split_fixed(df$sim_id, "[-|_]", n=7)
  stopifnot(all(fnc[,4] == "vardelta"))

  df$response_type <- fnc[,1] %>% as.factor()
  df$n_subjects <- fnc[,2] %>% as.integer()
  df$n_responses_bysubj <- fnc[,3] # %>% as.integer()
  df$delta_asym <- fnc[,5] %>% as.double()
  df$delta_invrate <- fnc[,6] %>% as.double()
  df$delta_incp <- fnc[,7] %>% as.double()
  
  df
}

compute_t_tests <- function(estimates, use_dynamics = F)
{
  estimates$est_invrate1 <- 1/estimates$est_rate1
  estimates$est_invrate2 <- 1/estimates$est_rate2
  
  t_test_results <- estimates %>% 
          group_by(simulation_id, model) %>%
          dplyr::summarise( avg_delta_icpt21 = mean(est_intercept2-est_intercept1), 
                            t_icpt = t.test(est_intercept2-est_intercept1)$statistic[['t']],
                            
                            avg_delta_invrate21 = mean(est_invrate2-est_invrate1),
                            t_invrate = t.test(est_invrate2-est_invrate1)$statistic[['t']],
                            
                            avg_delta_asym21 = mean(est_asymptote2-est_asymptote1),
                            t_asym = t.test(est_asymptote2-est_asymptote1)$statistic[['t']],
                            
                            avg_delta_dynamics21 = mean( (est_intercept2+est_invrate2) - (est_intercept1+est_invrate1) ),
                            t_dynamics = t.test( (est_intercept2+est_invrate2) - (est_intercept1+est_invrate1) )$statistic[['t']],
                            
                            .groups = "drop"
                          )

  if (use_dynamics) {
    t_test_results %<>% dplyr::select(-avg_delta_invrate21, -t_invrate, -avg_delta_icpt21, -t_icpt)
    t_test_results %<>% subset( model %in% c("2-2-2", "2-1-1", "1-2-2", "1-1-1"))
    t_test_results$model %<>% gsub("-[12]$", "", .)
    
  } else {
    t_test_results %<>% dplyr::select(-avg_delta_dynamics21, -t_dynamics)
  }
  
  t_test_results 
}

select_stepwise <- function(df_t, df_adjR2 = NULL, forward = NA)
{
  stopifnot( !is.na(forward) )
  stopifnot( nrow(df_t) == 8 )

  use_adjR2 <- !is.null(df_adjR2)
  
  df_t_asym <- df_t %>% dplyr::select(model, t_asym) %>% tidyr::pivot_wider(names_from = "model", values_from = "t_asym")
  df_t_icpt <- df_t %>% dplyr::select(model, t_icpt) %>% tidyr::pivot_wider(names_from = "model", values_from = "t_icpt")
  df_t_invrate <- df_t %>% dplyr::select(model, t_invrate) %>% tidyr::pivot_wider(names_from = "model", values_from = "t_invrate")
  
  model_id <- function(n_asym, n_invrate, n_icpt) paste( c(n_asym, n_invrate, n_icpt), collapse = "-")
  
  percentage_m2_higher <- function(model_id1, model_id2) { if (use_adjR2) { df_adjR2[[paste(model_id1, model_id2, sep ="_")]] > 0 } else { T } }

  if (forward) { n_icpt = n_invrate = 1 } 
  else { n_icpt = n_invrate = 2 }

  model_id1 = model_id(1, n_invrate, n_icpt); model_id2 = model_id(2, n_invrate, n_icpt)
  consider_m2 <- percentage_m2_higher(model_id1, model_id2)
  t_asym = df_t_asym[[model_id2]]
  n_asym = ifelse( consider_m2 && abs(t_asym) > t_critical, 2, 1)
  model_delta_asym = map_t_val(t_asym * consider_m2)

  model_id1 = model_id(n_asym, n_invrate, 1); model_id2 = model_id(n_asym, n_invrate, 2)
  consider_m2 <- percentage_m2_higher(model_id1, model_id2)
  t_icpt = df_t_icpt[[model_id2]]
  n_icpt = ifelse(consider_m2 && abs(t_icpt) > t_critical, 2, 1)
  model_delta_icpt = map_t_val(t_icpt * consider_m2)

  model_id1 = model_id(n_asym, 1, n_icpt); model_id2 = model_id(n_asym, 2, n_icpt)
  consider_m2 <- percentage_m2_higher(model_id1, model_id2)
  t_invrate = df_t_invrate[[model_id2]]
  n_invrate = ifelse(consider_m2 && abs(t_invrate) > t_critical, 2, 1)
  model_delta_invrate = map_t_val(t_invrate * consider_m2)
  
  data.frame(model_delta_asym, model_delta_invrate, model_delta_icpt)
}

select_stepwise_dynamics <- function(df_t, forward = NA)
{
  stopifnot( !is.na(forward) )
  stopifnot( nrow(df_t) == 4 )

  # use_adjR2 <- !is.null(df_adjR2)
  
  df_t_asym <- df_t %>% dplyr::select(model, t_asym) %>% tidyr::pivot_wider(names_from = "model", values_from = "t_asym")
  df_t_dynamics <- df_t %>% dplyr::select(model, t_dynamics) %>% tidyr::pivot_wider(names_from = "model", values_from = "t_dynamics")
  
  model_id <- function(n_asym, n_dynamics) paste( c(n_asym, n_dynamics), collapse = "-")
  percentage_m2_higher <- function(model_id1, model_id2) { T }
  
  if (forward) { n_dynamics = 1 
  } else { n_dynamics = 2 }
  
  model_id1 = model_id(1, n_dynamics); model_id2 = model_id(2, n_dynamics)
  t_asym = df_t_asym[[model_id2]]
  n_asym = ifelse( abs(t_asym) > t_critical, 2, 1)
  model_delta_asym = map_t_val(t_asym )
  
  model_id1 = model_id(n_asym, 1); model_id2 = model_id(n_asym, 2)
  t_dynamics = df_t_dynamics[[model_id2]]
  n_dynamics = ifelse( abs(t_dynamics) > t_critical, 2, 1)
  model_delta_dynamics = map_t_val(t_dynamics)
  
  data.frame(model_delta_asym, model_delta_dynamics)
}

calculate_avg_delta_adjR2_by_simulation <- function(estimates)
{
  # compile existing ids for all simulations, participants, and models 
  uniq_simulation_id <- unique(estimates$simulation_id)
  uniq_participant_id <- unique(estimates$participant_id)
  uniq_model <- unique(estimates$model)
  
  # create an otherwise empty data frame with all participant ids and model1 and model2 combinations for all simulation ids
  model_comparison <- tidyr::expand_grid(simulation_id=uniq_simulation_id, participant_id=uniq_participant_id, model1=uniq_model, model2=uniq_model)

  # merge in model comparison statistics from model1 and model2
  join_by = c("simulation_id", "participant_id")
  model_comparison %<>% left_join( select(estimates,  simulation_id, participant_id, model1=model, adjR2_1=adjR2, logLik_1=logLik, Deviance_1=Deviance, AIC_1=AIC, BIC_1=BIC), by = c(join_by, "model1") ) # R2_1=R2, SSE_1=SSE, RMSEfit_1=RMSEfit, npar_1=npar,
  model_comparison %<>% left_join( select(estimates,  simulation_id, participant_id, model2=model, adjR2_2=adjR2, logLik_2=logLik, Deviance_2=Deviance, AIC_2=AIC, BIC_2=BIC), by = c(join_by, "model2") )

  # compute adjusted R^2 differency between each pair of models by simulation id, by participant 
  model_comparison_adjR2 <- model_comparison %>%
                dplyr::select(simulation_id, participant_id, model1, model2, adjR2_1, adjR2_2) %>% 
                dplyr::mutate(delta_adjR2 = adjR2_2 - adjR2_1) %>%
                dplyr::select(-adjR2_1, -adjR2_2)

  # aggregate adusted R^2 differences by averaging over participants 
  model_comparison_adjR2_by_simulation <-
    model_comparison_adjR2 %>% group_by(simulation_id,model1,model2) %>%
                dplyr::summarise( avg_delta_adjR2 = mean(delta_adjR2), .groups = "drop") %>% 
                tidyr::pivot_wider(names_from = c("model1","model2"), values_from = "avg_delta_adjR2")
  
  model_comparison_adjR2_by_simulation %<>% as.data.frame()
  rownames(model_comparison_adjR2_by_simulation) <- model_comparison_adjR2_by_simulation$simulation_id
  
  model_comparison_adjR2_by_simulation
}

# calculate_BIC_metrics <- function(estimates)
# {
#   View(estimates)
# 
#   t_tests <- compute_t_tests(estimates)
# 
#   estimates %>% select(simulation_id, model, logLik, BIC) %>% head(8)
# 
#   head(t_tests %>% subset(model == "2-2-2"))
# 
#   weights_BIC_sum <-
#     estimates %>% group_by(simulation_id, model) %>%
#     dplyr::summarise( logLik = sum(logLik), n = sum(n), npar = sum(npar) ) %>%
#     group_by(simulation_id) %>% dplyr::mutate(
#       BIC = n*log(npar) - 2*logLik,
#       delta_BIC = BIC - min(BIC),
#       weight_BIC = exp(-0.5*delta_BIC)/sum(exp(-0.5*delta_BIC)),
#       best = weight_BIC == max(weight_BIC)
#       )
# 
#   best_models <- weights_BIC_sum %>% subset(best == T)
# 
#   xtabs(~best_models$model)
# 
#   weights_BIC_sum <-
#       estimates %>% group_by(simulation_id, model) %>%
#       dplyr::summarise( BIC_sum = sum(BIC) ) %>%
#       group_by(simulation_id) %>% dplyr::mutate(
#                         delta_BIC_sum = BIC_sum - min(BIC_sum),
#                         weight_BIC_sum = exp(-0.5*delta_BIC_sum)/sum(exp(-0.5*delta_BIC_sum)),
#                         best_by_BIC_sum = weight_BIC_sum == max(weight_BIC_sum) )
# 
#   best_models <- weights_BIC_sum %>% subset(best_by_BIC_sum == T)
# 
#   xtabs(~best_models$model)
# 
# }

bootstrap_model_quality <- function(fname, n_participants, n_boot_samples, use_dynamics = F )
{ 
  # read in model fit statistics
  estimates <- readRDS( fname )
  
  # create participant ids uniqu to the data set
  estimates %<>% ungroup() %>% mutate(compound_participant_id = simulation_id * (max(participant_id)+1) + participant_id )
  compound_participant_ids <- unique(estimates$compound_participant_id)
    
  df_boot_samples <-
  plyr::ldply(1:n_boot_samples, function (i_boot_sample) {
    df_iter <- data.frame(i_boot_sample = i_boot_sample,
                          boot_iter_participant_id = 1:n_participants, 
                          compound_participant_id = sample(compound_participant_ids, size = n_participants, replace = T)) 
    df_iter %<>% left_join( estimates, by = "compound_participant_id" )
  })
  df_boot_samples %<>% dplyr::select(-simulation_id, -participant_id, -compound_participant_id)  
  df_boot_samples %<>% dplyr::rename(simulation_id=i_boot_sample, participant_id=boot_iter_participant_id)
  df_boot_samples %<>% ungroup()
  
  if (use_dynamics) {
      compute_model_quality_dynamics(df_boot_samples)
  } else {
      compute_model_quality(df_boot_samples)
  }
}

compute_model_quality <- function(estimates)
{  
    # prepare empty data frame for model selection results
    model_quality <- tidyr::expand_grid(simulation_id = unique(estimates$simulation_id), 
                                        model_delta_asym = c("+1","-1","0"),
                                        model_delta_invrate = c("+1","-1","0"),
                                        model_delta_icpt = c("+1","-1","0"))
    
    # Compute t-values by simulation id
    df_t <- compute_t_tests(estimates)
    
    # prepare information about model comparison based on average differences between adjusted R^2 values
    df_avg_delta_adjR2_by_simulation <- calculate_avg_delta_adjR2_by_simulation(estimates) 
    avg_delta_adjR2 <- function(simulation_id) df_avg_delta_adjR2_by_simulation[as.character(simulation_id),]
    
    ### Method 1: simple t-tests on a 2-2-2 model
    df_t_222 <- df_t %>% subset(model == "2-2-2")
    best_t_val222 <- df_t_222 %>%
                      group_by(simulation_id) %>%
                      dplyr::summarise(model_delta_asym = map_t_val(t_asym),
                                       model_delta_invrate = map_t_val(t_invrate),
                                       model_delta_icpt = map_t_val(t_icpt),
                                       .groups = "drop")
    model_quality %<>% dplyr::left_join( best_t_val222 %>% mutate(t_val222 = T) )
    model_quality$t_val222[is.na(model_quality$t_val222)] <- FALSE
    
    # ### Method 2: forward model selection, using t-values only
    best_forward_selection_t <- df_t %>%
              group_by(simulation_id) %>%
              dplyr::do( select_stepwise(., forward = T) )
    model_quality %<>% dplyr::left_join( best_forward_selection_t %>% mutate(best_forward_selection_t = T) )
    model_quality$best_forward_selection_t[is.na(model_quality$best_forward_selection_t)] <- FALSE
    
    ### Method 3: backward model selection, using t-values only
    best_backward_selection_t <- df_t %>%
              group_by(simulation_id) %>%
              dplyr::do( select_stepwise(., forward = F) )
    model_quality %<>% dplyr::left_join( best_backward_selection_t %>% mutate(best_backward_selection_t = T) )
    model_quality$best_backward_selection_t[is.na(model_quality$best_backward_selection_t)] <- FALSE


    ### Method 4: forward model selection, using adjR2 and t-values
    best_forward_selection_adjR2_t <- df_t %>%
              group_by(simulation_id) %>%
              dplyr::do( select_stepwise(df_t = ., df_adjR2 = avg_delta_adjR2(.$simulation_id[1]), forward = T) )
    model_quality %<>% dplyr::left_join( best_forward_selection_adjR2_t %>% mutate(best_forward_selection_adjR2_t = T) )
    model_quality$best_forward_selection_adjR2_t[is.na(model_quality$best_forward_selection_adjR2_t)] <- FALSE

    ### Method 5: backward model selection, using adjR2 and t-values
    best_backward_selection_adjR2_t <- df_t %>%
              group_by(simulation_id) %>%
              dplyr::do( select_stepwise(df_t = ., df_adjR2 = avg_delta_adjR2(.$simulation_id[1]), forward = F) )
    model_quality %<>% dplyr::left_join( best_backward_selection_adjR2_t %>% mutate(best_backward_selection_adjR2_t = T) )
    model_quality$best_backward_selection_adjR2_t[is.na(model_quality$best_backward_selection_adjR2_t)] <- FALSE

    model_quality %>% subset(t_val222 | best_backward_selection_t | best_forward_selection_t | best_forward_selection_adjR2_t | best_backward_selection_adjR2_t)
}

compute_model_quality_dynamics <- function(estimates)
{  
  # prepare empty data frame for model selection results
  model_quality <- tidyr::expand_grid(simulation_id = unique(estimates$simulation_id), 
                                      model_delta_asym = c("+1","-1","0"),
                                      model_delta_dynamics = c("+1","-1","0")
                                      )
  
  # Compute t-values by simulation id
  df_t <- compute_t_tests(estimates, use_dynamics = T)

  ### Method 1: simple t-tests on a 2-2(-2) model
  df_t_222 <- df_t %>% subset(model == "2-2")
  best_t_val222 <- df_t_222 %>%
    group_by(simulation_id) %>%
    dplyr::summarise(model_delta_asym = map_t_val(t_asym),
                     model_delta_dynamics = map_t_val(t_dynamics),
                     .groups = "drop")
  model_quality %<>% dplyr::left_join( best_t_val222 %>% mutate(t_val222 = T) )
  model_quality$t_val222[is.na(model_quality$t_val222)] <- FALSE
  
  # ### Method 2: forward model selection, using t-values only
  best_forward_selection_t <- df_t %>%
    group_by(simulation_id) %>%
    dplyr::do( select_stepwise_dynamics(., forward = T) )
  model_quality %<>% dplyr::left_join( best_forward_selection_t %>% mutate(best_forward_selection_t = T) )
  model_quality$best_forward_selection_t[is.na(model_quality$best_forward_selection_t)] <- FALSE
  
  ### Method 3: backward model selection, using t-values only
  best_backward_selection_t <- df_t %>%
    group_by(simulation_id) %>%
    dplyr::do( select_stepwise_dynamics(., forward = F) )
  model_quality %<>% dplyr::left_join( best_backward_selection_t %>% mutate(best_backward_selection_t = T) )
  model_quality$best_backward_selection_t[is.na(model_quality$best_backward_selection_t)] <- FALSE

  model_quality %>% subset(t_val222 | best_backward_selection_t | best_forward_selection_t )
}



library(parallel)
library(doMC)
doMC::registerDoMC(parallel::detectCores())



t_critical <- 1.96 # critical value of t
map_t_val <- function(t) { ifelse(t >= t_critical, "+1", ifelse(t <= -t_critical, "-1", "0")) }


# Note: model n-i-k notation stands for n asymptotes, i rates, k intercepts

dir_in <- "../simulation_output/2_model_fits/"
# dir_out <- "../simulation_output/3_model_comparison/"

fnames <- dir(dir_in, pattern = "*.rds$", full.names = T, recursive = T)
df_fnames <- fname2sim_type(fnames)


df_fnames %>% plyr::ddply(c("fname"), function(df_fname)
{
    print(df_fname)
    n_boot_samples = 1000
  
    for (n_subj in c(10, 20, 30, 40, 50)) # 
    {
        fname_out_vanilla <- df_fname$fname %>% gsub("2_model_fits", paste0("3_model_comparison_vanilla/n", n_subj), .)
        fname_out_dynamics <- df_fname$fname %>% gsub("2_model_fits", paste0("3_model_comparison_dynamics/n", n_subj), .)
        
        sim_info <- df_fname %>% dplyr::select(sim_id:delta_incp) %>% dplyr::rename(simulation_name = sim_id)
        sim_info$n_subjects <- n_subj
        sim_info$n_boot_samples <- n_boot_samples
        
        if (!file.exists(fname_out_dynamics)) 
        {
          model_quality <- bootstrap_model_quality(df_fname$fname, n_participants=n_subj, n_boot_samples=n_boot_samples, use_dynamics = T )

          model_quality %<>% dplyr::bind_cols(sim_info, . )
          
          saveRDS(model_quality, file = fname_out_dynamics)
        }
        
        if (!file.exists(fname_out_vanilla)) 
        {
          model_quality <- bootstrap_model_quality(df_fname$fname, n_participants=n_subj, n_boot_samples=n_boot_samples )
          
          model_quality %<>% dplyr::bind_cols(sim_info, . )
          
          saveRDS(model_quality, file = fname_out_vanilla)
        }
        
        
    }
  
    NULL
}, .parallel = run_parallel)

