
source("./simulations_functions.R")
library(mrsat)
library(doMC)
doMC::registerDoMC(parallel::detectCores(logical=FALSE)-1)

# run_estimation =  T
# estimate_single_subj_models = F
# run_parallel = run_estimation

base_path <- "../simulation_output"
fname_template_data_mr <- file.path(base_path, "1_responses", "mr", "sim_%d-%d_%s.rds")
fname_template_data_sr <- file.path(base_path, "1_responses", "sr", "sim_%d-%d_%s.rds")


#options(warn=2)

simulate_experiments <- function(n_simulations, n_participants, n_trials_per_interval,
                                  time, mrsat, sim_id, fn_get_subject_parameters) 
{
  if (mrsat) {
      fname_template <- fname_template_data_mr
      fn_gen_responses <- satf_gen_mr
  } else {
      fname_template <- fname_template_data_sr
      fn_gen_responses <- satf_gen_sr
  }
  
  fname <- sprintf(fname_template, n_participants, n_trials_per_interval, sim_id)
  if (file.exists(fname) || !create_lockfile(fname)) {
      return (NULL)
  }

  simulate_experiment <- function(n_trials_per_interval, time, fn_get_subject_parameters, simulation_id)
  {
      responses_lst <- plyr::llply(1:n_participants, function(i) {
          par <- fn_get_subject_parameters(sample = TRUE)
          responses <- generate_responses(n = n_trials_per_interval, time = time, 
                                          intercepts = par$intercepts, rates = par$rates, 
                                          asymptotes = par$asymptotes, conditions = 2, 
                                          fn_responses = fn_gen_responses)
          responses_with_id <- cbind(participant_id=i, responses)
          attr(responses_with_id, "pars") <- c(list(n = n_trials_per_interval, time = time), par)
          responses_with_id
      })

      responses <- plyr::ldply(responses_lst, function(elem) elem )
      attr(responses, "pars") <- plyr::llply(responses_lst, function(elem) attr(elem, "par") )

      responses
  }

  responses_lst <- plyr::llply(1:n_simulations, function(j) {
      responses_j <- simulate_experiment(n_trials_per_interval = n_trials_per_interval, time = time, fn_get_subject_parameters, j)
      responses_with_id <- cbind(simulation_id=j, responses_j)
      attr(responses_with_id, "pars") <- attr(responses_j, "pars")
      responses_with_id
  }, .progress = "text")
  
  print(c("Finished", date()))
  saveRDS(responses_lst, file = fname)
  remove_lockfile(fname)
}


###################################################################################################################################

# log-normal reparameterization:
#   mu = log(mean / sqrt(1 + sd^2/mean^2))
#   sigma = log(1 + sd^2/mean^2)
dlnorm2 <- function(x, mean, sd) dlnorm(x=x, meanlog = log(mean / sqrt(1 + sd^2/mean^2)), sdlog = log(1 + sd^2/mean^2))
qlnorm2 <- function(p, mean, sd) qlnorm(p=p, meanlog = log(mean / sqrt(1 + sd^2/mean^2)), sdlog = log(1 + sd^2/mean^2))
rlnorm2 <- function(n, mean, sd) rlnorm(n=n, meanlog = log(mean / sqrt(1 + sd^2/mean^2)), sdlog = log(1 + sd^2/mean^2))

# plot(function(x) dlnorm2(x, mean=.8, sd=.35), xlim = c(0, 2) )
# qlnorm2(c(.01, .99), mean=.6, sd=.35)


get_subject_parameters <- function(mean_delta_asymptote, mean_delta_invrate, mean_delta_intercept,
                                   sample = TRUE, sample_fixdelta = FALSE)
{
  mean_asymptote =  2.25; sd_asymptote = 1.25; # lognormal: 98% in approx. [1.05; 3.68]
  mean_invrate   = .8;    sd_invrate = .35;    # lognormal: 98% in approx. [0.49; 1.02]
  mean_intercept = .8;    sd_intercept = .35;  # lognormal: 98% in approx. [0.49; 1.02]
  
  parfn_asymptote = function(n) rlnorm2(n = 1, mean = mean_asymptote, sd = sd_asymptote)
  parfn_invrate   = function(n) rlnorm2(n = 1, mean = mean_invrate, sd = sd_invrate)
  parfn_intercept = function(n) rlnorm2(n = 1, mean = mean_intercept, sd = sd_intercept)
  
  if (sample) {
      avg_asymptote = parfn_asymptote()
      avg_invrate = parfn_invrate()
      avg_intercept = parfn_intercept()
  }
  
  # Always returns 0 when deltas are 0, because the kind of lack of effect of experimental condition we are typically
  # interested translate to the lack of an effect in all (or most?) participants, and not to a distribution of effects
  # over participants which sums to 0.
  sd_delta_asymptote = mean_delta_asymptote/2
  sd_delta_invrate = mean_delta_invrate/2
  sd_delta_intercept = mean_delta_intercept/2
  parfn_delta_asymptote = function() { ifelse(mean_delta_asymptote == 0, 0, rlnorm2(n = 1, mean = mean_delta_asymptote, sd = sd_delta_asymptote)) }
  parfn_delta_invrate   = function() { ifelse(mean_delta_invrate == 0, 0, rlnorm2(n = 1, mean = mean_delta_invrate, sd = sd_delta_invrate)) }
  parfn_delta_intercept = function() { ifelse(mean_delta_intercept == 0, 0, rlnorm2(n = 1, mean = mean_delta_intercept, sd = sd_delta_intercept)) }
  
  if (sample) {
      if (sample_fixdelta) {
        delta_asymptote = mean_delta_asymptote
        delta_invrate = mean_delta_invrate
        delta_intercept = mean_delta_intercept
        
      } else {
        delta_asymptote = parfn_delta_asymptote()
        delta_invrate = parfn_delta_invrate()
        delta_intercept = parfn_delta_intercept()
        
      }
    
      intercepts = avg_intercept + c(-.5, .5) * delta_intercept
      invrates = avg_invrate + c(-.5, .5) * delta_invrate
      asymptotes = avg_asymptote + c(-.5, .5) * delta_asymptote
      stopifnot(all( c(intercepts, invrates, asymptotes) > 0 ))
  }
  
  mean_intercepts = mean_intercept + c(-.5, .5) * mean_delta_intercept
  mean_invrates = mean_invrate + c(-.5, .5) * mean_delta_invrate
  mean_asymptotes = mean_asymptote + c(-.5, .5) * mean_delta_asymptote
  stopifnot(all( c(mean_intercepts, mean_invrates, mean_asymptotes) > 0 ))
  
  ret <- list( par_asymptote = c(mean=mean_asymptote, sd=sd_asymptote), 
               par_invrate = c(mean=mean_invrate, sd=sd_invrate), 
               par_intercepts = c(mean=mean_intercept, sd=sd_intercept),
               par_delta_asymptote = c(mean=mean_delta_asymptote, sd=sd_delta_asymptote), 
               par_delta_invrate = c(mean=mean_delta_invrate, sd=sd_delta_invrate), 
               par_delta_intercepts = c(mean=mean_delta_intercept, sd=sd_delta_intercept),
               mean_asymptotes = mean_asymptotes,
               mean_invrates = mean_invrates,
               mean_intercepts = mean_intercepts
              )
  
  if (sample) {
    ret$intercepts <- intercepts
    ret$rates <- 1/invrates
    ret$asymptotes <- asymptotes
  }
  ret
}



time = seq(0, 5.6, 0.35)
n_simulations = 100
n_participants = 20

for (n_obs_per_interval_per_cond in c( 20, 80, 50) %>% sample()) {
  for (simulate_mrsat in (c(T, F) %>% sample()) ) {
    for (delta_asymptote in (c(0, .75, .5, .25, .125) %>% sample()) ) {
      for (delta_invrate in (c(.0, .15, .1, .05, .025) %>% sample()) ) {
        for(delta_intercept in (c(0, .15, .1, .05, .025) %>% sample()) ) {

            if ( ((delta_asymptote>0) + (delta_invrate>0) + (delta_intercept>0)) > 1 ) {
                next;
            }
          
            fn_get_subject_parameters <- function(...)
                  get_subject_parameters(mean_delta_asymptote = delta_asymptote,
                                         mean_delta_invrate = delta_invrate,
                                         mean_delta_intercept = delta_intercept,
                                         ...)

            simulation_id = sprintf("vardelta_%0.3f-%0.3f-%0.3f", delta_asymptote, delta_invrate, delta_intercept)

            print(c(n_obs_per_interval_per_cond, simulation_id))

            simulate_experiments(n_simulations = n_simulations,
                                 n_participants = n_participants,
                                 n_trials_per_interval = n_obs_per_interval_per_cond,
                                 time = time, mrsat = simulate_mrsat, sim_id = simulation_id,
                                 fn_get_subject_parameters = function(...)
                                          fn_get_subject_parameters(sample_fixdelta = FALSE, ...)
                                  )
          }
        }
      }
    }
  }
