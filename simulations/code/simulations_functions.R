library(tidyverse)
library(magrittr)

n_unique <- function(x) length(unique(x))

satf <- function(t, intercept, rate, asymptote) {
  ifelse(t > intercept, asymptote*(1-exp(-rate*(t-intercept))), 0)
}

p_yes <- function(mu, criterion) {
  1 - pnorm(criterion, mean = mu, sd = 1)
}

r_yes <- function(n, p_yes) {
  stopifnot(length(n) == 1 && length(p_yes) == 1)
  runif(n) <= p_yes
}

satf_gen_cond_mr <- function(mu, criterion, time, n, label) {
  stopifnot(length(mu) == length(criterion) && length(criterion) == length(time))
  trial = rep(1:n, each = length(time))
  interval <- 1:length(time)
  
  df <- data.frame(condition = label, interval = interval, time = time, trial = trial, mu = mu, 
                   criterion = criterion, delta_evidence = c(0, diff(mu)) )
  
  # generate start position (for index=1), and delta due to noise (for index > 1) 
  df$evidence_rel_position <- rnorm(nrow(df), mean = 0, sd = 1)
  # determine the relative evidence positions (modulo the provided mus)
  df %<>% group_by(trial) %>% dplyr::mutate(evidence_rel_position = cumsum(evidence_rel_position)/sqrt(interval),
                                            evidence_position = evidence_rel_position + mu#,
                                            #cor = cor(evidence_rel_position, lag(evidence_rel_position), use = "complete.obs")
  )
  # NOTE: The above assumptions about the generative process for MR SAT data results in the serial correlation between
  #       evidence_rel_position values being abound ~0.6. Other assumptions (lower or higher SD) would result in a different
  #       correlation coefficient.
  
  df$response <- df$evidence_position > criterion
  df %>% ungroup() %>% dplyr::select(condition, interval, time, trial, response)
}

satf_gen_cond_sr <- function(mu, criterion, time, n, label) {
  stopifnot(length(mu) == length(criterion) && length(criterion) == length(time))
  trial = rep(1:n, each = length(time))
  interval <- 1:length(time)
  df <- data.frame(condition = label, interval = interval, time = time, trial = trial, mu = mu, criterion = criterion)
  df$evidence_position <- rnorm(nrow(df), mean = mu, sd = 1)
  df$response <- df$evidence_position > criterion
  df %>% dplyr::select(condition, interval, time, trial, response)
}

satf_gen <- function(time, n, intercept, rate, asymptote, label = "condition1", fn_satf_gen_cond = NULL) {
  dprime <- satf(time, intercept, rate, asymptote)
  criterion <- 0.5*dprime
  df0 <- fn_satf_gen_cond(mu = dprime*0, criterion = criterion, time = time, n = n, label = label)
  df1 <- fn_satf_gen_cond(mu = dprime, criterion = criterion, time = time, n = n, label = label)
  df0$is_signal <- 0
  df1$is_signal <- 1
  df1$trial = df1$trial + max(df0$trial)
  rbind(df0, df1)
}

satf_gen_mr <- function(time, n, intercept, rate, asymptote, label) {
  satf_gen(time = time, n = n, intercept = intercept, rate = rate, asymptote = asymptote, label = label, fn_satf_gen_cond = satf_gen_cond_mr)
}
satf_gen_sr <- function(time, n, intercept, rate, asymptote, label = "condition1", fn_satf_gen_cond = satf_gen_cond) {
  satf_gen(time = time, n = n, intercept = intercept, rate = rate, asymptote = asymptote, label = label, fn_satf_gen_cond = satf_gen_cond_sr)
}

generate_responses <- function(n, time, intercepts, rates, asymptotes, conditions = 2, fn_responses)
{
  stopifnot( is.function(fn_responses) )
  
  responses1 <- fn_responses(time = time, n = n, intercept = intercepts[1], rate = rates[1], 
                             asymptote = asymptotes[1], label = "condition1")
  if (conditions == 1) {
    responses <- responses1
    
  } else if (conditions == 2) {
    responses2 <- fn_responses(time = time, n = n, intercept = intercepts[2], rate = rates[2], 
                               asymptote = asymptotes[2], label = "condition2")
    responses <- rbind(responses1, responses2)
  }
  
  responses
}

compute_accuracy_stats <- function(data)
{
  calc <- data %>% group_by(condition, interval, is_signal) %>% 
                   dplyr::summarize(p_yes = mean(response), time = mean(time), n = length(response))
  calc$correction <- 'none'

  # 'correct' extreme values
  is_zero <- calc$p_yes == 0
  is_one <- calc$p_yes == 1
  
  if (any(is_zero)) {
    calc$p_yes[is_zero] <- with(calc[is_zero,], 0.5/n)
    calc$correction[is_zero] <- 'extreme'
  }
  if (any(is_one)) {
    calc$p_yes[is_one] <- with(calc[is_one,], (n-0.5)/n)
    calc$correction[is_one] <- 'extreme'
  }
  
  calc_signal <- calc %>% filter(is_signal == 1) %>% 
                    dplyr::select(condition, interval, time, hit = p_yes, 
                                  n_signal = n, hit.correction = correction) %>% 
                    dplyr::mutate(miss = 1-hit)
  calc_noise <- calc %>% filter(is_signal == 0) %>% 
                    dplyr::select(condition, interval, time, fa = p_yes, 
                                  n_noise = n, fa.correction = correction) %>% 
                    dplyr::mutate(creject = 1-fa)

  responses <-
  left_join(calc_signal, calc_noise, by = c("condition", "interval", "time")) %>% 
            dplyr::select(condition,
                          interval, time, hit, miss, hit.correction, n_signal, 
                          fa, creject, fa.correction, n_noise, 
                          )
  
  n_signal <- with(responses, hit + miss)
  n_noise <- with(responses, fa + creject)
  p_hit <- responses$hit/n_signal
  p_fa <- responses$fa/n_noise
  
  responses$dprime <- qnorm(p_hit) - qnorm(p_fa)
  responses$criterion <- -0.5*(qnorm(p_hit) + qnorm(p_fa))

  responses
}


summarize_fit <- function(fit, pc, n_optim_attempts)
{
  convergence <- fit$fit[c('iterations', 'counts', 'convergence')] %>% unlist()
  
  model_fit <- summary.SATcurve(fit)
  model_fit %<>% cbind( t(convergence) ) %>% cbind(n_optim_attempts = n_optim_attempts) 
  
  model_fit$model <- sapply(pc, n_unique) %>% paste(collapse = "-")
  model_fit %<>% dplyr::rename(est_intercept1 = incp1, est_rate1 = rate1, est_asymptote1 = asym1)
  if ("incp2" %in% colnames(model_fit)) {
      model_fit %<>% dplyr::rename(est_intercept2 = incp2, est_rate2 = rate2, est_asymptote2 = asym2)
  }
  model_fit$method %<>% as.character()
  
  model_fit
}

create_par_constraints <- function(pc, start_values = NULL)
{
  epsilon <- 0.1
  bounds <- list( asym = c(0.1, 5), rate = c(0.1, 10), incp = c(0.1, 2))
  
  if (is.null(start_values)) {
      start_values <- c(asym = runif(1, min = bounds$asym[1]+epsilon, max = bounds$asym[2]-epsilon), 
                        rate = runif(1, min = bounds$rate[1]+epsilon, max = bounds$rate[2]-epsilon),
                        incp = runif(1, min = bounds$incp[1]+epsilon, max = bounds$incp[2]-epsilon))
  } else {
      constrain <- function(val, lower, upper) ifelse(val <= lower, lower+epsilon,
                                                      ifelse(val >= upper, upper-epsilon, val)
                                                      )
      start_values[['asym']] %<>% constrain(bounds$asym[1], bounds$asym[2])
      start_values[['rate']] %<>% constrain(bounds$rate[1], bounds$rate[2])
      start_values[['incp']] %<>% constrain(bounds$incp[1], bounds$incp[2])
  }
  
  par_constraints <- get.param(pc, auto.asym = FALSE,
                                asym = c(start_values[["asym"]], bounds$asym[1], bounds$asym[2]),
                                rate = c(start_values[["rate"]], bounds$rate[1], bounds$rate[2]),
                                incp = c(start_values[["incp"]], bounds$incp[1], bounds$incp[2]) )
  par_constraints
}

mrsat_fitcurve <- function(data, pc = list(asym = c(1, 2), rate = c(1, 2), incp = c(1, 2)), 
                           n_optim_attempts = 12, start_values = NULL)
{
    data <- data %>%  dplyr::select(bin = interval, hit, hit.correction,
                                    hit.denom = n_signal, fa, fa.correction, fa.denom = n_noise, 
                                    lags = time, dprimes = dprime, condition)
  
    # repeat optimization with random start values until it converges
    optim_algos = c("optim") # "acp", "nlminb" (optim to perform better than these two, and is faster than 'acp') / "nlm" (no bounding constraints)
    optim_algos <- rep(optim_algos, length.out = n_optim_attempts)
  
    best_fit <- NULL
    for (i in 1:length(optim_algos))
    {
        if (i == 1) {
            par_constraints <- create_par_constraints(pc, start_values)
        } else {
            par_constraints <- create_par_constraints(pc)
        }
      
        cur_fit <-
          tryCatch({ 
            fit.SATcurve(data, par.cond = pc, params = par_constraints, 
                         maxit = 10^6, opt = optim_algos[i], rep = 1) 
          }, error = function(e) NULL)
        
        if ( is.null(cur_fit) ) {
            n_optim_attempts <- n_optim_attempts - 1
            
        } else if ( is.null(best_fit) || cur_fit$R2 > best_fit$R2 ) {
            best_fit <- cur_fit
        }
    }

    summarize_fit(best_fit, pc, n_optim_attempts)
}

fit_model_one_condition <- function(responses) 
{
  acc_stat <- compute_accuracy_stats(data = responses)

  res <- mrsat_fitcurve(acc_stat, pc = list(asym = c(1), rate = c(1), incp = c(1)) )
  res %<>% cbind(data.frame(intercept1 = intercept[1], rate1 = rate[1], asymptote1 = asymptote[1]))
  res %<>% dplyr::select(model, intercept1:asymptote1, est_asymptote1:est_intercept1, R2:n_optim_attempts)
  stopifnot(all(res$model == "1-1-1"))
  
  res
}


fit_models_two_conditions <- function(responses, check_R2_ordering)
{
  `%<=%` <- function(lhs, rhs) {
      sapply(1:nrow(lhs), function(i) {
        all( round(lhs$R2[i], digits = 2) <= round(rhs$R2, digits = 2) )
      }) %>% all()
  }

  acc_stat <- compute_accuracy_stats(data = responses)

  n_max_iter = 4
  cur_iter = 1
  while (TRUE)
  {
    r_start <- mrsat_fitcurve(acc_stat, pc = list(asym = c(1, 1), rate = c(1, 1), incp = c(1, 1)) )
    start_values <- with(r_start, c(asym = est_asymptote1, rate = est_rate1, incp = est_intercept1))
    
    cur_mrsat_fitcurve <- function(pc) mrsat_fitcurve(acc_stat, pc = pc, n_optim_attempts = cur_iter, start_values = start_values)
    r111 <- cur_mrsat_fitcurve( pc = list(asym = c(1, 1), rate = c(1, 1), incp = c(1, 1)) )
    r112 <- cur_mrsat_fitcurve( pc = list(asym = c(1, 1), rate = c(1, 1), incp = c(1, 2)) )
    r121 <- cur_mrsat_fitcurve( pc = list(asym = c(1, 1), rate = c(1, 2), incp = c(1, 1)) )
    r122 <- cur_mrsat_fitcurve( pc = list(asym = c(1, 1), rate = c(1, 2), incp = c(1, 2)) )
    r211 <- cur_mrsat_fitcurve( pc = list(asym = c(1, 2), rate = c(1, 1), incp = c(1, 1)) )
    r221 <- cur_mrsat_fitcurve( pc = list(asym = c(1, 2), rate = c(1, 2), incp = c(1, 1)) )
    r212 <- cur_mrsat_fitcurve( pc = list(asym = c(1, 2), rate = c(1, 1), incp = c(1, 2)) )
    r222 <- cur_mrsat_fitcurve( pc = list(asym = c(1, 2), rate = c(1, 2), incp = c(1, 2)) )
    
    # make sure that the estimates are sensible, or re-estimate all models again
    # -> the purpose of re-estimation is to exploit all models' R2s to minimize the number of re-runs inside mrsat_fitcurve(),
    #    since 4 should be enough most of the time, at least if we double-check that the results are all sensible
    estimates_consistent <- r111 %<=% rbind(r112, r121, r211) &&
      r112 %<=% rbind(r122, r212) && r122 %<=% r222 && r212 %<=% r222 &&
      r121 %<=% rbind(r221, r122) && r221 %<=% r222 && r122 %<=% r222 &&
      r211 %<=% rbind(r221, r212) && r221 %<=% r222 && r212 %<=% r222
    if (estimates_consistent || cur_iter >= n_max_iter ) {
      break;
    }
    cur_iter <- cur_iter + 1
  }
  res <- dplyr::bind_rows(r111, r112, r121, r122, r211, r221, r212, r222)

  # res %<>% cbind(data.frame(intercept1 = intercept[1], rate1 = rate[1], asymptote1 = asymptote[1], 
  #                           intercept2 = intercept[2], rate2 = rate[2], asymptote2 = asymptote[2]))
  res %<>% dplyr::select(model, est_asymptote1:est_intercept2, R2:n_optim_attempts)

  res
}



create_lockfile <- function(fname) {
  fname <- paste0(fname, "_lock")
  if( !file.exists(fname) || difftime( Sys.time(), file.info(fname)$mtime, units = "days") > 1 ) {
      write("x", file = fname)
      return (TRUE)
  }
  return (FALSE)
}

remove_lockfile <- function(fname) {
  fname <- paste0(fname, "_lock")
  if( file.exists(fname) ) {
    unlink(fname)
  }
}
