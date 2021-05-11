
library(tidyverse)
library(magrittr)
library(ggplot2)

options(dplyr.summarise.inform=FALSE)


summarize_percentages <- function(df)
{
  df %>% 
    dplyr::mutate(model_delta_asym_positive = model_delta_asym == "+1", 
                  model_delta_invrate_positive = model_delta_invrate == "+1", 
                  model_delta_icpt_positive = model_delta_icpt == "+1") %>%
    group_by(simulation_name, response_type, n_subjects, n_responses_bysubj, 
             delta_asym, delta_invrate, delta_incp) %>%
    dplyr::summarise(perc_model_delta_asym_positive = mean(model_delta_asym_positive),
                     perc_model_delta_invrate_positive = mean(model_delta_invrate_positive),
                     perc_model_delta_icpt_positive = mean(model_delta_icpt_positive),
                     N = n(),
                     .groups = "keep")
}

summarize_percentages_dynamics <- function(df)
{
  df %>% 
    dplyr::mutate(delta_dynamics = delta_invrate + delta_incp,
                  model_delta_asym_positive = model_delta_asym == "+1", 
                  model_delta_dynamics_positive = model_delta_dynamics == "+1") %>%
    group_by(simulation_name, response_type, n_subjects, n_responses_bysubj, 
             delta_asym, delta_dynamics) %>%
    dplyr::summarise(perc_model_delta_asym_positive = mean(model_delta_asym_positive),
                     perc_model_delta_dynamics_positive = mean(model_delta_dynamics_positive),
                     N = n(),
                     .groups = "keep")
}



# fnames <- dir("/media/pavel3/STORAGE/satf_project_simulations/simulation_output/model_comparison/results_model_comparison/", pattern = "*.rds", full.names = T)
# fnames <- dir("../simulations/simulation_output/3_model_comparison_vanilla", pattern = "*.rds", full.names = T, recursive = T)
# df <- plyr::ldply(fnames, function(fname) readRDS(fname))

df_dirs_vanilla <- tidyr::expand_grid(response_type = c("mr", "sr"), n_subjects = c(10,20,30,40,50) )
df_dirs_vanilla %<>% mutate(dir = sprintf("../simulations/simulation_output/3_model_comparison_vanilla/n%d/%s", n_subjects, response_type) )
df_vanilla <- plyr::ddply(df_dirs_vanilla, "dir", function(df_dir) {
  #print(df_dir)
  fnames <- dir(df_dir$dir, pattern = "*.rds", full.names = T)
  if (length(fnames) == 0)
    return (NULL)
  df <- plyr::ldply(fnames, function(fname) readRDS(fname))
  df$response_type <- df_dir$response_type
  df$n_subjects <- df_dir$n_subjects
  df
})
df_vanilla %<>% dplyr::select(-dir)

df_dirs_dynamics <- tidyr::expand_grid(response_type = c("mr", "sr"), n_subjects = c(10,20,30,40,50) )
df_dirs_dynamics %<>% mutate(dir = sprintf("../simulations/simulation_output/3_model_comparison_dynamics/n%d/%s", n_subjects, response_type) )
df_dynamics <- plyr::ddply(df_dirs_dynamics, "dir", function(df_dir) {
  fnames <- dir(df_dir$dir, pattern = "*.rds", full.names = T)
  if (length(fnames) == 0)
    return (NULL)
  df <- plyr::ldply(fnames, function(fname) readRDS(fname))
  df$response_type <- df_dir$response_type
  df$n_subjects <- df_dir$n_subjects
  df
})
df_dynamics %<>% dplyr::select(-dir)

#head(df)
#tail(df)

best_t_val222 <- df_vanilla %>% subset(t_val222) %>% summarize_percentages() %>% bind_cols(method = "t_val222")
best_forward_selection_t <- df_vanilla %>% subset(best_forward_selection_t) %>% summarize_percentages() %>% bind_cols(method = "forward_selection_t")
best_backward_selection_t <- df_vanilla %>% subset(best_backward_selection_t) %>% summarize_percentages() %>% bind_cols(method = "backward_selection_t")
best_forward_selection_adjR2_t <- df_vanilla %>% subset(best_forward_selection_adjR2_t) %>% summarize_percentages() %>% bind_cols(method = "forward_selection_adjR2_t")
best_backward_selection_adjR2_t <- df_vanilla %>% subset(best_backward_selection_adjR2_t) %>% summarize_percentages() %>% bind_cols(method = "backward_selection_adjR2_t")

summary_by_method <-
      best_t_val222 %>% 
      bind_rows(best_forward_selection_t) %>% 
      bind_rows(best_backward_selection_t) %>% 
      bind_rows(best_forward_selection_adjR2_t) %>% 
      bind_rows(best_backward_selection_adjR2_t) %>%
      ungroup()
summary_by_method$one_param_nonzero <- with(summary_by_method, delta_asym!=0 | delta_invrate!=0 | delta_incp!=0)


best_t_val222_dynamics <- df_dynamics %>% subset(t_val222) %>% summarize_percentages_dynamics() %>% bind_cols(method = "t_val222")
best_forward_selection_t_dynamics <- df_dynamics %>% subset(best_forward_selection_t) %>% summarize_percentages_dynamics() %>% bind_cols(method = "forward_selection_t")
best_backward_selection_t_dynamics <- df_dynamics %>% subset(best_backward_selection_t) %>% summarize_percentages_dynamics() %>% bind_cols(method = "backward_selection_t")

summary_by_method_dynamics <-
  best_t_val222_dynamics %>% 
  bind_rows(best_forward_selection_t_dynamics) %>% 
  bind_rows(best_backward_selection_t_dynamics) %>% 
  ungroup()
summary_by_method_dynamics$one_param_nonzero <- with(summary_by_method_dynamics, delta_asym!=0 | delta_dynamics !=0 )


method_internal_labels <- c("t_val222", 
                            "forward_selection_t", "forward_selection_adjR2_t",
                            "backward_selection_t", "backward_selection_adjR2_t")
method_labels <- c("t-tests on 2-2-2 model estimates", 
                   "forward selection with t-tests", "forward selection with adjusted R^2 and t-tests",
                   "backward selection with t-tests", "backward selection with adjusted R^2 and t-tests")

format_plot <- function(p) {
  p <-
    p + scale_color_manual(name = "",
                           breaks = method_internal_labels,
                           values = c("#F8766D", "#00BA38", "#00BA38", "#619CFF", "#619CFF"),
                           labels = method_labels)
  
  p <-
    p + scale_linetype_manual(name = "",
                              breaks = method_internal_labels,
                              values = c("solid", "solid", "dashed", "solid", "dashed"),
                              labels = method_labels)
  
  p <-
    p + scale_shape_manual(name = "",
                           breaks = method_internal_labels,
                           values = c("circle", "circle", "triangle", "circle", "triangle"),
                           labels = method_labels
    )
  
  p <-
    p + theme_bw() + theme(legend.position = "top", legend.box="vertical", legend.margin=margin()) + #
    scale_y_continuous(labels = scales::percent, limits = c(0,1)) + #, breaks = c(0,.25,.5,.75,1)) + #, 
    #xlab("Δ parameter") +
    ylab("Positive difference proportion")
}



for (opt_response_type in c("mr", "sr"))
{

{
summary_by_method_long <-
        summary_by_method %>% 
            pivot_longer(names_to = "delta_name", values_to = "delta_val", cols = c("delta_asym", "delta_invrate", "delta_incp")) %>%
            pivot_longer(names_to = "positive_rate_name", values_to = "positive_rate", cols = c("perc_model_delta_asym_positive", "perc_model_delta_invrate_positive", "perc_model_delta_icpt_positive"))
summary_by_method_long$delta_name %<>% dplyr::recode("delta_asym"="asymptote", "delta_invrate"="invrate", "delta_incp"="intercept")
summary_by_method_long$positive_rate_name %<>% dplyr::recode("perc_model_delta_asym_positive"="asymptote", "perc_model_delta_invrate_positive"="invrate", "perc_model_delta_icpt_positive"="intercept")
summary_by_method_long %<>% subset(delta_name == positive_rate_name)
summary_by_method_long %<>% subset( !(delta_val == 0 & one_param_nonzero) )
summary_by_method_long$positive_rate_name %<>% factor(levels = c("asymptote","invrate","intercept"))
summary_by_method_long$param_name <- summary_by_method_long$positive_rate_name

summary_by_method_long$n_subjects %<>% sprintf("%s subjects", .)
summary_by_method_long$n_responses_bysubj %<>% sprintf("N = %s", .)
}

{
  summary_by_method_long_dynamics <-
    summary_by_method_dynamics %>% 
    pivot_longer(names_to = "delta_name", values_to = "delta_val", cols = c("delta_asym", "delta_dynamics")) %>%
    pivot_longer(names_to = "positive_rate_name", values_to = "positive_rate", cols = c("perc_model_delta_asym_positive", "perc_model_delta_dynamics_positive"))
  summary_by_method_long_dynamics$delta_name %<>% dplyr::recode("delta_asym"="asymptote", "delta_dynamics"="dynamics")
  summary_by_method_long_dynamics$positive_rate_name %<>% dplyr::recode("perc_model_delta_asym_positive"="asymptote", "perc_model_delta_dynamics_positive"="dynamics")
  summary_by_method_long_dynamics %<>% subset(delta_name == positive_rate_name)
  summary_by_method_long_dynamics %<>% subset( !(delta_val == 0 & one_param_nonzero) )
  summary_by_method_long_dynamics$positive_rate_name %<>% factor(levels = c("asymptote","dynamics"))
  summary_by_method_long_dynamics$param_name <- summary_by_method_long_dynamics$positive_rate_name
  
  summary_by_method_long_dynamics$n_subjects %<>% sprintf("%s subjects", .)
  summary_by_method_long_dynamics$n_responses_bysubj %<>% sprintf("N = %s", .)
}
  
ps <- list()
for (par in unique(summary_by_method_long$param_name)) {
    ps[[par]] <-
        summary_by_method_long %>% 
          subset(response_type == opt_response_type & param_name == par) %>%
          mutate(t_delta_val = ifelse(param_name == "asymptote", delta_val, delta_val*1000)) %>%
          ggplot(aes(t_delta_val, positive_rate, color = method, linetype = method, shape = method)) + geom_point() + geom_line() +
          facet_grid(n_responses_bysubj~n_subjects , scales = "free") + 
          theme(panel.grid.minor = element_blank())
    ps[[par]] %<>% format_plot()
}

ps[["asymptote"]] <- ps[["asymptote"]] + xlab("Δ asymptote (d' units)")
ps[["invrate"]] <- ps[["invrate"]] + xlab("Δ 1/rate (ms)")
ps[["intercept"]] <- ps[["intercept"]] + xlab("Δ intercept (ms)")


ps[["dynamics"]] <-
summary_by_method_long_dynamics %>% 
  subset(response_type == opt_response_type & param_name == "dynamics") %>%
  mutate(t_delta_val = ifelse(param_name == "asymptote", delta_val, delta_val*1000)) %>%
  group_by(t_delta_val, method, n_responses_bysubj, n_subjects) %>%
  dplyr::summarise(max_positive_rate = max(positive_rate)) %>%
  ggplot(aes(t_delta_val, max_positive_rate, color = method, linetype = method, shape = method)) + geom_point() + geom_line(aes(y=max_positive_rate)) +
  facet_grid(n_responses_bysubj~n_subjects , scales = "free") + 
  theme(panel.grid.minor = element_blank())
ps[["dynamics"]] %<>% format_plot()
ps[["dynamics"]] <- ps[["dynamics"]] + xlab("Δ intercept + 1/rate (ms)")



library(cowplot)

px <- 
  plot_grid(ps[["asymptote"]] + theme(legend.position="none"),
            ps[["invrate"]] + theme(legend.position="none"), 
            ps[["intercept"]] + theme(legend.position="none"), 
            ps[["dynamics"]] + theme(legend.position="none"))

legend <- get_legend(ps[["asymptote"]] + theme(legend.position="bottom"))

p <- plot_grid( legend, px, ncol = 1, rel_heights = c(.075, 1))

p

z = 2*7
ggsave(p, file = sprintf("../figures/results/results_%s.pdf", opt_response_type), width=z, height = 0.6*z, device = cairo_pdf)

}

