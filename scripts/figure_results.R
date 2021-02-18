
library(tidyverse)
library(magrittr)
library(ggplot2)

options(dplyr.summarise.inform=FALSE)


summarize_percentages <- function(df)
{
  # df$one_param_nonzero <- with(df, delta_asym | delta_invrate | delta_incp)

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



fnames <- dir("/media/pavel3/STORAGE/satf_project_simulations/simulation_output/model_comparison/results_model_comparison/", pattern = "*.rds", full.names = T)
df <- plyr::ldply(fnames, function(fname) readRDS(fname))

head(df)
tail(df)

best_t_val222 <- df %>% subset(t_val222) %>% summarize_percentages() %>% bind_cols(method = "t_val222")
best_forward_selection_t <- df %>% subset(best_forward_selection_t) %>% summarize_percentages() %>% bind_cols(method = "forward_selection_t")
best_backward_selection_t <- df %>% subset(best_backward_selection_t) %>% summarize_percentages() %>% bind_cols(method = "backward_selection_t")
best_forward_selection_adjR2_t <- df %>% subset(best_forward_selection_adjR2_t) %>% summarize_percentages() %>% bind_cols(method = "forward_selection_adjR2_t")
best_backward_selection_adjR2_t <- df %>% subset(best_backward_selection_adjR2_t) %>% summarize_percentages() %>% bind_cols(method = "backward_selection_adjR2_t")

summary_by_method <-
      best_t_val222 %>% 
      bind_rows(best_forward_selection_t) %>% 
      bind_rows(best_backward_selection_t) %>% 
      bind_rows(best_forward_selection_adjR2_t) %>% 
      bind_rows(best_backward_selection_adjR2_t) %>%
      ungroup()

summary_by_method$one_param_nonzero <- with(summary_by_method, delta_asym!=0 | delta_invrate!=0 | delta_incp!=0)

# summary_by_method$delta_incp %<>% multiply_by(1000)
# summary_by_method$delta_invrate %<>% multiply_by(1000)


for (opt_response_type in c("mr", "sr"))
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


method_labels <- c("t-tests on\n2-2-2 model\nestimates", 
                    "forward selection\n with t-tests", "forward selection\nwith adjusted R^2\nand t-tests",
                    "backward selection\n with t-tests", "backward selection\nwith adjusted R^2\nand t-tests")

summary_by_method_long$param_name <- summary_by_method_long$positive_rate_name
# %>% 
#         dplyr::recode("asymptote"="asymptote (d' units)",
#                       "invrate"="1/rate (ms)",
#                       "intercept"="intercept (ms)") %>%
#         factor(levels = c("asymptote (d' units)", "1/rate (ms)", "intercept (ms)") )

summary_by_method_long$n_responses_bysubj %<>% sprintf("N = %s", .)


# p <-
#   summary_by_method_long %>% subset(response_type == opt_response_type) %>% 
#   ggplot(aes(delta_val, positive_rate, color = method, linetype = method, shape = method)) + geom_point() + geom_line() + 
#   facet_wrap(param_name~n_responses_bysubj, scales = "free")
# p

p <-
summary_by_method_long %>% subset(response_type == opt_response_type & param_name == "asymptote") %>%
  ggplot(aes(delta_val, positive_rate, color = method, linetype = method, shape = method)) + geom_point() + geom_line() +
          facet_grid(~n_responses_bysubj, scales = "free")

p

p <-
p + scale_color_manual(name = "",
                        breaks = c("t_val222", 
                                    "forward_selection_t", "forward_selection_adjR2_t",
                                    "backward_selection_t", "backward_selection_adjR2_t"),
                         values = c("#F8766D", "#00BA38", "#00BA38", "#619CFF", "#619CFF"),
                         labels = method_labels)

p <-
p + scale_linetype_manual(name = "",
                          breaks = c("t_val222", 
                                  "forward_selection_t", "forward_selection_adjR2_t",
                                  "backward_selection_t", "backward_selection_adjR2_t"),
                          values = c("solid", "solid", "dashed", "solid", "dashed"),
                          labels = method_labels)

p <-
p + scale_shape_manual(name = "",
                       breaks = c("t_val222", 
                                   "forward_selection_t", "forward_selection_adjR2_t",
                                   "backward_selection_t", "backward_selection_adjR2_t"),
                      values = c("circle", "circle", "triangle", "circle", "triangle"),
                      labels = method_labels
                      )

p <-
p + theme_bw() + theme(legend.position = "top", legend.box="vertical", legend.margin=margin()) + #
    scale_y_continuous(labels = scales::percent, breaks = c(0,.25,.5,.75,1)) + #limits = c(0,1), 
    xlab("Δ parameter") +
    ylab("Est. probability of detection of positive difference")

p

z = 7
ggsave(p + theme(panel.grid.minor = element_blank()), file = sprintf("../figures/results/results_%s.pdf", opt_response_type), width=z, height = 0.6*z, device = cairo_pdf)
}

