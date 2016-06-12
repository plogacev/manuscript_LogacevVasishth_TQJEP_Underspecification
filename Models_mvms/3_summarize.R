
library(plyr)
library(dplyr)
library(magrittr)
library(mvmstan)
library(loo)
library(doMC)
library(ggplot2)
library(tidyr)
library(rstan)
library(TeachingDemos)
library(coda)


doMC::registerDoMC(6)
rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())
source("./aux_functions.R")

(load("./fit_nouspec_mlvl.rda"))
(load("./fit_uspec_amb_mlvl.rda"))
#(load("./fit_uspec_amb2groups_mlvl.rda"))


model_waic <- function(model) extract_log_lik(model) %>% waic
cmp <- list()
cmp$'nouspec_misparse-nouspec' <- compare(model_waic(nouspec_fit_misparse), model_waic(nouspec_fit_full))
cmp$'nouspec_misretrieval-nouspec' <- compare(model_waic(nouspec_fit_misretrieval), model_waic(nouspec_fit_full))
cmp$'nouspec-pspec_amb' <- compare(model_waic(nouspec_fit_full), model_waic(pspec_amb_fit))
cmp$'nouspec-nspec_amb' <- compare(model_waic(nouspec_fit_full), model_waic(nspec_amb_fit))
cmp$'nspec_amb-pspec_amb' <- compare(model_waic(nspec_amb_fit), model_waic(pspec_amb_fit))
save(cmp, file="./results/comparisons.rda")



(load("./SwetsDataForModelling.rda"))
d.qrc.stan <- subset(d.qrc.stan, response_RT > 700 & reading_time/crit_region_cnt < 2000)
d <- d.qrc.stan

extract_fit_params <- function(est, parameters)
{
  selected_cols <- c(mean='mean', lower='2.5%', upper='97.5%')
  selected_rows <- rownames(est) %in% parameters
  est %<>% .[selected_rows, selected_cols] %>% as.data.frame()
  colnames(est) <- names(selected_cols)
  fmt <- with(est, sprintf("%0.2f [%0.2f-%0.2f]", mean, lower, upper)) %>% gsub("0\\.", ".", .)
  names(fmt) <- rownames(est)
  fmt <- fmt[parameters]
  names(fmt) <- parameters
  fmt
}


extract_fit_info <- function(mean_region_cnt, fit)
{
  extract_prob <- function(name) {
    if(name %in% fit@model_pars) {
      samples <- (extract(fit, name)[[1]])[,1]
      avg <- mean(samples)
      #hpd <- emp.hpd(samples)
      hpd <- HPDinterval(as.mcmc(samples))
      fmt <- sprintf("%0.2f [%0.2f-%0.2f]", avg, hpd[1], hpd[2]) %>% gsub("0\\.", ".", .)
    } else {
      fmt <- NA
    }
    fmt
  }
  
  # extract probability estimates for the group
  parameters <- c("p_uspec_hyper","p_arbitrary_attachment_hyper",
                  "p_retrieval_fail_hyper", "p_att_n1_hyper", 
                  "p_guess_yes_hyper")
  est_fmt <- sapply(parameters, extract_prob)
  
  # extract time estimates *by participant* (this is because the mean RT for the group can't simply be computed from mean scale and mean shape for the group)
  scale_subj <- extract(fit, "scale_subj")[[1]]
  shape_reading_base_subj <- extract(fit, "shape_reading_base_subj")[[1]]
  shape_delta_attach_N2_subj <- extract(fit, "shape_delta_attach_N2_subj")[[1]]
  shape_delta_attach_N1_subj <- extract(fit, "shape_delta_attach_N1_subj")[[1]]
  shape_response_base_subj <- extract(fit, "shape_response_base_subj")[[1]]
  shape_response_delta_guess_subj <- extract(fit, "shape_response_delta_guess_subj")[[1]]
  
  # extract time estimates
  fmt_time <- function(time) sprintf("%0.0f [%0.0f-%0.0f]", mean(time), quantile(time, .025), quantile(time, .975))
  
  time_uspec <- (apply(scale_subj*shape_reading_base_subj, 1, mean) * mean_region_cnt)
  fmt_time_uspec <- time_uspec %>% fmt_time
  fmt_time_n2 <- (time_uspec + apply(scale_subj*(shape_delta_attach_N2_subj), 1, mean)) %>% fmt_time
  fmt_time_n1 <- (time_uspec + apply(scale_subj*(shape_delta_attach_N2_subj+shape_delta_attach_N1_subj), 1, mean)) %>% fmt_time
  fmt_time_inf_resp <- apply(scale_subj*shape_response_base_subj, 1, mean) %>% fmt_time
  fmt_time_guess <- apply(scale_subj*(shape_response_base_subj + shape_response_delta_guess_subj), 1, mean) %>% fmt_time
  times <- c(time_uspec=fmt_time_uspec, time_n2=fmt_time_n2, time_n1=fmt_time_n1, time_inf_resp=fmt_time_inf_resp, time_guess=fmt_time_guess)
  est_fmt %<>% c(., times)
  
  est_waic <- extract_log_lik(fit) %>% waic
  est_fmt %<>% c(., with(est_waic, c(#elpd_waic = sprintf("%0.0f (SE=%0.0f)", elpd_waic, se_elpd_waic),
    #p_waic = sprintf("%0.0f (SE=%0.0f)", p_waic, se_p_waic),
    waic = sprintf("%0.1f (SE=%0.1f)", waic, se_waic))))
  est_fmt
}


mean_region_cnt <- mean(d$crit_region_cnt)
est_nouspec_fit_misparse <- extract_fit_info(mean_region_cnt, nouspec_fit_misparse)
est_nouspec_fit_misretrieval <- extract_fit_info(mean_region_cnt, nouspec_fit_misretrieval)
est_nouspec_fit_full <- extract_fit_info(mean_region_cnt, nouspec_fit_full)
est_pspec_amb_fit <- extract_fit_info(mean_region_cnt, pspec_amb_fit)
est_nspec_amb_fit <- extract_fit_info(mean_region_cnt, nspec_amb_fit)
#est_pspec_amb2groups_fit <- extract_fit_info(d.qrc.stan, pspec_amb2groups_fit)
#est_nspec_amb2groups_fit <- extract_fit_info(d.qrc.stan, nspec_amb2groups_fit)

table_models <- cbind(est_nouspec_fit_misparse, est_nouspec_fit_misretrieval,
                      est_nouspec_fit_full, est_pspec_amb_fit, est_nspec_amb_fit)

table_models['time_uspec', 1:3] <- NA
table_models['time_guess', 'est_nouspec_fit_misparse'] <- NA

save(table_models, file="./results/table_models.rda")




d_nouspec <- plot_predictions(d.qrc.stan, nouspec_fit_full) %>% attr(., "data")
d_pspec <- plot_predictions(d.qrc.stan, pspec_amb_fit) %>% attr(., "data")
d_nspec <- plot_predictions(d.qrc.stan, nspec_amb_fit) %>% attr(., "data")

d_nouspec$model <- "nouspec"
d_pspec$model <- "pspec"
d_nspec$model <- "nspec"
d_predictions <- rbind(d_nouspec, d_pspec, d_nspec)
save(d_predictions, file="./results/predictions.rda")

extract_log_lik(nouspec_fit_misparse) %>% waic
f <- function(x) (waic(extract_log_lik(x)))

compare(f(nouspec_fit_misparse), f(nouspec_fit_misretrieval), f(nouspec_fit_full), f(pspec_amb_fit), f(nspec_amb_fit))
#                          waic     se_waic  elpd_waic se_elpd_waic p_waic   se_p_waic weights 
# f(nouspec_misretrieval)  44532.0    143.2 -22266.0      71.6        225.3     11.3       0.0
# f(nouspec_misparse)      44821.0    150.0 -22410.5      75.0        166.0      8.3       0.0
# f(nouspec_full)          44248.1    131.6 -22124.1      65.8        207.6     12.1       1.0
# f(non-specification)     44259.6    131.7 -22129.8      65.8        231.6     11.5       0.0
# f(partial specification) 44263.1    131.6 -22131.6      65.8        229.8     11.6       0.0
