
library(plyr)
library(dplyr)
library(magrittr)
library(mvmstan)
library(loo)
library(doMC)
library(ggplot2)
library(tidyr)
library(rstan)
doMC::registerDoMC(6)
rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())
source("./aux_functions.R")

# these values are required for switching between different models as values of 0 or 1 for probabilites cause problems with the back-transformation
epsilon <- 1e-06
zero <- epsilon
one <- 1-epsilon

#################
### Load data ###
#################

(load("./SwetsDataForModelling.rda"))

# spot potential outliers:
#d$response_correct <- with(d, ifelse(iv_cond==0, "NA", ifelse(iv_cond==1, ifelse(iv_questN1==response_yes, "corr", "wrong"), ifelse(iv_questN1!=response_yes, "corr", "wrong")  ))) %>% as.factor
#ggplot(d, aes(x=as.factor(subj), y=reading_time/crit_region_cnt, color=response_correct)) + geom_point() #geom_violin() #+ coord_cartesian(ylim=c(10000,20000))

d.qrc.stan <- subset(d.qrc.stan, response_RT > 700 & reading_time/crit_region_cnt < 2000)
d <- d.qrc.stan
as.stan.data.list <- function(d) c(c(n_obs=nrow(d), n_subj=length(unique(d$subj)), n_item=length(unique(d$item))), as.list(d))


#######################################################################################################
### Fit model without any type of underspecification for each participant individually ################
#######################################################################################################

uebermodel <- stan_model(model_name="uebermodel", model_code = read_model("./models/uebermodel_functions.stan","./models/uebermodel_mvm.stan","./models/uebermodel_generated.stan"))
par_fixed_nouspec <- c(p_uspec=zero, p_pspec=zero, p_uspec_unamb=zero, p_arbitrary_attachment=zero, p_arbitrary_attachment_late=zero)

#cur_d <- subset(d, subj == 1)
#(res <- sampling(uebermodel, data=c(as.stan.data.list(cur_d), par_fixed_nouspec), verbose=T, chains=1, cores=1, iter=1000, seed=1234 ))

nouspec_fits <- dlply(d, .(subj), function(cur_d) {
  sampling(uebermodel, data=c(as.stan.data.list(cur_d), par_fixed_nouspec), verbose=T, chains=4, cores=1, iter=1000)
}, .parallel = T)
save(nouspec_fits, file="fits_nouspec_single.rda")


##################################################################################
### Fit hierarchival models without any type of underspecification ################
##################################################################################

# load by-participant fits
(load(file="fits_nouspec_single.rda"))

# define variables with hyperparameters in this model
lognormal_vars <- c("scale", "shape_reading_base", "shape_delta_attach_N2", "shape_delta_attach_N1",
                    "shape_response_base", "shape_response_delta_guess")
prob_vars <- c("p_retrieval_fail", "p_att_n1", "p_guess_yes", "p_arbitrary_attachment")
#"p_uspec", "p_uspec_unamb"

# determine starting values for hyperparameters and random effects from by-participant fits
init_uebermodel_nouspec <- generate_init(nouspec_fits, lognormal_vars, prob_vars)
init_uebermodel_nouspec$shape_delta_attach_N1_hyper %<>% pmax(epsilon)

# disable underspecification and other model features which are not used in this model
n_subj <- length(unique(d$subj))
par_fixed_nouspec <- list(p_pspec=one, 
                          p_uspec_subj=rep(zero,n_subj), p_uspec_unamb_subj=rep(zero,n_subj),
                          p_arbitrary_attachment_subj=rep(zero,n_subj),
                          p_retrieval_fail_subj=rep(zero,n_subj))

# compile multilevel models without underspecification
uebermodel_nouspec_misparse_mlvl <- stan_model(model_code = read_model("./models/uebermodel_functions.stan", "./models/uebermodel_mlvl_nouspec_misparse.stan","./models/uebermodel_generated.stan"))
uebermodel_nouspec_misretrieval_mlvl <- stan_model(model_code = read_model("./models/uebermodel_functions.stan", "./models/uebermodel_mlvl_nouspec_misretrieval.stan","./models/uebermodel_generated.stan"))
uebermodel_nouspec_full_mlvl <- stan_model(model_code = read_model("./models/uebermodel_functions.stan", "./models/uebermodel_mlvl_nouspec_full.stan","./models/uebermodel_generated.stan"))

# prepare data and init values
data_all <- as.stan.data.list(d)
init_list <- list(init_uebermodel_nouspec, init_uebermodel_nouspec, init_uebermodel_nouspec, init_uebermodel_nouspec)

# make sure that model estimation can be run
nouspec_fit_misparse <- sampling(uebermodel_nouspec_misparse_mlvl, data=c(data_all, par_fixed_nouspec), verbose=T, chains=1, iter=10, cores=1, seed=1234)
nouspec_fit_misretrieval <- sampling(uebermodel_nouspec_misretrieval_mlvl, data=c(data_all, par_fixed_nouspec), verbose=T, chains=1, iter=10, cores=1, seed=1234)
nouspec_fit_full <- sampling(uebermodel_nouspec_full_mlvl, data=c(data_all, par_fixed_nouspec), verbose=T, chains=1, iter=10, cores=1, seed=1234)

show_params(nouspec_fit_misparse) %>% round(2)
show_params(nouspec_fit_misretrieval) %>% round(2)
show_params(nouspec_fit_full) %>% round(2)

# fit the model and save the results
n_chains <- 4
n_iter <- 2000
nouspec_fit_misparse <- sampling(uebermodel_nouspec_misparse_mlvl, data=c(data_all, par_fixed_nouspec), verbose=T, 
                                 chains=n_chains, iter=n_iter, cores=n_chains, seed=1234, init=init_list)
nouspec_fit_misretrieval <- sampling(uebermodel_nouspec_misretrieval_mlvl, data=c(data_all, par_fixed_nouspec), verbose=T, 
                                     chains=n_chains, iter=n_iter, cores=n_chains, seed=1234, init=init_list)
nouspec_fit_full <- sampling(uebermodel_nouspec_full_mlvl, data=c(data_all, par_fixed_nouspec), verbose=T, 
                             chains=n_chains, iter=n_iter, cores=n_chains, seed=1234, init=init_list)
  
save(nouspec_fit_misparse, nouspec_fit_misretrieval, nouspec_fit_full, file="fit_nouspec_mlvl.rda")


### Evaluate results ###
(extract_log_lik(nouspec_fit_misparse) %>% waic)
(extract_log_lik(nouspec_fit_misretrieval) %>% waic)
(extract_log_lik(nouspec_fit_full) %>% waic)
show_params(nouspec_fit_full)
# Conclusion: The full model is superior to the models which posit only misparses or only retrieval failures.


###################################################################
### Fit hierarchical models with underspecification ###############
###################################################################

(load(file="fit_nouspec_mlvl.rda"))

# prepare data and init values
n_chains <- 4
n_iterations <- 2000
n_subj <- length(unique(d$subj))
data_all <- as.stan.data.list(d)
init_uebermodel_nouspec <- show_params(nouspec_fit_full)[,'50%'] %>% .[!grepl("lp__", names(.))] %>% par2list()
init_list <- list(init_uebermodel_nouspec, init_uebermodel_nouspec, init_uebermodel_nouspec, init_uebermodel_nouspec)

init_best_fit <- llply(1:4, function(idx_chain) {
  max_idx <- (get_logposterior(x))[[idx_chain]] %>% which.max
  nouspec_fit_full@sim$samples[[idx_chain]] %>% names
  sapply(nouspec_fit_full@sim$samples[[idx_chain]], function(x) x[max_idx]) %>% par2list()
})
init_list <- init_best_fit

# simple underspecification models
uebermodel_uspec_amb_mlvl <- stan_model(model_code = read_model("./models/uebermodel_functions.stan", "./models/uebermodel_mlvl_uspec_amb.stan", "./models/uebermodel_generated.stan"))
par_fixed_uspec <- list(p_uspec_unamb_subj=rep(zero,n_subj))

# sample from the model
nspec_amb_fit <- sampling(uebermodel_uspec_amb_mlvl, data=c(p_pspec=zero, data_all, par_fixed_uspec), verbose=T, iter=n_iterations, 
                      chains=n_chains, cores=n_chains, init=init_list, seed=1234)
pspec_amb_fit <- sampling(uebermodel_uspec_amb_mlvl, data=c(p_pspec=one, data_all, par_fixed_uspec), verbose=T, iter=n_iterations, 
                          chains=n_chains, cores=n_chains, init=init_list, seed=1234)
save(pspec_amb_fit, nspec_amb_fit, file="./fit_uspec_amb_mlvl.rda")

summary(nouspec_fit_full)$summary %>% tail
summary(pspec_amb_fit)$summary %>% tail
summary(nspec_amb_fit)$summary %>% tail

##########################################################################
### Fit hierarchical models with underspecification mixture for groups ###
##########################################################################

(load(file="./fit_uspec_amb_mlvl.rda"))

# prepare data and init values
n_chains <- 4
n_iterations <- 2000
n_subj <- length(unique(d$subj))
data_all <- as.stan.data.list(d)
init_pspec <- show_params(pspec_amb_fit)[,'50%'] %>% .[!grepl("lp__", names(.))] %>% par2list()
init_nspec <- show_params(nspec_amb_fit)[,'50%'] %>% .[!grepl("lp__", names(.))] %>% par2list()

# add default values for the second mixture component of the p_uspec_hyper
init_pspec$p_uspec_hyper %<>% c(.8, .8, .5)
init_nspec$p_uspec_hyper %<>% c(.8, .8, .5)

# set those underspecification probabilites which were set to 0, to a very small number instead
init_nspec$p_uspec_subj %<>% pmax(epsilon)
init_pspec$p_uspec_subj %<>% pmax(epsilon)

# populate initialization lists
init_pspec_list <- list(init_pspec, init_pspec, init_pspec, init_pspec)
init_nspec_list <- list(init_nspec, init_nspec, init_nspec, init_nspec)

# simple underspecification models
uebermodel_uspec_amb2groups_mlvl <- stan_model(model_code = read_model("./models/uebermodel_functions.stan", "./models/uebermodel_mlvl_uspec_amb_2groups.stan", "./models/uebermodel_generated.stan"))
par_fixed_uspec <- list(p_uspec_unamb_subj=rep(zero,n_subj))

# sample from the model
nspec_amb2groups_fit <- sampling(uebermodel_uspec_amb2groups_mlvl, data=c(p_pspec=zero, data_all, par_fixed_uspec), verbose=T, iter=n_iterations, 
                          chains=n_chains, cores=n_chains, init=init_nspec_list, seed=1234)
pspec_amb2groups_fit <- sampling(uebermodel_uspec_amb2groups_mlvl, data=c(p_pspec=one, data_all, par_fixed_uspec), verbose=T, iter=n_iterations, 
                          chains=n_chains, cores=n_chains, init=init_pspec_list, seed=1234)
save(pspec_amb2groups_fit, nspec_amb2groups_fit, file="./fit_uspec_amb2groups_mlvl.rda")


###################################################
### Fit hierarchical models with all parameters ###
###################################################


(load(file="./fit_uspec_amb_mlvl.rda"))

# prepare data and init values
n_chains <- 4
n_iterations <- 2000
n_subj <- length(unique(d$subj))
data_all <- as.stan.data.list(d)
init_pspec <- show_params(pspec_amb_fit)[,'50%'] %>% .[!grepl("lp__", names(.))] %>% par2list()
init_nspec <- show_params(nspec_amb_fit)[,'50%'] %>% .[!grepl("lp__", names(.))] %>% par2list()

# set those underspecification probabilites which were set to 0, to a very small number instead
init_nspec$p_uspec_subj %<>% pmax(epsilon)
init_pspec$p_uspec_subj %<>% pmax(epsilon)

# populate initialization lists
init_pspec_list <- list(init_pspec, init_pspec, init_pspec, init_pspec)
init_nspec_list <- list(init_nspec, init_nspec, init_nspec, init_nspec)

# simple underspecification models
uebermodel_uspec_full_mlvl <- stan_model(model_code = read_model("./models/uebermodel_functions.stan", "./models/uebermodel_mlvl_uspec_full.stan", "./models/uebermodel_generated.stan"))

# sample from the model
nspec_full_fit <- sampling(uebermodel_uspec_full_mlvl, data=c(p_pspec=zero, data_all), verbose=T, iter=n_iterations, 
                                 chains=n_chains, cores=n_chains, init=init_nspec_list, seed=1234)
pspec_full_fit <- sampling(uebermodel_uspec_full_mlvl, data=c(p_pspec=one, data_all), verbose=T, iter=n_iterations, 
                                 chains=n_chains, cores=n_chains, init=init_pspec_list, seed=1234)
save(nspec_full_fit, pspec_full_fit, file="./fit_uspec_full_mlvl.rda")

#######################################################################################################################################
#######################################################################################################################################
#######################################################################################################################################

if (FALSE)
{
show_params(nspec_full_fit)

(extract_log_lik(nouspec_fit_full) %>% waic)
(extract_log_lik(nspec_amb_fit) %>% waic)
(extract_log_lik(nspec_amb2groups_fit) %>% waic)
(extract_log_lik(nspec_full_fit) %>% waic)
show_params(nspec_amb_fit)
show_params(pspec_amb_fit)



plot_predictions(d, pspec_amb_fit)
(plot_predictions(d, nspec_amb_fit))




uebermodel_uspec_amb2g_mlvl <- stan_model(model_code = read_model("./models/uebermodel_functions.stan", "./models/uebermodel_mlvl_uspec_amb_2groups.stan"))

nspec_amb_2groups_fit <- sampling(uebermodel_uspec_amb_mlvl, data=c(p_pspec=zero, data_all, par_fixed_uspec), verbose=T, iter=n_iterations, 
                          chains=n_chains, cores=n_chains, init=init_list, seed=1234)
pspec_amb_2groups_fit <- sampling(uebermodel_uspec_amb_mlvl, data=c(p_pspec=one, data_all, par_fixed_uspec), verbose=T, iter=n_iterations, 
                                  chains=n_chains, cores=n_chains, init=init_list, seed=1234)


(extract_log_lik(nouspec_fit_full) %>% waic)
(extract_log_lik(nspec_amb_fit) %>% waic)
(extract_log_lik(pspec_amb_fit) %>% waic)
show_params(nspec_amb_fit)
show_params(pspec_amb_fit)


# underspecification models with underspecification in unambiguous conditions as well
uebermodel_uspec_unamb_mlvl <- stan_model(model_code = read_file("./models/uebermodel_mlvl_uspec_unamb.stan"))
par_fixed_uspec_unamb <- list(p_arbitrary_attachment_subj=rep(zero,n_subj), p_arbitrary_attachment_late_subj=rep(zero,n_subj))

pspec_unamb_fit <- sampling(uebermodel_uspec_unamb_mlvl, data=c(p_pspec=one, data_all, par_fixed_uspec_unamb), verbose=T, iter=n_iterations, 
                            chains=n_chains, cores=n_chains, init=init_list, seed=1234)
nspec_unamb_fit <- sampling(uebermodel_uspec_unamb_mlvl, data=c(p_pspec=zero, data_all, par_fixed_uspec_unamb), verbose=T, iter=n_iterations,
                            chains=n_chains, cores=n_chains, init=init_list, seed=1234)
save(pspec_unamb_fit, nspec_unamb_fit, file="./fit_uspec_unamb_mlvl.rda")

# underspecification models with arbitrary attachment in unambiguous conditions
uebermodel_uspec_arb_mlvl <- stan_model(model_code = read_file("./models/uebermodel_mlvl_uspec_arb.stan"))
par_fixed_uspec_arb <- list(p_uspec_unamb_subj=rep(zero,n_subj),
                            p_arbitrary_attachment_subj=rep(one,n_subj))

pspec_arb_fit <- sampling(uebermodel_uspec_arb_mlvl, data=c(p_pspec=one, data_all, par_fixed_uspec_arb), verbose=T, iter=n_iterations, 
                            chains=n_chains, cores=n_chains, init=init_list, seed=1234)
nspec_arb_fit <- sampling(uebermodel_uspec_arb_mlvl, data=c(p_pspec=zero, data_all, par_fixed_uspec_arb), verbose=T, iter=n_iterations,
                            chains=n_chains, cores=n_chains, init=init_list, seed=1234)
save(pspec_arb_fit, nspec_arb_fit, file="./fit_uspec_arb_mlvl.rda")


# underspecification models with arbitrary attachment and underspecification in unambiguous conditions
uebermodel_uspec_unamb_arb_mlvl <- stan_model(model_code = read_file("./models/uebermodel_mlvl_uspec_unamb_arb.stan"))
par_fixed_uspec_unamb_arb <- list(p_arbitrary_attachment_subj=rep(one,n_subj))

pspec_unamb_arb_fit <- sampling(uebermodel_uspec_unamb_arb_mlvl, data=c(p_pspec=one, data_all, par_fixed_uspec_unamb_arb), verbose=T, iter=n_iterations, 
                                chains=n_chains, cores=n_chains, init=init_list, seed=1234)
nspec_unamb_arb_fit <- sampling(uebermodel_uspec_unamb_arb_mlvl, data=c(p_pspec=zero, data_all, par_fixed_uspec_unamb_arb), verbose=T, iter=n_iterations,
                                chains=n_chains, cores=n_chains, init=init_list, seed=1234)
save(pspec_unamb_arb_fit, nspec_unamb_arb_fit, file="./fit_uspec_unamb_arb_mlvl.rda")

# underspecification models with arbitrary attachment to varying degrees for early and postponed attachment and underspecification in unambiguous conditions
uebermodel_full <- stan_model(model_code = read_file("./models/uebermodel_mlvl_full.stan"))

pspec_full <- sampling(uebermodel_full, data=c(p_pspec=one, data_all), verbose=T, iter=n_iterations, 
                                chains=n_chains, cores=n_chains, init=init_list, seed=1234)
nspec_full <- sampling(uebermodel_full, data=c(p_pspec=zero, data_all), verbose=T, iter=n_iterations,
                                chains=n_chains, cores=n_chains, init=init_list, seed=1234)
save(pspec_full, nspec_full, file="./fit_uspec_full.rda")



(extract_log_lik(nouspec_fit) %>% waic)

(extract_log_lik(pspec_fit) %>% waic)
(extract_log_lik(nspec_fit) %>% waic)

(extract_log_lik(pspec_unamb_fit) %>% waic)
(extract_log_lik(nspec_unamb_fit) %>% waic)

(extract_log_lik(pspec_arb_fit) %>% waic) #!
(extract_log_lik(nspec_arb_fit) %>% waic) #!

(extract_log_lik(pspec_unamb_arb_fit) %>% waic)
(extract_log_lik(nspec_unamb_arb_fit) %>% waic)

(extract_log_lik(pspec_full) %>% waic)
(extract_log_lik(nspec_full) %>% waic)

summary(pspec_arb_fit)$summary %>% .[grepl("hyper", rownames(.)),] %>% round(2)
summary(nspec_arb_fit)$summary %>% .[grepl("hyper", rownames(.)),] %>% round(2)

summary(pspec_arb_fit)$summary %>% .[grepl("p_uspec_subj\\[", rownames(.)),] %>% round(2)
 
source("./optimize.R")



nspec_fit <- sampling(uebermodel_mlvl, data=c(data_all, p_pspec=.00001), verbose=T, chains=1, iter=50, cores=1, seed=1234)
pspec_fit <- sampling(uebermodel_mlvl, data=c(data_all, p_pspec=1-.00001), verbose=T, chains=1, iter=50, cores=1, seed=1234)

summary(nspec_fit)$summary["lp__",] %>% round(2)
summary(pspec_fit)$summary["lp__",] %>% round(2)


summary(nspec_fit)$summary[1:507,] %>% round(2)
summary(pspec_fit)$summary[1:507,] %>% round(2)



library(doMC)
doMC::registerDoMC(8)

data_all <- as.stan.data.list(d)
nspec_fit <- sampling(uebermodel_mlvl, data=c(data_all, p_pspec=.00001), verbose=T, chains=6, iter=2000, cores=6, seed=1234)
pspec_fit <- sampling(uebermodel_mlvl, data=c(data_all, p_pspec=1-.00001), verbose=T, chains=6, iter=2000, cores=6, seed=1234)

save(pspec_fit, nspec_fit, file="./uebermodel_full_.rda")

(pspec_waic <- extract_log_lik(pspec_fit) %>% waic)
(nspec_waic <- extract_log_lik(nspec_fit) %>% waic)

compare(pspec_waic, nspec_waic)

summary(nspec_fit)$summary[1:400,]


extract_log_lik(pspec_fit) %>% waic
extract_log_lik(nspec_fit) %>% waic

summary(pspec_fit)

pspec_log_lik <- extract_log_lik(pspec_fit)
waic(pspec_log_lik)

init <- list(init_nspec, init_nspec, init_nspec, init_nspec)
nspec_fit <- sampling(nspec_simple_all, data=as.stan.data.list(d), verbose=T, chains=4, iter=150, cores=4, init=init)



library(doMC)
library(plyr)
doMC::registerDoMC(8)

all_data_lst <- as.stan.data.list(d)
models <- list(nspec_simple_all, nspec_simple_all, nspec_simple_all, nspec_simple_all,
               pspec_simple_all, pspec_simple_all, pspec_simple_all, pspec_simple_all)
sflist <- llply(1:8, function(i) sampling(models[[i]], data=all_data_lst, verbose=T, chains=1, iter=2000,
                                          sample_file=sprintf("./uebermodel_uspec_multi.samples_%d", i),
                                          diagnostic_file=sprintf("./uebermodel_uspec_multi.diagnostic_%d", i)
), .parallel=T)
save(sflist, file="models_simple_all.rda")

sflist[[4]]

nspec_fit <- sflist2stanfit(sflist[c(1,2,4)])
pspec_fit <- sflist2stanfit(sflist[c(6,7,8)])

nspec_fit

sampling(nspec_simple_all, data=as.stan.data.list(d), verbose=T, chains=1, iter=100)

init=list(init), 

# 100 iterations: 
# - 84.915 sec, no hyperparams
# - ~75 sec, hyperparams: beta1
# - 80 sec, hyperparams: beta0,1
# - 100 sec, hyperparams: alpha0,1,2, beta0,1
# - 133 sec, hyperparams: scale, alpha0,1,2, beta0,1




pspec_single_fits <- 
  dlply(d, .(subj), function(d) {
    cur_data <- as.stan.data.list(d)
    sampling(pspec_simple, data=cur_data, verbose=T, chains=4, iter=2000)
  }, .parallel = T)

# collect all the by-participant estimates
pspec_init <- sapply(pspec_single_fits, function(fit) { 
  summary(fit)$summary[,'50%']
})

pspec_init <- as.list(as.data.frame(t(pspec_init)))
pspec_init$lp__ <- NULL





lapply(list, function)
apply(pspec_init, 1, function(x) list(x))



pspec_lp <- sapply(pspec_single_fits, function(fit) { 
  summary(fit)$summary['lp__', 'mean']
})
nspec_lp <- sapply(nspec_single_fits, function(fit) { 
  summary(fit)$summary['lp__', 'mean']
})

summary(pspec_single_fits[[37]])$summary['lp__', 'mean']

pspec_single_fits[[37]]
pspec_lp
nspec_lp

plot(pspec_lp - nspec_lp)

ddply(d, .(subject), function(d) {
  sflist <- mclapply(1:4, mc.cores = 4, 
                     function(i) sampling(pspec_simple, data=d_lst, verbose=T, chains=1, iter=500
                     ))
  fit <- sflist2stanfit(sflist)
})



library(parallel)
rng_seed <- "964";
ddply(d, .(subject), function(d) {
  sflist <- mclapply(1:4, mc.cores = 4, 
                     function(i) sampling(pspec_simple, data=d_lst, verbose=T, chains=1, iter=500
                     ))
  fit <- sflist2stanfit(sflist)
})

fit_tst <- sampling(pspec_simple, data=data_all, verbose=T, chains=2, iter=50)
fit_tst




source("./optimize.R")
pspec_mlvl <- stan_model(model_name="uebermodel_multi_pspec", model_code = read_file("./mvm_pspec_mlvl.stan"))
nspec_mlvl <- stan_model(model_name="uebermodel_multi_pspec", model_code = read_file("./mvm_nspec_mlvl.stan"))

#parnames_fixeff <- c("p_uspec", "p_att_n1", "p_retrieval_fail", "p_guess_yes", "alpha0", "alpha1", "alpha2", "beta0", "beta1", "scale")
#init_fixeff <- rnorm(length(parnames_fixeff))
#names(init_fixeff) <- parnames_fixeff

transformed_vars <- c() # c("scale_subj", "alpha0_subj", "alpha1_subj", "alpha2_subj", "beta0_subj", "beta1_subj")
# init_prob <- "runif(0.01,0.99)"
init_fixeff <- list(p_uspec="runif(.1,.3)", p_att_n1="runif(.2,.5)", p_retrieval_fail="runif(.4,.6)", p_guess_yes=.5, 
                    alpha0="runif(3, 6)", alpha1=1, alpha2=1,
                    beta0="runif(4, 7)", beta1="runif(3, 7)", 
                    scale="runif(200, 400)")
init_raneff_sd <- list(p_uspec_subj_sd=0.1, p_att_n1_subj_sd=0.1, p_retrieval_fail_subj_sd=0.1, p_guess_yes_subj_sd=0.1,
                       scale_subj_sd=0.5, alpha0_subj_sd=0.5, alpha1_subj_sd=0.1, alpha2_subj_sd=0.1, beta0_subj_sd=0.5, beta1_subj_sd=0.1
)
init_raneff <- list(p_uspec_subj=0, p_att_n1_subj=0, p_retrieval_fail_subj=0, p_guess_yes_subj=0, 
                    scale_subj=0, alpha0_subj=0, alpha1_subj=0, 
                    alpha2_subj=0, beta0_subj=0, beta1_subj=0)

data_all <- as.stan.data.list(d)
init <- do.call(c, list(init_fixeff, init_raneff_prob, init_raneff_sd, init_raneff))
op1 <- optimizer(model_functions(pspec_mlvl, data=data_all, transformed_vars), parnames_free=names(init), init=init)


fit_tst <- sampling(pspec_mlvl, data=data_all, verbose=T, init=list(op1$par_constr(), op1$par_constr()),
                    chains=2, iter=50)
fit_tst


# optimize over all data, just fixed effects
data_all <- as.stan.data.list(d)
init <- do.call(c, list(init_fixeff, init_raneff_prob, init_raneff_sd, init_raneff))
op1 <- optimizer(model_functions(pspec_mlvl, data=data_all, transformed_vars), parnames_free=names(init), init=init)

(res1 <- op1$optimize(method="L-BFGS-B", maxit=10^5, verbose=TRUE))
(res1 <- op1$optimize(method="BFGS", maxit=10^5, verbose=TRUE))
par1 <- op1$par_constr(transformed = TRUE)

(res1b <- op1$optimize(method="Nelder-Mead", maxit=10^6, verbose=TRUE))
par1b <- op1$par_constr(transformed = TRUE)

# res1$value
# res1b$value
# #[1] -19651.29
# 
# 
# cbind(unlist(par1), unlist(par1b), unlist(par1b)-unlist(par1)) %>% round(2)
# 
# 
# save(res1, init, op1, file="./tmp.rda")
# 
# op2 <- optimizer(model_functions(nspec_mlvl, data=data_all, transformed_vars), parnames_free=names(init), init=init)
# (res2 <- op2$optimize(method="L-BFGS-B", maxit=10^5))
# (res2 <- op2$optimize(method="Nelder-Mead", maxit=10^8))
# 
# op1b <- optimizer(model_functions(pspec_mlvl, data=data_all, transformed_vars), parnames_free=names(init), init=op2$par_constr())
# (res1b <- op1b$optimize(method="L-BFGS-B", maxit=10^5))
# 
# op1$par_constr()[parnames_fixeff] %>% unlist
# op1b$par_constr()[parnames_fixeff] %>% unlist
# 
# sort(log( res1$par[parnames_fixeff] / res1b$par[parnames_fixeff] ))




# library(ggplot2)
# cur_subj <- 22
# pspec_single_fits[[cur_subj]]
# cur_data <- subset(d.qrc.stan, subj==cur_subj)
# cur_data$ypos <- rnorm(nrow(cur_data), sd=10^-5)
# cur_data$response_correct <- with(cur_data, ifelse(iv_cond==0, "NA", ifelse(iv_cond==1, ifelse(iv_questN1==response_yes, "corr", "wrong"), ifelse(iv_questN1!=response_yes, "corr", "wrong")  ))) %>% as.factor
# medians <- summary(pspec_single_fits[[cur_subj]])$summary[,'50%']
# 
# head(cur_data)
# 
# p <- ggplot(cur_data, aes(x=response_RT)) + geom_density(linetype="dotted") + 
#   geom_point(size=3, alpha=.5, aes(y=ypos, color=as.factor(iv_cond), shape=response_correct))
# p <- p + theme_bw() + scale_x_continuous(limits=c(0, max(cur_data$response_RT)))
# p <- p + #stat_function(fun=function(x) dgamma(x, shape=medians['beta0'], scale=medians['scale'])) + 
#          stat_function(fun=function(x) dgamma(x, shape=sum(medians[c('beta0','alpha1','alpha2')]), scale=medians['scale']), color="green") +
#          stat_function(fun=function(x) dgamma(x, shape=sum(medians[c('beta0','alpha1')]), scale=medians['scale']), color="blue") +
#          stat_function(fun=function(x) dgamma(x, shape=sum(medians[c('beta0','beta1')]), scale=medians['scale']))
# p
}
