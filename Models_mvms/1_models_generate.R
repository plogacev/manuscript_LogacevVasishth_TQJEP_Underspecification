# Potentially relevant info:
# http://stackoverflow.com/questions/27202395/how-do-i-get-standard-errors-of-maximum-likelihood-estimates-in-stan

library(plyr)
library(dplyr)
library(magrittr)
library(mvmstan)
library(rstan)
rstan_options(auto_write = TRUE)

(load("./SwetsDataForModelling.rda"))
d <- d.qrc.stan
d_lst <- c(c(n_obs=nrow(d), n_subj=length(unique(d$subj)), n_item=length(unique(d$item))), as.list(d))

# read the extended models
uebermodel <- read_mvm(sprintf("./models/uebermodel.mvm", dir))

# x <- subset(d.qrc.stan, subj==3)
# x$cond <- nmap(x$iv_cond, c('0'='AMB', '1'='N1', '2'='N2'))
# x$questN1 <- nmap(x$iv_questN1, c('0'='qN2', '1'='qN1'))
# with(x, tapply(response_RT, list(questN1, response_yes,
#                                cond), function(x) paste(round(median(x)), length(x)) ))

# plot model without and without variable assignment
plot_mvm(uebermodel, show_prob = T, show_code = F, fname_dot="./models/nspec.dot", start_xdot = T)
plot_mvm(uebermodel, show_prob = T, show_code = T, fname_dot="./models/nspec.dot", start_xdot = F)

# TODO: The checks fail if a data column is included which is not used in the model. This should be a warning, not a error.
# independent variables
iv_vars = c(iv_cond='int<lower=0,upper=2>',iv_questN1='int<lower=0,upper=1>',crit_region_cnt='int<lower=1>')

# dependent variables
dv_vars = c(response_yes='int<lower=0,upper=1>', reading_time='real<lower=0>', response_RT='real<lower=0>')

# parameters to be estimated
par_vars = c(p_pspec='real<lower=0,upper=1>', p_uspec='real<lower=0, upper=1>', p_uspec_unamb='real<lower=0, upper=1>', 
                  p_arbitrary_attachment='real<lower=0, upper=1>',
                  p_retrieval_fail='real<lower=0, upper=1>', p_att_n1='real<lower=0, upper=1>',
                  p_guess_yes='real<lower=0, upper=1>',
                  scale='real<lower=0, upper=5000>', # TODO: 'scale' was missing here, and the code was still generated. Fix this behavior in mvmstan.
                  shape_reading_base='real<lower=0, upper=50>', shape_delta_attach_N2='real<lower=0, upper=50>', shape_delta_attach_N1='real<lower=0, upper=50>',
                  shape_response_base='real<lower=0, upper=50>', shape_response_delta_guess='real<lower=0, upper=50>')

# relationship between variables assigned in the model and the log-likelihood of the DVs
logLik = c(response='bernoulli_log(response_yes, pY)', 
           reading_time='gamma_log(reading_time, shape_reading, 1/scale)', 
           response_RT='gamma_log(response_RT, shape_response, 1/scale)')

mvm_generate_code(uebermodel, iv_vars=iv_vars, par_vars=par_vars, dv_vars=dv_vars, logLik=logLik, raneff=c(), 
                  file="./models/uebermodel_mvm.stan")


######################################################################
######### model with random effects ##################################
######################################################################

# now let's add random effects
raneff = c(p_uspec = 'inv_logit( logit(p_uspec) + subj)', 
           p_att_n1 = 'inv_logit( logit(p_att_n1) + subj)',
           p_retrieval_fail = 'inv_logit( logit(p_retrieval_fail) + subj)',
           p_guess_yes = 'inv_logit( logit(p_guess_yes) + subj)',
           scale  = 'exp( log(scale ) + subj)',
           alpha0 = 'exp( log(alpha0 ) + subj)',
           alpha1 = 'exp( log(alpha1 ) + subj)',
           alpha2 = 'exp( log(alpha2 ) + subj)',
           beta0  = 'exp( log(beta0 ) + subj)',
           beta1  = 'exp( log(beta1 ) + subj)'
)

logLik[['reading_time']] = 'gamma_log(reading_time, reading_shape, 1/cur_scale)'
logLik[['response_RT']] = 'gamma_log(response_RT, response_shape, 1/cur_scale)'
mvm_generate_code(m_pspec, iv_vars=iv_vars, par_vars=par_vars, dv_vars=dv_vars, logLik=logLik, raneff=raneff, file="./mvm_pspec_mlvl.stan")
mvm_generate_code(m_nspec, iv_vars=iv_vars, par_vars=par_vars, dv_vars=dv_vars, logLik=logLik, raneff=raneff, file="./mvm_nspec_mlvl.stan")

# Note: Code changed by hand
pspec_mlvl <- stan_model(model_name="uebermodel_multi_pspec", model_code = read_file("./mvm_pspec_mlvl.stan"))
nspec_mlvl <- stan_model(model_name="uebermodel_multi_pspec", model_code = read_file("./mvm_nspec_mlvl.stan"))

set.seed(1234)
res1 <- optimizing(pspec_mlvl, data=d_lst, iter=10^4, algorithm="LBFGS", verbose=TRUE)
set.seed(1234)
res2 <- optimizing(nspec_mlvl, data=d_lst, iter=10^4, algorithm="LBFGS", verbose=TRUE)

curve(dbeta(x, 20, 900), xlim=c(-.1,1.1))


res_pspec_mlvl <- run_optimizing(pspec_mlvl, data=d_lst, iter=10^5, algorithm="LBFGS", n=20)
res_nspec_mlvl <- run_optimizing(nspec_mlvl, data=d_lst, iter=10^5, algorithm="LBFGS", n=20)
res_pspec_mlvl[1:50] %>% round(2)
res_nspec_mlvl[1:30] %>% round(2)


population <- rbind(res_pspec_mlvl[,-1], res_nspec_mlvl[,-1])
population <- sapply(1:nrow(population), function(i) unconstrain_pars(model_fit, par2list(population[i,])) )
population <- t(population)
dim(population)

library(GA)
?GA::ga
fn_population <- function(...) return(as.matrix(population))
fn_population() %>% dim()

GA <- ga(type = "real-valued", fitness = fn, min = rep(-10, n_args), max = rep(10, n_args),
         population = fn_population, popSize=nrow(population))

summary(GA)
object <- GA

gaControl("real-valued")$population
gareal_Population

min <- object@min
max <- object@max
nvars <- length(min)
population <- matrix(as.double(NA), nrow = object@popSize, 
                     ncol = nvars)
for (j in 1:nvars) {
  population[, j] <- runif(object@popSize, min[j], max[j])
}
return(population)

#curve(dlnorm(x, 6, 1), xlim=c(0,4000))
#pars = c(names(res_uebermodel_multi$par)[1:11], names(res_uebermodel_multi$par[grep("_sd", names(res_uebermodel_multi$par))]))

fit_tst <- sampling(pspec_mlvl, data=d_lst, verbose=T,
                     chains=2, iter=50)
fit_tst

 
fit_tst2 <- sampling(pspec_mlvl, data=d_lst, verbose=T,
                    chains=2, iter=150)
fit_tst2


library(parallel)
rng_seed <- "964";
sflist <- mclapply(1:4, mc.cores = 4, 
                   function(i) sampling(pspec_mlvl, data=d_lst, verbose=T, #pars=pars, 
                                        chains=1, iter=500
                                        #sample_file=sprintf("./uebermodel_uspec_multi.samples_%d", i)
                                        #diagnostic_file=sprintf("./uebermodel_uspec_multi.diagnostic_%d", i)
                   ))
fit <- sflist2stanfit(sflist)
save(sflist, fit, file="./pspec_mlvl_fit_iter500.rda")
print(fit)

library(shinyStan)



(load(file="./pspec_mlvl_fit_iter1000.rda"))

library(shinystan)

shinystan::launch_shinystan(fit)


fit_tst <- sampling(uebermodel_multi_pspec, data=d_lst, verbose=T, #pars=pars, 
                    chains=2, iter=50,
                    sample_file="./uebermodel_uspec_multi.samples",
                    diagnostic_file="./uebermodel_uspec_multi.diagnostic")
fit_tst

save(fit_tst, file="./uebermodel_uspec_multi_fit_tst.rda")




fit_tst <- sampling(uebermodel_multi, data=d_lst, verbose=T, #pars=pars, 
                    chains=1, iter=1500,
                 sample_file="./uebermodel_uspec_multi.samples",
                 diagnostic_file="./uebermodel_uspec_multi.diagnostic")
save(fit_tst, file="./uebermodel_uspec_multi_fit_tst.rda")


(load(file="./uebermodel_uspec_multi_fit_tst.rda"))

fit_tst

hist( extract( fit_tst , "p_uspec" , inc_warmup=FALSE )$p_uspec )


Dbar <- mean( extract( mstan , "dev" , inc_warmup=FALSE )$dev,na.rm=TRUE )
Epars <- apply( as.data.frame( fit_tst@sim[[1]][[1]])[501:1000,] , 2 , mean )
Ebeta <- matrix( 0 , nrow=J , ncol=2 )
Ebeta[,1] <- Epars[ 8:107 ]
Ebeta[,2] <- Epars[ 108:207 ]
glm <- Ebeta[id,1] + Ebeta[id,2]*x
Dhat <- -2*sum( dnorm( y , mean=glm , sd=Epars[1] , log=TRUE ) )
(pD <- Dbar - Dhat)
(DIC <- Dbar + pD)



library(parallel)
rng_seed <- "964";
sflist <- mclapply(1:4, mc.cores = 4, 
           function(i) sampling(uebermodel_multi, data=d_lst, verbose=T, pars=pars, 
                                sample_file=sprintf("./uebermodel_uspec_multi.samples_%d", i), 
                                diagnostic_file=sprintf("./uebermodel_uspec_multi.diagnostic_%d", i)))
fit <- sflist2stanfit(sflist)
save(sflist, fit, file="./uebermodel_uspec_multi_fit.rda")


###########################################################################
# 
#  fit <- sampling(sm, data=d_lst, chains=0)
#  fn <- function(par) log_prob(fit, par)
#  gr <- function(par) grad_log_prob(fit, par)
#  
#  p <- c(0, 0, 0, 0, 0, 0)
#  fn(p)
#  gr(p)
#  
#  p2logodds <- function(p) log(p/(1-p))
#  logodds2p <- function(lodds) exp(lodds)/(1+exp(lodds))
#  
#  res <- optim(p, fn, gr, method="SANN", control=list(fnscale=-1, maxit=10^3))
#  res
#  
# logodds2p(res$par)
# 
# res2 <- optim(p, fn, gr, method="BFGS", control=list(fnscale=-1), hessian=TRUE)
# 
# logodds2p(res2$par)
# 
# sapply(seq(0,1,.05), function(p1) sapply(seq(0,1,.05), function(p2) fn(c(p1,2))))

