library(plyr)
library(dplyr)

par2list <- function(par) {
  names <- gsub("\\[\\d+\\]$", "", names(par)) %>% gsub("\\.\\d+\\.$", "", .)
  df <- data.frame(name=names, val=unlist(par))
  dlply(df, .(name), function(df) {
    df$val
  })
}


model_functions <- function(model, data, transformed_vars) {
  model_fit <- sampling(model, data=data, chains=0, iter=1)
  n_args <- model_fit@.MISC$stan_fit_instance$num_pars_unconstrained() #do.call(sum, model_fit@par_dims)
  names_args <- model_fit@model_pars %>% setdiff("lp__")
  fn <- function(par) rstan::log_prob(model_fit, par)
  gr <- function(par) rstan::grad_log_prob(model_fit, par)
  constrain <- function(par, transformed=FALSE) {
    pars <- constrain_pars(model_fit, par)
    if (!transformed) {
      pars <- pars[setdiff(names_args, transformed_vars)]
    }
    pars
  }
  unconstrain <- function(par) unconstrain_pars(model_fit, par)
  list(fn=fn, gr=gr, n_args=n_args, names_args=names_args, constrain=constrain, unconstrain=unconstrain, mfit=model_fit)
}


smooth_extreme_values <- function(values, as_inf=10^4) {
  values[values == -Inf] <- -as_inf
  values[values == Inf] <- as_inf
  values
}


params_structure <- function(model_fn) {
  parval_unconstr <- rep(NA, model_fn$n_args)
  parval <- model_fn$constrain(parval_unconstr)
  if (model_fn$n_args != length(unlist(parval)) ) {
    stop("Missed some transformed variables.")
  }
  parval
}

params_initialize <- function(parval, init)
{
  for (name in names(parval))
  {
    n_vals <- length(parval[[name]]) 
    if(!name %in% names(init)) {
      stop(sprintf("%s is not in init.", name))
    }
    vals <- init[[name]]
    suppressWarnings(is_numeric <- all(!is.na(as.numeric( vals ))))
    if (is_numeric) {
      if (length(vals) < n_vals ) {
        vals <- as.numeric(vals) %>% rep(n_vals)
      } else {
        vals <- as.numeric(vals)
      }
    } else {
      vals <- eval(parse(text=paste("length(parval[[name]]) %>%", init[[name]])))
      stopifnot(length(vals) == n_vals)
    }
    parval[[name]] <- vals
  }
  parval
}

indicator_param_free <- function(parstruct, free) {
  ind_free <- c()
  for (name in names(parstruct)) {
    n_vals <- length(parstruct[[name]]) 
    ind_free <- c(ind_free, rep((name %in% free), n_vals))
  }
  ind_free
}


optimizer <- function(model_fn, parnames_free, init=c()) 
{
  parval_constr <- params_structure(model_fn)
  parval_constr <- params_initialize(parval_constr, init)
  parval_unconstr <- model_fn$unconstrain(parval_constr) %>% smooth_extreme_values()
  ind_parfree <- indicator_param_free(parval_constr, parnames_free)

  free_parnames <- function(names) {
    parnames_free <<- unique(c(parnames_free, names))
    ind_parfree <<- indicator_param_free(parval_constr, parnames_free)  
  }
  fix_parnames <- function(names) {
    parnames_free <<- setdiff(parnames_free, names)
    ind_parfree <<- indicator_param_free(parval_constr, parnames_free)  
  }

  successively_free_parnames <- function(names_list, method, maxit, verbose=TRUE) {
    if(!verbose) cat <- function(...) {}
    for(names in names_list) {
      cat("Freeing ", paste(names, collapse=","))
      free_parnames(names)
      res <- optimize(method=method, maxit=maxit)
      cat("\t\tValue: ", res$value, "\n")
    }
    cat("---\n")
    res
  }
  
  n_free_par <- function() sum(ind_parfree)
  
  expand_par <- function(par) {
    stopifnot(length(par) == n_free_par())
    parval_unconstr[ind_parfree] <- par
    parval_unconstr
  }

  update_parval <- function(par) {
    parval_unconstr <<- expand_par(par)
    parval_constr <<- model_fn$constrain(parval_unconstr)
  }

  optimize <- function(method, maxit, verbose=FALSE) {
    init <- get_freeval()
    res <- optim(init, fn = restricted_fn, gr = restricted_gr,
                 verbose=verbose,
                 method=method,
                 control = list(fnscale=-1, maxit=maxit))
    update_parval(res$par)
    res
  }
  
  get_freeval <- function() {
    init_unconstr <- parval_unconstr[ind_parfree]
    names(init_unconstr) <- names(unlist(parval_constr))[ind_parfree]
    init_unconstr
  }
  
  restricted_fn <- function(par, verbose=FALSE) {
    #print(model_fn$constrain(expand_par(par)))
    res <- model_fn$fn(expand_par(par))
    if (verbose)
      cat(res, "\n")
    res
  }
  restricted_gr <- function(par, verbose=FALSE) {
    res <- model_fn$gr(expand_par(par))[ind_parfree]
    #print(res)
    res
  }

  list(get_freeval = get_freeval,
       par_constr = function(transformed=FALSE) model_fn$constrain(parval_unconstr, transformed),
       par_unconstr = function() parval_unconstr, 
       parnames_free = function() parnames_free,
       free_parnames = free_parnames,
       fix_parnames = fix_parnames,
       successively_free_parnames = successively_free_parnames,
       fn = restricted_fn, 
       gr = restricted_gr,
       optimize = optimize)  
}


run_optimizing <- function(model, data, n, n_restart=1, verbose=FALSE, init_range=c(-10,10), maxit=3000, ...) {
  
}

rerun_optimizing <- function(model, data, n, n_restart=1, verbose=FALSE, init_range=c(-10,10), maxit=3000, ...)
{
  model_fn <- model_functions(model, data=data)
  
  model_fn$constrain(rep(0, model_fn$n_args))
  
  model_fit@.MISC$stan_fit_instance$.pointer
  
  
  init <- runif(n_args, min=init_range[1], max=init_range[2])
  
  p <- progress_text()
  p$init(n)
  res_df <- data.frame()
  for( i in 1:n ) {
    # 
    res <- rstan::optimizing(model, data=data, verbose=verbose, ...)
    for(i in 1:n_restart) {
      res <- rstan::optimizing(model, init=par2list(res$par), data=data, verbose=verbose,  ...)
    }
    #res <- rstan::optimizing(model, data=data, verbose=verbose, algorithm="BFGS")
    #res <- optim(init, fn=fn, gr=gr, method="L-BFGS-B", control = list(fnscale=-1, maxit=maxit))
    #res
    res <- data.frame(value=res$value, t(res$par))
    res_df <- rbind(res_df, res)
    p$step()
  }
  p$term()
  res_df %>% arrange(desc(value))
}

