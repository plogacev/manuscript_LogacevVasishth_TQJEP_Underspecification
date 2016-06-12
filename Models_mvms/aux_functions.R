
read_model <- function (functions, rest, generated) {
  functions <- read_file(functions)
  rest <- read_file(rest)
  generated <- read_file(generated)
  paste(functions, rest, generated, sep="\n")
}

par2list <- function(par) {
  names <- gsub("\\[\\d+\\]$", "", names(par)) %>% gsub("\\.\\d+\\.$", "", .)
  df <- data.frame(name=names, val=unlist(par))
  dlply(df, .(name), function(df) {
    df$val
  })
}

generate_init <- function(model_fits, lognormal_vars, prob_vars) {
  # collect all the by-participant estimates
  init <- sapply(model_fits, function(fit) { 
    summary(fit)$summary[,'50%']
  })
  init <- as.list(as.data.frame(t(init)))
  init$lp__ <- NULL
  names(init) %<>% paste0(., "_subj")
  
  for (varname in lognormal_vars) {
    varname_subj <- paste0(varname,"_subj")
    varname_hyper <- paste0(varname,"_hyper")
    init[[varname_hyper]] <- c(median(log(init[[varname_subj]])), sd(log(init[[varname_subj]])))
  }
  for (varname in prob_vars) {
    varname_subj <- paste0(varname,"_subj")
    varname_hyper <- paste0(varname,"_hyper")
    init[[varname_hyper]] <- c(mean(init[[varname_subj]]), sd(init[[varname_subj]]))
  }
  init
}

show_params <- function(fit) summary(fit)$summary %>% .[!grepl("log_lik", rownames(.)),]%>% .[!grepl("^pred", rownames(.)),] %>% round(2)

plot_predictions <- function(d, model_fit)
{
  pred_response_yes <- extract_log_lik(model_fit, parameter_name="pred_response") %>% t %>% cbind(iv_questN1=d$iv_questN1, iv_cond=d$iv_cond) %>% 
    as.data.frame %>% gather(iteration, val, -iv_questN1, -iv_cond)
  pred_reading_time <- extract_log_lik(model_fit, parameter_name="pred_reading_time") %>% t %>% cbind(iv_questN1=d$iv_questN1, iv_cond=d$iv_cond) %>% 
    as.data.frame %>% gather(iteration, val, -iv_questN1, -iv_cond)
  pred_response_RT <- extract_log_lik(model_fit, parameter_name="pred_response_RT") %>% t %>% cbind(iv_questN1=d$iv_questN1, iv_cond=d$iv_cond) %>% 
    as.data.frame %>% gather(iteration, val, -iv_questN1, -iv_cond)
  
  pred_response_N1 <- pred_response_yes %>% group_by(iv_cond, iteration) %>% summarize(val=mean(iv_questN1 == val)) %>% 
                          summarize(lower=quantile(val, .025), upper=quantile(val, .975), val=mean(val), type="predicted", DV="% N1 Responses")
  pred_reading_time <- pred_reading_time %>% group_by(iv_cond, iteration) %>% summarize(val=mean(val)) %>% 
                          summarize(lower=quantile(val, .025), upper=quantile(val, .975), val=mean(val), type="predicted", DV="Reading Time")
  pred_response_RT <- pred_response_RT %>% group_by(iv_cond, iteration) %>% summarize(val=mean(val)) %>% 
                          summarize(lower=quantile(val, .025), upper=quantile(val, .975), val=mean(val), type="predicted", DV="Response Time")
  
  obs_response_N1 <- d %>% group_by(iv_cond) %>% summarize(lower=NA, upper=NA, val=mean(iv_questN1 == response_yes), type="observed", DV="% N1 Responses")
  obs_reading_time <- d %>% group_by(iv_cond) %>% summarize(lower=NA, upper=NA, val=mean(reading_time), type="observed", DV="Reading Time")
  obs_response_RT <- d %>% group_by(iv_cond) %>% summarize(lower=NA, upper=NA, val=mean(response_RT), type="observed", DV="Response Time")
  
  tmp <- rbind(pred_response_N1, pred_reading_time, pred_response_RT, obs_response_N1, obs_reading_time, obs_response_RT)
  tmp$DV %<>% ordered(c("% N1 Responses", "Reading Time", "Response Time"))
  tmp$iv_cond %<>% factor(labels=c("ambiguous", "N1 attachment", "N2 attachment")) %>% as.character %>% ordered(levels=c("N1 attachment", "N2 attachment", "ambiguous"))
  
  p <- tmp %>% ggplot(aes(x=iv_cond, y=val, color=type, group=type)) + facet_wrap(~DV, ncol=1, scales="free_y") + geom_point() + geom_line() + geom_errorbar(aes(ymin=lower, ymax=upper), width=.0)
  attr(p, "data") <- tmp
  p
}
