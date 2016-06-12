functions
{
	real beta_alpha(real mu, real varprop) {
		real sigmasq;
		sigmasq <- varprop*mu*(1-mu);
		return mu*(mu*(1-mu)/sigmasq - 1);
	}

	real beta_beta(real mu, real varprop) {
		real sigmasq;
		sigmasq <- varprop*mu*(1-mu);
		return (1-mu)*(mu*(1-mu)/sigmasq - 1);
	}
  
	real beta_log_hyper_mix(real value, real[] par_hyper) {
	  real log_lik[2];
	  real alpha;
	  real beta;

	  alpha <- beta_alpha(par_hyper[1], par_hyper[2]);
	  beta <- beta_beta(par_hyper[1], par_hyper[2]);
		log_lik[1] <- log(par_hyper[5]) + beta_log(value, alpha, beta);
		
	  alpha <- beta_alpha(par_hyper[1]*par_hyper[3], par_hyper[4]);
	  beta <- beta_beta(par_hyper[1]*par_hyper[3], par_hyper[4]);
		log_lik[2] <- log(1-par_hyper[5]) + beta_log(value, alpha, beta);
		
		return log_sum_exp(log_lik);
	}

	real beta_log_hyper(real value, real[] par_hyper) {
	  real alpha;
	  real beta;
	  alpha <- beta_alpha(par_hyper[1], par_hyper[2]);
	  beta <- beta_beta(par_hyper[1], par_hyper[2]);
		return beta_log(value, alpha, beta);
	}

	real model_log_lik(int iv_cond, int crit_region_cnt, int iv_questN1, int response_yes, real reading_time, real response_RT,
					   real p_pspec, real p_uspec, real p_uspec_unamb, real p_arbitrary_attachment,
					   real p_att_n1, real p_retrieval_fail, real p_guess_yes,
					   real scale, real shape_reading_base, real shape_delta_attach_N2, real shape_delta_attach_N1, 
					   real shape_response_base, real shape_response_delta_guess) 
	{
		int i;
		real prob_trial_type;
		real logProb_path;
		real shape_reading;
		real shape_response;
		real pY;
		real final_logLik;
		
    if(iv_cond == 0)
    {
      real logLik[7];
      
      // path: start0-start-cAMB-N1-GUESS;
      shape_reading <- shape_reading_base * crit_region_cnt;
      shape_response <- shape_response_base;
      shape_reading <- shape_reading + shape_delta_attach_N2 + shape_delta_attach_N1;
      pY <- p_guess_yes;
      shape_response <- shape_response + shape_response_delta_guess;
      logProb_path <- log(((1 - p_uspec) * p_att_n1) * (p_retrieval_fail));
      logLik[1] <- logProb_path + bernoulli_log(response_yes, pY) + gamma_log(reading_time, shape_reading, 1/scale) + gamma_log(response_RT, shape_response, 1/scale);
      
      // path: start0-start-cAMB-N1-resp_N1;
      shape_reading <- shape_reading_base * crit_region_cnt;
      shape_response <- shape_response_base;
      shape_reading <- shape_reading + shape_delta_attach_N2 + shape_delta_attach_N1;
      pY <- (iv_questN1 == 1);
      logProb_path <- log(((1 - p_uspec) * p_att_n1) * (1 - p_retrieval_fail));
      logLik[2] <- logProb_path + bernoulli_log(response_yes, pY) + gamma_log(reading_time, shape_reading, 1/scale) + gamma_log(response_RT, shape_response, 1/scale);
      
      // path: start0-start-cAMB-N2-GUESS;
      shape_reading <- shape_reading_base * crit_region_cnt;
      shape_response <- shape_response_base;
      shape_reading <- shape_reading + shape_delta_attach_N2;
      pY <- p_guess_yes;
      shape_response <- shape_response + shape_response_delta_guess;
      logProb_path <- log(((1 - p_uspec) * (1 - p_att_n1)) * (p_retrieval_fail));
      logLik[3] <- logProb_path + bernoulli_log(response_yes, pY) + gamma_log(reading_time, shape_reading, 1/scale) + gamma_log(response_RT, shape_response, 1/scale);
      
      // path: start0-start-cAMB-N2-resp_N2;
      shape_reading <- shape_reading_base * crit_region_cnt;
      shape_response <- shape_response_base;
      shape_reading <- shape_reading + shape_delta_attach_N2;
      pY <- (iv_questN1 == 0);
      logProb_path <- log(((1 - p_uspec) * (1 - p_att_n1)) * (1 - p_retrieval_fail));
      logLik[4] <- logProb_path + bernoulli_log(response_yes, pY) + gamma_log(reading_time, shape_reading, 1/scale) + gamma_log(response_RT, shape_response, 1/scale);
      
      // path: start0-start-cAMB-USPEC_AMB-GUESS;
      shape_reading <- shape_reading_base * crit_region_cnt;
      shape_response <- shape_response_base;
      pY <- p_guess_yes;
      shape_response <- shape_response + shape_response_delta_guess;
      logProb_path <- log((p_uspec) * (p_pspec * p_retrieval_fail + (1 - p_pspec)));
      logLik[5] <- logProb_path + bernoulli_log(response_yes, pY) + gamma_log(reading_time, shape_reading, 1/scale) + gamma_log(response_RT, shape_response, 1/scale);
      
      // path: start0-start-cAMB-USPEC_AMB-late_N1-resp_N1;
      shape_reading <- shape_reading_base * crit_region_cnt;
      shape_response <- shape_response_base;
      shape_response <- shape_response + shape_delta_attach_N2 + shape_delta_attach_N1;
      pY <- (iv_questN1 == 1);
      logProb_path <- log((p_uspec) * (p_pspec * (1 - p_retrieval_fail) * p_att_n1));
      logLik[6] <- logProb_path + bernoulli_log(response_yes, pY) + gamma_log(reading_time, shape_reading, 1/scale) + gamma_log(response_RT, shape_response, 1/scale);
      
      // path: start0-start-cAMB-USPEC_AMB-late_N2-resp_N2;
      shape_reading <- shape_reading_base * crit_region_cnt;
      shape_response <- shape_response_base;
      shape_response <- shape_response + shape_delta_attach_N2;
      pY <- (iv_questN1 == 0);
      logProb_path <- log((p_uspec) * (p_pspec * (1 - p_retrieval_fail) * (1 - p_att_n1)));
      logLik[7] <- logProb_path + bernoulli_log(response_yes, pY) + gamma_log(reading_time, shape_reading, 1/scale) + gamma_log(response_RT, shape_response, 1/scale);
      
      final_logLik <- log_sum_exp(logLik);
    };
    
    if(iv_cond == 1)
    {
      real logLik[6];
      
      // path: start0-start-cN1-N1-GUESS;
      shape_reading <- shape_reading_base * crit_region_cnt;
      shape_response <- shape_response_base;
      shape_reading <- shape_reading + shape_delta_attach_N2 + shape_delta_attach_N1;
      pY <- p_guess_yes;
      shape_response <- shape_response + shape_response_delta_guess;
      logProb_path <- log(((1 - p_uspec_unamb) * ((1 - p_arbitrary_attachment) + p_arbitrary_attachment * p_att_n1)) * (p_retrieval_fail));
      logLik[1] <- logProb_path + bernoulli_log(response_yes, pY) + gamma_log(reading_time, shape_reading, 1/scale) + gamma_log(response_RT, shape_response, 1/scale);
      
      // path: start0-start-cN1-N1-resp_N1;
      shape_reading <- shape_reading_base * crit_region_cnt;
      shape_response <- shape_response_base;
      shape_reading <- shape_reading + shape_delta_attach_N2 + shape_delta_attach_N1;
      pY <- (iv_questN1 == 1);
      logProb_path <- log(((1 - p_uspec_unamb) * ((1 - p_arbitrary_attachment) + p_arbitrary_attachment * p_att_n1)) * (1 - p_retrieval_fail));
      logLik[2] <- logProb_path + bernoulli_log(response_yes, pY) + gamma_log(reading_time, shape_reading, 1/scale) + gamma_log(response_RT, shape_response, 1/scale);
      
      // path: start0-start-cN1-N2-GUESS;
      shape_reading <- shape_reading_base * crit_region_cnt;
      shape_response <- shape_response_base;
      shape_reading <- shape_reading + shape_delta_attach_N2;
      pY <- p_guess_yes;
      shape_response <- shape_response + shape_response_delta_guess;
      logProb_path <- log(((1 - p_uspec_unamb) * p_arbitrary_attachment * (1 - p_att_n1)) * (p_retrieval_fail));
      logLik[3] <- logProb_path + bernoulli_log(response_yes, pY) + gamma_log(reading_time, shape_reading, 1/scale) + gamma_log(response_RT, shape_response, 1/scale);
      
      // path: start0-start-cN1-N2-resp_N2;
      shape_reading <- shape_reading_base * crit_region_cnt;
      shape_response <- shape_response_base;
      shape_reading <- shape_reading + shape_delta_attach_N2;
      pY <- (iv_questN1 == 0);
      logProb_path <- log(((1 - p_uspec_unamb) * p_arbitrary_attachment * (1 - p_att_n1)) * (1 - p_retrieval_fail));
      logLik[4] <- logProb_path + bernoulli_log(response_yes, pY) + gamma_log(reading_time, shape_reading, 1/scale) + gamma_log(response_RT, shape_response, 1/scale);
      
      // path: start0-start-cN1-USPEC_N1-GUESS;
      shape_reading <- shape_reading_base * crit_region_cnt;
      shape_response <- shape_response_base;
      pY <- p_guess_yes;
      shape_response <- shape_response + shape_response_delta_guess;
      logProb_path <- log((p_uspec_unamb) * (p_pspec * p_retrieval_fail + (1 - p_pspec)));
      logLik[5] <- logProb_path + bernoulli_log(response_yes, pY) + gamma_log(reading_time, shape_reading, 1/scale) + gamma_log(response_RT, shape_response, 1/scale);
      
      // path: start0-start-cN1-USPEC_N1-late_N1-resp_N1;
      shape_reading <- shape_reading_base * crit_region_cnt;
      shape_response <- shape_response_base;
      shape_response <- shape_response + shape_delta_attach_N2 + shape_delta_attach_N1;
      pY <- (iv_questN1 == 1);
      logProb_path <- log((p_uspec_unamb) * (p_pspec * (1 - p_retrieval_fail)));
      logLik[6] <- logProb_path + bernoulli_log(response_yes, pY) + gamma_log(reading_time, shape_reading, 1/scale) + gamma_log(response_RT, shape_response, 1/scale);
      
      final_logLik <- log_sum_exp(logLik);
    };
    
    if(iv_cond == 2)
    {
      real logLik[6];
      
      // path: start0-start-cN2-N1-GUESS;
      shape_reading <- shape_reading_base * crit_region_cnt;
      shape_response <- shape_response_base;
      shape_reading <- shape_reading + shape_delta_attach_N2 + shape_delta_attach_N1;
      pY <- p_guess_yes;
      shape_response <- shape_response + shape_response_delta_guess;
      logProb_path <- log(((1 - p_uspec_unamb) * p_arbitrary_attachment * p_att_n1) * (p_retrieval_fail));
      logLik[1] <- logProb_path + bernoulli_log(response_yes, pY) + gamma_log(reading_time, shape_reading, 1/scale) + gamma_log(response_RT, shape_response, 1/scale);
      
      // path: start0-start-cN2-N1-resp_N1;
      shape_reading <- shape_reading_base * crit_region_cnt;
      shape_response <- shape_response_base;
      shape_reading <- shape_reading + shape_delta_attach_N2 + shape_delta_attach_N1;
      pY <- (iv_questN1 == 1);
      logProb_path <- log(((1 - p_uspec_unamb) * p_arbitrary_attachment * p_att_n1) * (1 - p_retrieval_fail));
      logLik[2] <- logProb_path + bernoulli_log(response_yes, pY) + gamma_log(reading_time, shape_reading, 1/scale) + gamma_log(response_RT, shape_response, 1/scale);
      
      // path: start0-start-cN2-N2-GUESS;
      shape_reading <- shape_reading_base * crit_region_cnt;
      shape_response <- shape_response_base;
      shape_reading <- shape_reading + shape_delta_attach_N2;
      pY <- p_guess_yes;
      shape_response <- shape_response + shape_response_delta_guess;
      logProb_path <- log(((1 - p_uspec_unamb) * ((1 - p_arbitrary_attachment) + p_arbitrary_attachment * (1 - p_att_n1))) * (p_retrieval_fail));
      logLik[3] <- logProb_path + bernoulli_log(response_yes, pY) + gamma_log(reading_time, shape_reading, 1/scale) + gamma_log(response_RT, shape_response, 1/scale);
      
      // path: start0-start-cN2-N2-resp_N2;
      shape_reading <- shape_reading_base * crit_region_cnt;
      shape_response <- shape_response_base;
      shape_reading <- shape_reading + shape_delta_attach_N2;
      pY <- (iv_questN1 == 0);
      logProb_path <- log(((1 - p_uspec_unamb) * ((1 - p_arbitrary_attachment) + p_arbitrary_attachment * (1 - p_att_n1))) * (1 - p_retrieval_fail));
      logLik[4] <- logProb_path + bernoulli_log(response_yes, pY) + gamma_log(reading_time, shape_reading, 1/scale) + gamma_log(response_RT, shape_response, 1/scale);
      
      // path: start0-start-cN2-USPEC_N2-GUESS;
      shape_reading <- shape_reading_base * crit_region_cnt;
      shape_response <- shape_response_base;
      pY <- p_guess_yes;
      shape_response <- shape_response + shape_response_delta_guess;
      logProb_path <- log((p_uspec_unamb) * (p_pspec * p_retrieval_fail + (1 - p_pspec)));
      logLik[5] <- logProb_path + bernoulli_log(response_yes, pY) + gamma_log(reading_time, shape_reading, 1/scale) + gamma_log(response_RT, shape_response, 1/scale);
      
      // path: start0-start-cN2-USPEC_N2-late_N2-resp_N2;
      shape_reading <- shape_reading_base * crit_region_cnt;
      shape_response <- shape_response_base;
      shape_response <- shape_response + shape_delta_attach_N2;
      pY <- (iv_questN1 == 0);
      logProb_path <- log((p_uspec_unamb) * (p_pspec * (1 - p_retrieval_fail)));
      logLik[6] <- logProb_path + bernoulli_log(response_yes, pY) + gamma_log(reading_time, shape_reading, 1/scale) + gamma_log(response_RT, shape_response, 1/scale);
      
      final_logLik <- log_sum_exp(logLik);
    };
		return final_logLik;
	}
	
	vector model_pred_rng(int iv_cond, int crit_region_cnt, int iv_questN1, int response_yes, real reading_time, real response_RT,
					   real p_pspec, real p_uspec, real p_uspec_unamb, real p_arbitrary_attachment,
					   real p_att_n1, real p_retrieval_fail, real p_guess_yes,
					   real scale, real shape_reading_base, real shape_delta_attach_N2, real shape_delta_attach_N1, 
					   real shape_response_base, real shape_response_delta_guess) 
	{
		real shape_reading;
		real shape_response;
		real pY;
		int selected_path;
		int i;
		vector[3] predictions;

    if(iv_cond == 0)
    {
  		vector[7] prob_path;
  		
      // path: start0-start-cAMB-N1-GUESS;
      prob_path[1] <- ((1 - p_uspec) * p_att_n1) * (p_retrieval_fail);
      // path: start0-start-cAMB-N1-resp_N1;
      prob_path[2] <- ((1 - p_uspec) * p_att_n1) * (1 - p_retrieval_fail);
      // path: start0-start-cAMB-N2-GUESS;
      prob_path[3] <- ((1 - p_uspec) * (1 - p_att_n1)) * (p_retrieval_fail);
      // path: start0-start-cAMB-N2-resp_N2;
      prob_path[4] <- ((1 - p_uspec) * (1 - p_att_n1)) * (1 - p_retrieval_fail);
      // path: start0-start-cAMB-USPEC_AMB-GUESS;
      prob_path[5] <- (p_uspec) * (p_pspec * p_retrieval_fail + (1 - p_pspec));
      // path: start0-start-cAMB-USPEC_AMB-late_N1-resp_N1;
      prob_path[6] <- (p_uspec) * (p_pspec * (1 - p_retrieval_fail) * p_att_n1);
      // path: start0-start-cAMB-USPEC_AMB-late_N2-resp_N2;
      prob_path[7] <- ((p_uspec) * (p_pspec * (1 - p_retrieval_fail) * (1 - p_att_n1)));
      
      selected_path <- categorical_rng(prob_path);

      if (selected_path == 1) { // path: start0-start-cAMB-N1-GUESS;
        shape_reading <- shape_reading_base * crit_region_cnt + shape_delta_attach_N2 + shape_delta_attach_N1;
        shape_response <- shape_response_base + shape_response_delta_guess;
        pY <- p_guess_yes;

      } else if (selected_path == 2) { // path: start0-start-cAMB-N1-resp_N1;
        shape_reading <- shape_reading_base * crit_region_cnt + shape_delta_attach_N2 + shape_delta_attach_N1;
        shape_response <- shape_response_base;
        pY <- (iv_questN1 == 1);

      } else if (selected_path == 3) { // path: start0-start-cAMB-N2-GUESS;
        shape_reading <- shape_reading_base * crit_region_cnt + shape_delta_attach_N2;
        shape_response <- shape_response_base + shape_response_delta_guess;
        pY <- p_guess_yes;

      } else if (selected_path == 4) { // path: start0-start-cAMB-N2-resp_N2;
        shape_reading <- shape_reading_base * crit_region_cnt + shape_delta_attach_N2;
        shape_response <- shape_response_base;
        pY <- (iv_questN1 == 0);

      } else if (selected_path == 5) { // path: start0-start-cAMB-USPEC_AMB-GUESS;
        shape_reading <- shape_reading_base * crit_region_cnt;
        shape_response <- shape_response_base + shape_response_delta_guess;
        pY <- p_guess_yes;

      } else if (selected_path == 6) { // path: start0-start-cAMB-USPEC_AMB-late_N1-resp_N1;
        shape_reading <- shape_reading_base * crit_region_cnt;
        shape_response <- shape_response_base + shape_delta_attach_N2 + shape_delta_attach_N1;
        pY <- (iv_questN1 == 1);

      } else if (selected_path == 7) { // path: start0-start-cAMB-USPEC_AMB-late_N2-resp_N2;
        shape_reading <- shape_reading_base * crit_region_cnt;
        shape_response <- shape_response_base + shape_delta_attach_N2;
        pY <- (iv_questN1 == 0);
      }
      
    };
    
    if(iv_cond == 1)
    {
  		vector[6] prob_path;

      // path: start0-start-cN1-N1-GUESS;
      prob_path[1] <- ((1 - p_uspec_unamb) * ((1 - p_arbitrary_attachment) + p_arbitrary_attachment * p_att_n1)) * (p_retrieval_fail);
      // path: start0-start-cN1-N1-resp_N1;
      prob_path[2] <- ((1 - p_uspec_unamb) * ((1 - p_arbitrary_attachment) + p_arbitrary_attachment * p_att_n1)) * (1 - p_retrieval_fail);
      // path: start0-start-cN1-N2-GUESS;
      prob_path[3] <- ((1 - p_uspec_unamb) * p_arbitrary_attachment * (1 - p_att_n1)) * (p_retrieval_fail);
      // path: start0-start-cN1-N2-resp_N2;
      prob_path[4] <- ((1 - p_uspec_unamb) * p_arbitrary_attachment * (1 - p_att_n1)) * (1 - p_retrieval_fail);
      // path: start0-start-cN1-USPEC_N1-GUESS;
      prob_path[5] <- (p_uspec_unamb) * (p_pspec * p_retrieval_fail + (1 - p_pspec));
      // path: start0-start-cN1-USPEC_N1-late_N1-resp_N1;
      prob_path[6] <- (p_uspec_unamb) * (p_pspec * (1 - p_retrieval_fail));

      selected_path <- categorical_rng(prob_path);

      if (selected_path == 1) { // path: start0-start-cN1-N1-GUESS;
        shape_reading <- shape_reading_base * crit_region_cnt;
        shape_response <- shape_response_base;
        shape_reading <- shape_reading + shape_delta_attach_N2 + shape_delta_attach_N1;
        pY <- p_guess_yes;
        shape_response <- shape_response + shape_response_delta_guess;

      } else if (selected_path == 2) { // path: start0-start-cN1-N1-resp_N1;
        shape_reading <- shape_reading_base * crit_region_cnt;
        shape_response <- shape_response_base;
        shape_reading <- shape_reading + shape_delta_attach_N2 + shape_delta_attach_N1;
        pY <- (iv_questN1 == 1);

      } else if (selected_path == 3) { // path: start0-start-cN1-N2-GUESS;
        shape_reading <- shape_reading_base * crit_region_cnt;
        shape_response <- shape_response_base;
        shape_reading <- shape_reading + shape_delta_attach_N2;
        pY <- p_guess_yes;
        shape_response <- shape_response + shape_response_delta_guess;

      } else if (selected_path == 4) { // path: start0-start-cN1-N2-resp_N2;
        shape_reading <- shape_reading_base * crit_region_cnt;
        shape_response <- shape_response_base;
        shape_reading <- shape_reading + shape_delta_attach_N2;
        pY <- (iv_questN1 == 0);

      } else if (selected_path == 5) { // path: start0-start-cN1-USPEC_N1-GUESS;
        shape_reading <- shape_reading_base * crit_region_cnt;
        shape_response <- shape_response_base;
        pY <- p_guess_yes;
        shape_response <- shape_response + shape_response_delta_guess;

      } else if (selected_path == 6) { // path: start0-start-cN1-USPEC_N1-late_N1-resp_N1;
        shape_reading <- shape_reading_base * crit_region_cnt;
        shape_response <- shape_response_base;
        shape_response <- shape_response + shape_delta_attach_N2 + shape_delta_attach_N1;
        pY <- (iv_questN1 == 1);
      }
    };
    
    if(iv_cond == 2)
    {
      vector[6] prob_path;

      // path: start0-start-cN2-N1-GUESS;
      prob_path[1] <- (((1 - p_uspec_unamb) * p_arbitrary_attachment * p_att_n1) * (p_retrieval_fail));
      // path: start0-start-cN2-N1-resp_N1;
      prob_path[2] <- (((1 - p_uspec_unamb) * p_arbitrary_attachment * p_att_n1) * (1 - p_retrieval_fail));
      // path: start0-start-cN2-N2-GUESS;
      prob_path[3] <- (((1 - p_uspec_unamb) * ((1 - p_arbitrary_attachment) + p_arbitrary_attachment * (1 - p_att_n1))) * (p_retrieval_fail));
      // path: start0-start-cN2-N2-resp_N2;
      prob_path[4] <- (((1 - p_uspec_unamb) * ((1 - p_arbitrary_attachment) + p_arbitrary_attachment * (1 - p_att_n1))) * (1 - p_retrieval_fail));
      // path: start0-start-cN2-USPEC_N2-GUESS;
      prob_path[5] <- ((p_uspec_unamb) * (p_pspec * p_retrieval_fail + (1 - p_pspec)));
      // path: start0-start-cN2-USPEC_N2-late_N2-resp_N2;
      prob_path[6] <- ((p_uspec_unamb) * (p_pspec * (1 - p_retrieval_fail)));
      
      selected_path <- categorical_rng(prob_path);

      if (selected_path == 1) { // path: start0-start-cN2-N1-GUESS;
        shape_reading <- shape_reading_base * crit_region_cnt;
        shape_response <- shape_response_base;
        shape_reading <- shape_reading + shape_delta_attach_N2 + shape_delta_attach_N1;
        pY <- p_guess_yes;
        shape_response <- shape_response + shape_response_delta_guess;

      } else if (selected_path == 2) { // path: start0-start-cN2-N1-resp_N1;
        shape_reading <- shape_reading_base * crit_region_cnt;
        shape_response <- shape_response_base;
        shape_reading <- shape_reading + shape_delta_attach_N2 + shape_delta_attach_N1;
        pY <- (iv_questN1 == 1);

      } else if (selected_path == 3) { // path: start0-start-cN2-N2-GUESS;
        shape_reading <- shape_reading_base * crit_region_cnt;
        shape_response <- shape_response_base;
        shape_reading <- shape_reading + shape_delta_attach_N2;
        pY <- p_guess_yes;
        shape_response <- shape_response + shape_response_delta_guess;

      } else if (selected_path == 4) { // path: start0-start-cN2-N2-resp_N2;
        shape_reading <- shape_reading_base * crit_region_cnt;
        shape_response <- shape_response_base;
        shape_reading <- shape_reading + shape_delta_attach_N2;
        pY <- (iv_questN1 == 0);

      } else if (selected_path == 5) { // path: start0-start-cN2-USPEC_N2-GUESS;
        shape_reading <- shape_reading_base * crit_region_cnt;
        shape_response <- shape_response_base;
        pY <- p_guess_yes;
        shape_response <- shape_response + shape_response_delta_guess;

      } else if (selected_path == 6) { // path: start0-start-cN2-USPEC_N2-late_N2-resp_N2;
        shape_reading <- shape_reading_base * crit_region_cnt;
        shape_response <- shape_response_base;
        shape_response <- shape_response + shape_delta_attach_N2;
        pY <- (iv_questN1 == 0);
      }
    }
    
    predictions[1] <- bernoulli_rng(pY);
    predictions[2] <- gamma_rng(shape_reading, 1/scale);
    predictions[3] <- gamma_rng(shape_response, 1/scale);
    return predictions;
	}
}
