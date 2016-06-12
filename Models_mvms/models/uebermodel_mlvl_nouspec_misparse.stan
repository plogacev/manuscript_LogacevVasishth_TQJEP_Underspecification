
data
{
	int<lower=1> n_obs;
	int<lower=1> n_subj;
	int<lower=1,upper=n_subj> subj[n_obs];

	int<lower=0,upper=2> iv_cond[n_obs];
	int<lower=0,upper=1> iv_questN1[n_obs];
	int<lower=1> crit_region_cnt[n_obs];
	int<lower=0,upper=1> response_yes[n_obs];
	real<lower=0> reading_time[n_obs];
	real<lower=0> response_RT[n_obs];

	real<lower=0, upper=1> p_pspec;
	real<lower=0, upper=1> p_uspec_subj[n_subj];
	real<lower=0, upper=1> p_uspec_unamb_subj[n_subj];
	real<lower=0, upper=1> p_retrieval_fail_subj[n_subj];
}

parameters
{
	real<lower=0, upper=1> p_arbitrary_attachment_subj[n_subj];
	real<lower=0, upper=1> p_att_n1_subj[n_subj];
	real<lower=0, upper=1> p_guess_yes_subj[n_subj];
	real<lower=0, upper=2000> scale_subj[n_subj];
	real<lower=0, upper=50> shape_reading_base_subj[n_subj];
	real<lower=0, upper=50> shape_delta_attach_N2_subj[n_subj];
	real<lower=0, upper=50> shape_delta_attach_N1_subj[n_subj];
	real<lower=0, upper=50> shape_response_base_subj[n_subj];
	real<lower=0, upper=50> shape_response_delta_guess_subj[n_subj];

	//real<lower=0, upper=1> p_uspec_hyper[2];
	//real<lower=0, upper=1> p_uspec_unamb_hyper[2];
	real<lower=0, upper=1> p_arbitrary_attachment_hyper[2];
	real<lower=0, upper=1> p_att_n1_hyper[2];
	//real<lower=0, upper=1> p_retrieval_fail_hyper[2];
	real<lower=0, upper=1> p_guess_yes_hyper[2];

	real<lower=0, upper=50> scale_hyper[2];
	real<lower=0, upper=50> shape_reading_base_hyper[2];
	real<lower=0, upper=50> shape_delta_attach_N2_hyper[2];
	real<lower=0, upper=50> shape_delta_attach_N1_hyper[2];
	real<lower=0, upper=50> shape_response_base_hyper[2];
	real<lower=0, upper=50> shape_response_delta_guess_hyper[2];
}

model
{
	// specify priors for hyperparameters
	scale_hyper[1] ~ normal(5, 5);
	scale_hyper[2] ~ normal(0, 2);
	shape_reading_base_hyper[1] ~ normal(0, 5);
	shape_reading_base_hyper[2] ~ normal(0, 2);
	shape_delta_attach_N2_hyper[1] ~ normal(0, 5);
	shape_delta_attach_N2_hyper[2] ~ normal(0, 2);
	shape_delta_attach_N1_hyper[1] ~ normal(0, 5);
	shape_delta_attach_N1_hyper[2] ~ normal(0, 2);
	shape_response_base_hyper[1] ~ normal(2, 5);
	shape_response_base_hyper[2] ~ normal(0, 2);
	shape_response_delta_guess_hyper[1] ~ normal(2, 5);
	shape_response_delta_guess_hyper[2] ~ normal(0, 2);

	// specify priors for subject-level effects
	for(j in 1:n_subj) {
		//p_uspec_subj ~ beta(1, 1);
		//p_uspec_unamb_subj ~ beta(1, 1);
		p_arbitrary_attachment_subj ~ beta(1, 1);
		p_att_n1_subj ~ beta(1, 1);
		//p_retrieval_fail_subj ~ beta(1, 1);
		p_guess_yes_subj ~ beta(1, 1);
		scale_subj[j] ~ lognormal(6, 1);
		shape_reading_base_subj[j] ~ lognormal(2.5, 1);
		shape_delta_attach_N2_subj[j] ~ lognormal(2.5, 1);
		shape_delta_attach_N1_subj[j] ~ lognormal(2.5, 1);
		shape_response_base_subj[j] ~ lognormal(2.5, 1);
		shape_response_delta_guess_subj[j] ~ lognormal(2.5, 1);
	}

	// specify distribution of subject-level effects
	for(j in 1:n_subj) {
		//p_uspec_subj[j] ~ beta( beta_alpha(p_uspec_hyper[1], p_uspec_hyper[2]),
		//						            beta_beta(p_uspec_hyper[1], p_uspec_hyper[2]) );
		//p_uspec_unamb_subj[j] ~ beta( beta_alpha(p_uspec_unamb_hyper[1], p_uspec_unamb_hyper[2]),
		//						 	                beta_beta(p_uspec_unamb_hyper[1], p_uspec_unamb_hyper[2]) );
		p_arbitrary_attachment_subj[j] ~ beta( beta_alpha(p_arbitrary_attachment_hyper[1], p_arbitrary_attachment_hyper[2]),
								 	  		                    beta_beta(p_arbitrary_attachment_hyper[1], p_arbitrary_attachment_hyper[2]) );
		p_att_n1_subj[j] ~ beta( beta_alpha(p_att_n1_hyper[1], p_att_n1_hyper[2]),
								             beta_beta(p_att_n1_hyper[1], p_att_n1_hyper[2]) );
		//p_retrieval_fail_subj[j] ~ beta( beta_alpha(p_retrieval_fail_hyper[1], p_retrieval_fail_hyper[2]),
		//								                  beta_beta(p_retrieval_fail_hyper[1], p_retrieval_fail_hyper[2]) );
		p_guess_yes_subj[j] ~ beta( beta_alpha(p_guess_yes_hyper[1], p_guess_yes_hyper[2]),
									              beta_beta(p_guess_yes_hyper[1], p_guess_yes_hyper[2]) );
		scale_subj[j] ~ lognormal(scale_hyper[1], scale_hyper[2]);
		shape_reading_base_subj[j] ~ lognormal(shape_reading_base_hyper[1], shape_reading_base_hyper[2]);
		shape_delta_attach_N2_subj[j] ~ lognormal(shape_delta_attach_N2_hyper[1], shape_delta_attach_N2_hyper[2]);
		shape_delta_attach_N1_subj[j] ~ lognormal(shape_delta_attach_N1_hyper[1], shape_delta_attach_N1_hyper[2]);
		shape_response_base_subj[j] ~ lognormal(shape_response_base_hyper[1], shape_response_base_hyper[2]);
		shape_response_delta_guess_subj[j] ~ lognormal(shape_response_delta_guess_hyper[1], shape_response_delta_guess_hyper[2]);
	}

	for(i_obs in 1:n_obs)
	{
		int cur_subj;
		real logLik;
		cur_subj <- subj[i_obs];
		logLik <- model_log_lik(iv_cond[i_obs], crit_region_cnt[i_obs], iv_questN1[i_obs], 
						   response_yes[i_obs], reading_time[i_obs], response_RT[i_obs],
						   p_pspec, p_uspec_subj[cur_subj], p_uspec_subj[cur_subj]*p_uspec_unamb_subj[cur_subj], 
						   p_arbitrary_attachment_subj[cur_subj],
						   p_att_n1_subj[cur_subj], p_retrieval_fail_subj[cur_subj], p_guess_yes_subj[cur_subj],
						   scale_subj[cur_subj], shape_reading_base_subj[cur_subj], 
						   shape_delta_attach_N2_subj[cur_subj], shape_delta_attach_N1_subj[cur_subj], 
						   shape_response_base_subj[cur_subj], shape_response_delta_guess_subj[cur_subj]);
		increment_log_prob(logLik);
	};
}

