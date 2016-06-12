
data
{
	int<lower=1> n_obs;
	int<lower=0,upper=2> iv_cond[n_obs];
	int<lower=0,upper=1> iv_questN1[n_obs];
	int<lower=1> crit_region_cnt[n_obs];
	int<lower=0,upper=1> response_yes[n_obs];
	real<lower=0> reading_time[n_obs];
	real<lower=0> response_RT[n_obs];
}

parameters
{
	real<lower=0,upper=1> p_pspec;
	real<lower=0, upper=1> p_arbitrary_attachment;
	real<lower=0, upper=1> p_retrieval_fail;
	real<lower=0, upper=1> p_att_n1;
	real<lower=0, upper=1> p_guess_yes;
	real<lower=0, upper=2000> scale;
	real<lower=0, upper=50> shape_reading_base;
	real<lower=0, upper=50> shape_delta_attach_N2;
	real<lower=0, upper=50> shape_delta_attach_N1;
	real<lower=0, upper=50> shape_response_base;
	real<lower=0, upper=50> shape_response_delta_guess;
}
transformed parameters {
	real p_uspec;
	real p_uspec_unamb;
	p_uspec <- 0.000001;
	p_uspec_unamb <- 0.000001;
}
model
{
	int i;
	real prob_trial_type;
	real logProb_path;
	real shape_reading;
	real shape_response;
	real pY;
	
	for(i_obs in 1:n_obs)
	{
		int cur_subj;
		real logLik;
		logLik <- model_log_lik(iv_cond[i_obs], crit_region_cnt[i_obs], iv_questN1[i_obs], 
						   response_yes[i_obs], reading_time[i_obs], response_RT[i_obs],
						   p_pspec, p_uspec, p_uspec*p_uspec_unamb, 
						   p_arbitrary_attachment,
						   p_att_n1, p_retrieval_fail, p_guess_yes,
						   scale, shape_reading_base,
						   shape_delta_attach_N2, shape_delta_attach_N1, 
						   shape_response_base, shape_response_delta_guess);
		increment_log_prob(logLik);
	};

}