
generated quantities
{
  vector[n_obs] log_lik;
  vector[n_obs] pred_response;
  vector[n_obs] pred_reading_time;
  vector[n_obs] pred_response_RT;
  
  for(i_obs in 1:n_obs)
  {
    vector[3] pred;

    log_lik[i_obs] <- model_log_lik(iv_cond[i_obs], crit_region_cnt[i_obs], iv_questN1[i_obs], 
                                   response_yes[i_obs], reading_time[i_obs], response_RT[i_obs],
                                   p_pspec, p_uspec, p_uspec*p_uspec_unamb, 
                                   p_arbitrary_attachment,
                                   p_att_n1, p_retrieval_fail, p_guess_yes,
                                   scale, shape_reading_base,
                                   shape_delta_attach_N2, shape_delta_attach_N1, 
                                   shape_response_base, shape_response_delta_guess);

    pred <- model_pred_rng(iv_cond[i_obs], crit_region_cnt[i_obs], iv_questN1[i_obs], 
                                   response_yes[i_obs], reading_time[i_obs], response_RT[i_obs],
                                   p_pspec, p_uspec, p_uspec*p_uspec_unamb, 
                                   p_arbitrary_attachment,
                                   p_att_n1, p_retrieval_fail, p_guess_yes,
                                   scale, shape_reading_base,
                                   shape_delta_attach_N2, shape_delta_attach_N1, 
                                   shape_response_base, shape_response_delta_guess);
    pred_response[i_obs] <- pred[1];
    pred_reading_time[i_obs] <- pred[2];
    pred_response_RT[i_obs] <- pred[3];
  }
}