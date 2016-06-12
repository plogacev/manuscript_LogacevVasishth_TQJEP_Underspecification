
generated quantities
{
  vector[n_obs] log_lik;
  vector[n_obs] pred_response;
  vector[n_obs] pred_reading_time;
  vector[n_obs] pred_response_RT;
  
  for(i_obs in 1:n_obs)
  {
    int cur_subj;
    vector[3] pred;
    cur_subj <- subj[i_obs];
    
    log_lik[i_obs] <- model_log_lik(iv_cond[i_obs], crit_region_cnt[i_obs], iv_questN1[i_obs], 
                                    response_yes[i_obs], reading_time[i_obs], response_RT[i_obs],
                                    p_pspec, p_uspec_subj[cur_subj], p_uspec_subj[cur_subj]*p_uspec_unamb_subj[cur_subj], 
                                    p_arbitrary_attachment_subj[cur_subj],
                                    p_att_n1_subj[cur_subj], p_retrieval_fail_subj[cur_subj], p_guess_yes_subj[cur_subj],
                                    scale_subj[cur_subj], shape_reading_base_subj[cur_subj], 
                                    shape_delta_attach_N2_subj[cur_subj], shape_delta_attach_N1_subj[cur_subj], 
                                    shape_response_base_subj[cur_subj], shape_response_delta_guess_subj[cur_subj]);
    
    pred <- model_pred_rng(iv_cond[i_obs], crit_region_cnt[i_obs], iv_questN1[i_obs], 
                           response_yes[i_obs], reading_time[i_obs], response_RT[i_obs],
                           p_pspec, p_uspec_subj[cur_subj], p_uspec_subj[cur_subj]*p_uspec_unamb_subj[cur_subj], 
                           p_arbitrary_attachment_subj[cur_subj],
                           p_att_n1_subj[cur_subj], p_retrieval_fail_subj[cur_subj], p_guess_yes_subj[cur_subj],
                           scale_subj[cur_subj], shape_reading_base_subj[cur_subj], 
                           shape_delta_attach_N2_subj[cur_subj], shape_delta_attach_N1_subj[cur_subj], 
                           shape_response_base_subj[cur_subj], shape_response_delta_guess_subj[cur_subj]);
    
    pred_response[i_obs] <- pred[1];
    pred_reading_time[i_obs] <- pred[2];
    pred_response_RT[i_obs] <- pred[3];
  }
}