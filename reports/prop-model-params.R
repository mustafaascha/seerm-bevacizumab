


#readlines("Press enter to proceed with matched analysis")
message("Matching subjects and doing propensity-score analysis")
prop_mdl_params <- 
  list(
    data = list(
            add_matched(srv_df(unmatched_data), "unmatched"), 
            matched_data
           ), 
    prds = c("bmab01", "bmab01*seer_br_mets", "seer_br_mets", prds$o_cox), 
    wph  = c(TRUE, FALSE), 
    mdls = TRUE
  ) %>% 
cross()
         
bm_prop_mdl_params <- 
  list(
    data = list(
             add_matched(srv_df(unmatched_bm_data), "unmatched"), 
             matched_bm_data
           ),
    prds =  c("bmab01", prds$bm_cox), 
    wph  = c(TRUE, FALSE),
    mdls = TRUE
  ) %>% 
  cross()

 
