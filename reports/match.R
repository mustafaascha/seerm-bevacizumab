
#readlines("Press enter to proceed with matching")

the_match <- 
  matchit(
    prop$frm, 
    data = (unmatched_data), 
    method = "nearest", 
    ratio = 1
  )

matched_data <- match.data(the_match) %>% srv_df() 
#m_ps <- xs_y(matched_data, frm = prop$frm)
matched_data[["propensity_score"]] <- 
  predict(glm(prop$frm, data = matched_data, family = "binomial"))
 # cv_glm_lasso(m_ps[["xs"]], m_ps[["y"]], family = "binomial") %>% 
 # (function(z) predict(z[["c_val"]]))

results[["overall_ps_plots"]] <- 
  map(list(unmatched = unmatched_data, matched = matched_data), 
      function(z) {
    ggplot(sw_mab(z), aes(x = propensity_score, fill = factor(bmab01))) + 
      geom_histogram() + 
      labs(x = "Propensity Score", y = "", fill = "Bev.") 
  })
attr(matched_data, "matched") <- "matched"
results[["ps_counts"]] <- 
  count_bmab_prds(prop$frm, unmatched_data, matched_data)

bm_prop_frm       <- update(prop$frm, . ~ . - d_ssg00 - RAD)
unmatched_bm_data <- 
  unmatched_data[with(unmatched_data, seer_br_mets == "SEER_Positive"),]
unmatched_bm_data[["propensity_score"]] <- 
  predict(glm(bm_prop_frm, data = unmatched_bm_data))

the_bm_match <- 
  matchit(bm_prop_frm, 
          data = (unmatched_bm_data), 
          method = "nearest", 
          ratio = 1)
matched_bm_data <- match.data(the_bm_match) %>% srv_df() 
attr(matched_bm_data, "matched") <- "matched"
matched_bm_data[["propensity_score"]] <- 
  predict(glm(bm_prop_frm, data = matched_bm_data, family = "binomial"))
results[["bm_ps_plots"]] <- 
  map(list(unmatched = unmatched_bm_data, matched = matched_bm_data), 
      function(z) {
    ggplot(sw_mab(z), aes(x = propensity_score, fill = factor(bmab01))) + 
      geom_histogram() + 
      labs(x = "Propensity Score", y = "", fill = "Bev.") 
  })

results[["bm_ps_counts"]] <- count_bmab_prds(bm_prop_frm, unmatched_bm_data, matched_bm_data)

write_rds(list(the_match,    matched_data),    "cache/propensity-data.rds")
write_rds(list(the_bm_match, matched_bm_data), "cache/propensity-bm-data.rds")


