
glm_df      <- cancers %>% reg_df() 
imputation  <- impute(glm_df, prds$candidates)

cox_df      <- 
  cancers %>% 
  reg_df() %>% 
  impute(a_df = ., predictors = c(prds$impute, "brain_mets"))

unmatched_data <- 
  cox_df[["imputed"]] %>% 
  modify_at(c("srv_time_mon", "srv_time_mab"), as.numeric)
unmatched_data$id <- 1:nrow(unmatched_data)
unmatched_data[["propensity_score"]] <- predict(glm(prop$frm, data = unmatched_data))
attr(unmatched_data, "matched") <- "unmatched"

keepers <- c("bmab01", "dexamethasone")

bm_glm      <- imputation$imputed %>% filter(brain_mets == "Yes")
gnet        <- xs_y(imputation$imputed, "bmab01", rm_pr("brain|bone|liver"), keepers)
gnet_bm     <- xs_y(bm_glm, "bmab01", rm_pr("bone|liver|brain|d_ssg"), keepers)

ov_cox      <- cox_df[["imputed"]] %>% srv_df()
bm_cox      <- ov_cox[with(ov_cox, seer_br_mets == "SEER_Positive"),] %>% srv_df()

mb_ov_cox <- srv2mbs(ov_cox, cancers) %>% filter(srv_time > 0 | is.na(srv_time))
mb_bm_cox <- srv2mbs(bm_cox, cancers) %>% filter(srv_time > 0 | is.na(srv_time))


