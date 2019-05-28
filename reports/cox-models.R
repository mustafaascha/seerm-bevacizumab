
# readline(prompt = "Press enter to proceed with Cox modeling")
message("Doing Cox models")

srv_fn <- function(df, ps) do_cox(df, ps, wph = FALSE, mdls = TRUE)

net_mv <- function(.df, prds, also = NULL) {
  cn <- cv_cox_lasso(.df, prds)
  mv_prds <- select_coefs(cn)
  uni <- map_dfr(c(prds, also), ~ do_cox(.df, .x, wph = FALSE, mdls = FALSE))
  mv <- srv_fn(.df, mv_prds)
  mv_prds <- mv_prds[!grepl("Stereo_RS|bmab01", mv_prds)]
  mv_int <- srv_fn(.df, c(mv_prds, "Stereo_RS*bmab01"))
  list(cn, mv_prds, mv_mdl = mv, uni, mv_int = mv_int)
}

message("coxnet_overall")
results[c("coxnet_overall", "ov_prds", "mv_cox", "uni_cox", "mv_cox_int")] <- 
  net_mv(ov_cox, c(prds$cox, "surgery"), also = "bmab01*Stereo_RS")
prds$o_cox <- results$ov_prds
prds$strt_o_cox <- c(grep("d_ssg", prds$o_cox, value = T, invert = T), "strata(d_ssg00)")
message("mv_cox_o_sss")
results[["mv_cox_o_sss"]]   <- srv_fn(ov_cox, prds$strt_o_cox)

message("coxnet_bm")
results[c("coxnet_bm", "bm_prds", "bm_mv_cox", "bm_uni_cox", "bm_mv_cox_int")] <- 
  net_mv(bm_cox, c(prds$cox[-3], "surgery"), also = "bmab01*Stereo_RS")
prds$bm_cox <- results$bm_prds

message("coxnet_ov_mbs")
results[c("coxnet_ov_mvs", "mbs_prds", "mv_cox_mbs", "uni_cox_mbs", "mv_cox_mbs_int")] <- 
  net_mv(mb_ov_cox, c(prds$cox, "surgery"), also = "bmab01*Stereo_RS")
prds$strt_o_mbs_cox <- c(grep("d_ssg", prds$o_cox, value = T, invert = T), "strata(d_ssg00)")
message("mv_cox_o_sss_mbs")
results[["mv_cox_o_sss_mbs"]] <- srv_fn(mb_ov_cox, prds$strt_o_mbs_cox)

message("coxnet_bm_mbs")
results[c("coxnet_bm_mbs", "bm_mbs_prds", "bm_mv_cox_mbs", "bm_uni_cox_mbs", "bm_mv_cox_int")] <- 
  net_mv(mb_bm_cox, c(prds$cox[-3], "surgery"), also = "bmab01*Stereo_RS")



cox.zph(coxph(Surv(srv_time, event = alive == "Dead") ~ bmab01, data = srv_df(unmatched_data)))
cox.zph(coxph(Surv(srv_time, event = alive == "Dead") ~ bmab01, data = srv_df(unmatched_bm_data)))


#results[["coxnet_overall"]] <- cv_cox_lasso(ov_cox, c(prds$cox, "surgery"))
#prds$o_cox                  <- select_coefs(results$coxnet_overall)
#message("mv_cox")
#results[["mv_cox"]]         <- srv_fn(ov_cox, prds$o_cox)
#results[["coxnet_bm"]]      <- cv_cox_lasso(bm_cox, c(prds$cox[-3], "surgery"))
#prds$bm_cox                 <- select_coefs(results$coxnet_bm)
#message("bm_mv_cox")
#results[["bm_mv_cox"]]      <- srv_fn(bm_cox, prds$bm_cox)
#  results[["coxnet_ov_mbs"]] <- cv_cox_lasso(mb_ov_cox, prds$cox)
#  prds$mbs_cox               <- select_coefs(results$coxnet_overall)
#  message("mv_cox_mbs")
#  results[["mv_cox_mbs"]]    <- srv_fn(mb_ov_cox, prds$o_cox)
#results[["coxnet_bm_mbs"]]  <- cv_cox_lasso(mb_bm_cox, prds$cox[-3])
#prds$bm_cox                 <- select_coefs(results$coxnet_bm)
#message("bm_mv_cox_mbs")
#results[["bm_mv_cox_mbs"]]  <- srv_fn(mb_bm_cox, prds$bm_cox)

