

results[["plots"]] <- list()

results[["plots"]][["uni_cox_bm"]] <-
  ggsurv_(
    srv_frm("bmab01"), 
    the_df = sw_mab(bm_cox), 
    pval = FALSE,
    pval.coord = c(42, 0.7)
  )

results[["plots"]][["uni_cox_ov"]] <- 
  ggsurv_(
    srv_frm("bmab01"), 
    the_df = sw_mab(ov_cox), 
    pval = FALSE
  )

results[["plots"]][["mv_cox_bm"]] <- 
  ggsurvadj_(
    srv_frm(prds$bm_cox), 
    the_df = sw_mab(bm_cox) 
    
  )

results[["plots"]][["mv_cox_ov"]] <- 
  ggsurvadj_(
    srv_frm(prds$o_cox), 
    the_df = sw_mab(ov_cox)
  )

#write_rds(results$plots, "results/plots.rds", compress = TRUE)


