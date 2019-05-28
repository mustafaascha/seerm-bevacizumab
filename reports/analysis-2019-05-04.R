
devtools::load_all("mabsr")
load_libs()

switches(
  #do_imputation = TRUE,
  do_reload       = FALSE
  , do_glms       = FALSE 
  , do_coxs       = TRUE 
  , do_bm_cox     = TRUE 
  , do_srv_plt    = TRUE
  , do_match      = TRUE 
  , do_mtch_mdls  = TRUE
  , srs_analysis  = TRUE
  , do_tbl_ones   = TRUE 
  , write_results = TRUE
  , .force = T
  )

set.seed(1001)

if (do_reload | !exists("cancers")) source("reports/assemble-cancers.R") 
else cancers <- read_csv("cache/cancers-bev.csv.gz")

if (!exists("results")) results <- list()

# setup imputed data for models
if (!exists("unmatched_data")) { source("reports/unmatched-data.R") } 
if (do_glms)                   { source("reports/glms.R")           }

if (do_coxs) {
 # readline(prompt = "Press enter to proceed with Cox modeling")
  message("Doing Cox models")

  srv_fn <- function(df, ps) do_cox(df, ps, wph = FALSE, mdls = TRUE)

  net_mv <- function(.df, prds) {
    cn <- cv_cox_lasso(.df, prds)
    mv_prds <- select_coefs(cn)
    list(cnet = cn, prds = mv_prds, mv_mdl = srv_fn(.df, mv_prds))
  }
  
  message("coxnet_overall")
  results[c("coxnet_overall", "ov_prds", "mv_cox")] <- 
    net_mv(ov_cox, c(prds$cox, "surgery"))
  prds$o_cox <- results$ov_prds
  #results[["coxnet_overall"]] <- cv_cox_lasso(ov_cox, c(prds$cox, "surgery"))
  #prds$o_cox                  <- select_coefs(results$coxnet_overall)
  #message("mv_cox")
  #results[["mv_cox"]]         <- srv_fn(ov_cox, prds$o_cox)
  prds$strt_o_cox <- c(grep("d_ssg", prds$o_cox, value = T, invert = T), "strata(d_ssg00)")
  message("mv_cox_o_sss")
  results[["mv_cox_o_sss"]]   <- srv_fn(ov_cox, prds$strt_o_cox)

  message("coxnet_bm")
  results[c("coxnet_bm", "bm_prds", "bm_mv_cox")] <- 
    net_mv(bm_cox, c(prds$cox[-3], "surgery"))
  prds$bm_cox <- results$bm_prds
  #results[["coxnet_bm"]]      <- cv_cox_lasso(bm_cox, c(prds$cox[-3], "surgery"))
  #prds$bm_cox                 <- select_coefs(results$coxnet_bm)
  #message("bm_mv_cox")
  #results[["bm_mv_cox"]]      <- srv_fn(bm_cox, prds$bm_cox)

  message("coxnet_ov_mbs")
  results[c("coxnet_ov_mvs", "mbs_prds", "mv_cox_mbs")] <- 
    net_mv(mb_ov_cox, c(prds$cox, "surgery"))
#  results[["coxnet_ov_mbs"]] <- cv_cox_lasso(mb_ov_cox, prds$cox)
#  prds$mbs_cox               <- select_coefs(results$coxnet_overall)
#  message("mv_cox_mbs")
#  results[["mv_cox_mbs"]]    <- srv_fn(mb_ov_cox, prds$o_cox)
  prds$strt_o_mbs_cox <- c(grep("d_ssg", prds$o_cox, value = T, invert = T), "strata(d_ssg00)")
  message("mv_cox_o_sss_mbs")
  results[["mv_cox_o_sss_mbs"]] <- srv_fn(mb_ov_cox, prds$strt_o_mbs_cox)

  message("coxnet_bm_mbs")
  results[c("coxnet_bm_mbs", "bm_mbs_prds", "bm_mv_cox_mbs")] <- 
    net_mv(mb_bm_cox, c(prds$cox[-3], "surgery"))
  #results[["coxnet_bm_mbs"]]  <- cv_cox_lasso(mb_bm_cox, prds$cox[-3])
  #prds$bm_cox                 <- select_coefs(results$coxnet_bm)
  #message("bm_mv_cox_mbs")
  #results[["bm_mv_cox_mbs"]]  <- srv_fn(mb_bm_cox, prds$bm_cox)
}

source("reports/figures.R")

#propensity score-matched survival analysis=======================

if (do_match) {
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
} else {
  if (!exists("matched_data")) {
     c(the_match, matched_data) %<-% read_rds("cache/propensity-data.rds")
     c(the_bm_match, matched_bm_data) %<-% read_rds("cache/propensity-bm-data.rds")
  }
}

# NOTE:
cox.zph(coxph(Surv(srv_time, event = alive == "Dead") ~ bmab01, data = srv_df(unmatched_data)))
cox.zph(coxph(Surv(srv_time, event = alive == "Dead") ~ bmab01, data = srv_df(unmatched_bm_data)))


drug_prds <- 
  c(paste("n_", prds$drugs,sep = ""), paste("bmab01*n_", prds$drugs, sep = ""))

if (do_mtch_mdls) { 
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
           
  results[["propensity_hrs"]] <- map(prop_mdl_params, ~ do_prop_mdls(.x))
   
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

  results[["bm_propensity_hrs"]] <- map(bm_prop_mdl_params, ~ do_prop_mdls(.x))  
} else {
}

results[["bm_mv_cox"]]      <- srv_fn(bm_cox, prds$bm_cox)
results[["mv_cox_o_sss"]]   <- srv_fn(ov_cox, prds$strt_o_cox)

if (srs_analysis) source("reports/srs.R")

#table ones===============================================

if (do_tbl_ones) {
  cancers <- 
    reg_df(cancers) %>% 
    modify_at(
      c(
       "race", "s_sex_v", "beh03v", "d_ssg00", 
       "csmetsdxb_pub_v", "csmetsdxliv_pub_v", "csmetsdxlung_pub_v"
       ),
      ~ ifelse(is.na(.x), "Missing", .x)
    )
  mtch_dfs <- 
    c("unmatched_data", 
      "matched_data", 
      "unmatched_bm_data", 
      "matched_bm_data")

  results[["to"]] <- 
    map(mabsr::prds$strata, ~ tbl_one(cancers, mabsr::prds$to, .x)) %>% 
    set_names(names(mabsr::prds$strata)) #%>% 
    #map2(list(1:2, 1:2, 1, 1:4), ~ .x[,.y])

  results[["tom"]] <- 
    map(mtch_dfs, 
      function(df_nm) {
           tbl_one(get(df_nm, envir = .GlobalEnv), strt = "bmab01")
    }) %>% 
    set_names(mtch_dfs)
  results[["histo_key"]] <- histo_key(cancers)
} 

write_rds(results$to, "results/table-ones.rds")

results <- modify_if(results, is.data.frame, compose(as.data.frame, ungroup))


#if (write_results) {
  save(results, file = "results/models-mabs.rdata")
  write_rds(results, "results/models-mabs.rds", compress = "gz")
  writeLines("", "cache/analysis")
#}

