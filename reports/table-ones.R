




cancers <- 
  reg_df(cancers) %>% 
  modify_at(
    c("race", "s_sex_v", "beh03v", "d_ssg00", 
     "csmetsdxb_pub_v", "csmetsdxliv_pub_v", "csmetsdxlung_pub_v"
     ),
    ~ ifelse(is.na(.x), "Missing", .x)
  )
mtch_dfs <- 
  c("unmatched_data", "matched_data", "unmatched_bm_data", "matched_bm_data")

results[["to"]] <- 
  map(mabsr::prds$strata, ~ tbl_one(cancers, mabsr::prds$to, .x)) %>% 
  set_names(names(mabsr::prds$strata)) #%>% 
  #map2(list(1:2, 1:2, 1, 1:4), ~ .x[,.y])

results[["tom"]] <- 
  map(mtch_dfs, ~ tbl_one(get(.x, envir = .GlobalEnv), strt = "bmab01") ) %>% 
  set_names(mtch_dfs)
results[["histo_key"]] <- histo_key(cancers)


