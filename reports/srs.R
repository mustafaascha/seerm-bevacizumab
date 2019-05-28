
dts <- 
  cancers %>% 
  mutate(dxd = lubridate::mdy(paste(dx_month, "01", dx_year, sep = "-"))) %>%
  select(patient_id, mab_start, str_start, dxd)

# CHANGE DATE SO THAT THE START TIME IS DATE OF FIRST SRS RX

an_dt <- function(a, b) as.numeric(difftime(a, b, units = "days"))
zero_default <- function(z) ifelse(z <= 0, 0.1, z)

surv_dfs <- 
  list(
    bm = bm_cox, 
    ov = ov_cox, 
    bm_matched = matched_bm_data, 
    ov_matched = matched_data
  ) %>% 
  map(
    function(.df) {
      left_join(.df, dts, by = "patient_id") %>% 
      filter(Stereo_RS == "SRS" & bmab01 == 1) %>% 
      mutate(
        after        = mab_start > str_start
        , tt_srs_bev   = an_dt(mab_start, str_start)
        , mab_near_srs = abs(tt_srs_bev) < 30
        , tt_dx_srs    = an_dt(dxd, str_start)
        , tt_dx_srs_m  = an_dt(dxd, str_start) / 30
        , tt_dx_mab    = an_dt(dxd, mab_start)
      ) %>% 
      left_join(
        select(cancers, 
               patient_id, n_str_crs, n_mab_crs, mab_durat, str_durat, n_str_crs
        ), 
        by = "patient_id"
      )
    }
  ) 

surv_dfs[paste(names(surv_dfs), "srsrv", sep = "_")] <- 
  lapply(surv_dfs, function(z) mutate(z, srv_time = zero_default(srv_time - tt_dx_srs_m)))

results[["srs_mab_timing"]] <- 
  surv_dfs %>% 
  map(
    ~ list(
        after          = do_cox(.x, "after"       , wph = F, mdls = T),
        tt_srs_bev     = do_cox(.x, "tt_srs_bev"  , wph = F, mdls = T), 
        bev_near_srs   = do_cox(.x, "mab_near_srs", wph = F, mdls = T),
        after_w        = do_cox(.x, "after"       , wph = T, mdls = T),
        tt_srs_bev_w   = do_cox(.x, "tt_srs_bev"  , wph = T, mdls = T), 
        bev_near_srs_w = do_cox(.x, "mab_near_srs", wph = T, mdls = T)
      )
  )

results[["srs_strict_mab_timing"]] <- 
  surv_dfs %>% 
  map(~ filter(.x, n_str_crs > 1 & n_mab_crs > 1)) %>% 
  map(
    ~ list(after        = do_cox(.x, "after"       , wph = F, mdls = T),
           tt_srs_bev   = do_cox(.x, "tt_srs_bev"  , wph = F, mdls = T), 
           bev_near_srs = do_cox(.x, "mab_near_srs", wph = F, mdls = T),
           after_w        = do_cox(.x, "after"       , wph = T, mdls = T),
           tt_srs_bev_w   = do_cox(.x, "tt_srs_bev"  , wph = T, mdls = T), 
           bev_near_srs_w = do_cox(.x, "mab_near_srs", wph = T, mdls = T)
      )
  )

results[["srs_mab_timeto_dfs"]] <- 
  surv_dfs %>% 
  map(~ select(.x, tt_srs_bev, Stereo_RS, bmab01, srv_time_mon, alive))

library("magrittr")
cancers %<>% 
  mutate(dxd = lubridate::mdy(paste(dx_month, "01", dx_year, sep = "-"))) %>%
  mutate(
    after        = mab_start > str_start
    , tt_srs_bev   = an_dt(mab_start, str_start)
    , mab_near_srs = abs(tt_srs_bev) < 30
    , tt_dx_srs    = an_dt(dxd, str_start)
    , tt_dx_mab    = an_dt(dxd, mab_start)
  )

results[["srs_mab_df"]] <- 
  cancers %>% 
  mutate(dxd = lubridate::mdy(paste(dx_month, "01", dx_year, sep = "-"))) %>%
  mutate(
    after        = mab_start > str_start
    , tt_srs_bev   = an_dt(mab_start, str_start)
    , mab_near_srs = abs(tt_srs_bev) < 30
    , tt_dx_srs    = an_dt(dxd, str_start)
    , tt_dx_mab    = an_dt(dxd, mab_start)
  ) %>% 
  select(tt_srs_bev, after, mab, n_mab_crs, mab_durat, str_durat, n_str_crs, 
        srv_time_mon, stat_rec_v, seer_br_mets, tt_dx_mab)

results[["srs_mab_df_mbs"]] <- 
  mb_ov_cox %>% 
  left_join(
    cancers %>% 
    mutate(dxd = lubridate::mdy(paste(dx_month, "01", dx_year, sep = "-"))) %>%
    select(patient_id, mab_start, mab_durat, str_durat, str_start, n_mab_crs, n_str_crs, dxd)
  ) %>% 
  mutate(dxd = lubridate::mdy(paste(dx_month, "01", dx_year, sep = "-"))) %>%
  mutate(
    after        = mab_start > str_start
    , tt_srs_bev   = an_dt(mab_start, str_start)
    , mab_near_srs = abs(tt_srs_bev) < 30
    , tt_dx_srs    = an_dt(dxd, str_start)
    , tt_dx_mab    = an_dt(dxd, mab_start)
  ) %>% 
  select(tt_srs_bev, after, bmab01, n_mab_crs, mab_durat, str_durat, n_str_crs, 
        srv_time_mon, stat_rec_v, seer_br_mets, tt_dx_mab)

#c_mb <- nest(cancers, -mab)

map_mabs <- function(.fn, vr, .gp) {
  .gp <- enquo(.gp)
  c_mb <- nest(cancers, -!!.gp)
  map(c_mb[["data"]], function(z) .fn(z[[vr]])) %>% set_names(c_mb[[quo_name(.gp)]])
}

results[
  c("n_mab_crs", "n_str_crs", "dur_mab_crs", "dur_str_crs"
  , "tt_dx_mab", "tt_dx_srs", "tt_srs_bev")
] <- 
  lapply(
  c("n_mab_crs",  "n_str_crs", "mab_durat", "str_durat"
    ,"tt_dx_mab", "tt_dx_srs", "tt_srs_bev"),
  function(z) summary(cancers[[z]])
  )

results[["dur_str_crs_strat"]] <- map_mabs(summary, "str_durat", mab)
results[["n_str_crs_strat"]]   <- map_mabs(summary, "n_str_crs", mab)
results[["tt_dx_srs_strat"]]   <- map_mabs(summary, "tt_dx_srs", mab)
results[["n_mab_crs_strat"]]   <- map_mabs(summary, "n_mab_crs", seer_br_mets)
results[["dur_mab_crs_strat"]] <- map_mabs(summary, "mab_durat", seer_br_mets)
results[["tt_dx_mab_strat"]]   <- map_mabs(summary, "tt_dx_mab", seer_br_mets)
results[["tt_srs_bev_strat"]]  <- map_mabs(summary, "tt_dx_srs", seer_br_mets)

