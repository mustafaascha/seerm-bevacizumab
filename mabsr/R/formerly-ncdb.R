

inline_num <- function(x) {
  the_x <- trimws(x,"both")
  if (all(grepl("[0-9]", unlist(strsplit(the_x, split = ""))))) {
    prettyNum(round(as.numeric(the_x),2), big.mark = ",")
  } else {
    x
  }
}

count_by <- function(.df, .by) {
  .df %>% 
    group_by_at(.vars = .by) %>% 
    summarise(the_count = n())
}

rename_count <- function(a_df, nm) {
  a_df[[nm]] <- a_df[["the_count"]]
  a_df[["the_count"]] <- NULL
  a_df
}

write_cache <- function(l, nm) {
  items <- grep(nm, names(l), value = T)
  l[["data"]] <- NULL
  write_rds(l[c("site", items)], 
            paste("cache/", nm, "-", Sys.Date(), ".rds", sep = "")
  )
  l
}

write_results <- function(obj) {
  nm <- deparse(substitute(obj))
  nm <- gsub("_", "-", nm)
  write_rds(obj, 
            paste("results/", nm, "-", Sys.Date(), ".rds", sep = "")
            )
}


autorecode <- function(a_df, nm) {
  stopifnot(all(a_df[[nm]] %in% c(0, 1, NA)))
  a_df[[paste(nm, "ar", sep = "_")]] <- 
    ifelse(a_df[[nm]] == 1, nm, paste("No", nm))
  a_df
}


prop_by <- function(.df, prop_over, .by) {
  po_nm <- paste(prop_over, "prop", sep = "_")
  
  numerator   <- 
    count_by(.df = .df, c(prop_over, .by)) %>% 
    rename_count("numerator")
  
  denominator <- 
    count_by(.df = .df, .by) %>%
    rename_count("denominator")
  
  #browser()
  left_join(numerator, denominator) %>% 
    mutate(!!po_nm := numerator / denominator)
  
}

rm_rnm_research <- function(z) {
  z <- z[,-1]
  rownames(z) <- gsub("\\/Research", "", rownames(z))
  z
}

manuscript_tbl_one <- function(a_df, frmla) {

  

  the_tbl <- 
    table_one(
      a_df, 
      stratum = vr, 
      table_vars = to_vars  
      #which_cols = 1:6
     )
  
  rownames(the_tbl) <- 
    reduce2(
      c("Cs Mets Dx", 
        "Sequence 00 = ", 
        "12 ", 
        "Cdcc Total Best", 
        " Cancer Program", 
        "Brain Rad", 
        "East", 
        "West", 
        " = Present"
        ), 
      c("Mets: ", 
        "", 
        "", 
        "Charlson Score", 
        "", 
        "Brain", 
        "E.", 
        "W.", 
        ""
        ), 
      function(a, b, d) gsub(b, d, a), 
      .init = rownames(the_tbl)
    )
  
  paste_cell <- function(a, r_nm) {
    if (length(a) > 1 & 
        as.numeric(a[1]) < 11 & 
        as.numeric(a[1]) > 0 & 
        !(any(grepl("Charlson|Rad", r_nm)))
        ) {
      "Censored"
    }
    else if (length(a) > 1) {
      paste(inline_num(as.numeric(a[1])), 
            paste(a[-1], collapse = ""))
    } else {
      paste(a, collapse = "")
    }
  }
  
  for (i in seq(ncol(the_tbl) - 1)) {
    the_tbl[,i] <- 
      map2_chr(strsplit(trimws(the_tbl[,i]), "\\ "), 
               rownames(the_tbl), 
               paste_cell)
  }
  the_tbl
  
}


post_recode <- function(ncdb_df, to_one_hot = NULL) {
  if (missing(to_one_hot)) to_one_hot <- ncdbp::to_one_hot
  c(prvs_nms, ncdb_df) %<-% one_hot_ncdb(ncdb_df, to_one_hot)
  
  ncdb_df[["year_of_diagnosis"]] <- 
    factor(as.character(ncdb_df[["year_of_diagnosis"]]), 
           levels = as.character(2010:2015)
           #, ordered = TRUE
    )
  
  ncdb_df[["facility_type"]] <- 
    reduce2(c(AR = "Academic/Research Program",
              CC = "Comprehensive Community Cancer Program",
              IN = "Integrated Network Cancer Program", 
              C = "Community Cancer Program"), 
            c("^AR$", "^CC$", "^IN$", "^C$"),
            function(x, y, z) gsub(z, y, x),
            .init = as.character(ncdb_df[["facility_type_cd_v"]]))
  
  ncdb_df[["age"]] <- 
    cut(as.numeric(ncdb_df[["age"]]), 
        breaks = c(18, 30 + 0:4 * 10, 90), 
        labels = c("18-30", 
                   "30-40", "40-50", 
                   "50-60", "60-70", 
                   "70-90"))
  
  ncdb_df <- 
    reduce(
      c("Stereotactic_RS", "Gamma_Knife"),
      autorecode, 
      .init = ncdb_df
    )
  
  ncdb_df
  
}

count_bm_by <- function(a_df, gp) {
  gp <- enquo(gp)
  a_df %>% 
    group_by(site, brain_mets) %>% 
    nest() %>% 
    pmap(function(site, brain_mets, data) {
      data %>%
        group_by(year_of_diagnosis, !!gp) %>%
        summarise(count = n()) %>%
        mutate(site = site,
               brain_mets = brain_mets)
    })
  
}


xs_y <- function(model_df, outcome, prds, pf_prd = "", ..., frm = NULL) {
  if (!missing(frm)) {
    outcome <- all.vars(update(frm, . ~ 1))
    prds <- all.vars(update(frm, 1 ~ .))
  }
  mdf_fct  <- map_df(drop_na(model_df), as.factor)
  these <- function(vs) data.matrix(select(mdf_fct, one_of(vs)))
  r <- list(xs = these(prds), y  = these(outcome))
  r[["pf"]] <- !(grepl(paste(pf_prd, collapse = "|"), colnames(r$xs)))
  r
}

mk_frm <- function(prds, ot) {
  as.formula(
    paste(ot, "~", paste(prds, collapse = "+"))
  )
}

mv_mdl <- function(fn) {
  function(dta, model_predictors) {
    dta <- relevel_predictors(dta)
    fn(mk_frm(model_predictors, "Stereotactic_RS"), data = dta) %>%
      tidy_glm() %>%
      mutate(model = "Multivariable")
  }
}

mv_glm <- mv_mdl(glm)
mv_speedglm <- mv_mdl(speedglm)
mv_biglm    <- mv_mdl(biglm)

uv_mdl <- function(fn) {
  function(dta, model_predictors, rlvl = TRUE) {
    if(rlvl) dta <- relevel_predictors(dta)
    map_dfr(model_predictors,
        function(prd) {
          map_dfr("Stereotactic_RS", #ncdbp::treatments, 
              function(otcm) {
                fn(mk_frm(prd, otcm), data = dta) %>%
                  tidy_glm() %>%
                  mutate(model = "Univariable", 
                         outcome = otcm)
              })
        })
  }
}

uv_glm      <- uv_mdl(glm)
uv_speedglm <- uv_mdl(speedglm)
uv_biglm    <- uv_mdl(biglm)

mk_others <- function(a_df) {
  modify_if(a_df, 
            is.character, 
            ~ forcats::fct_lump(.x, prop = 0.05, other_level = "Other"))
}

impute_l <- function(the_dfs, the_nms, model_predictors) {
  dfs <- 
    map(the_dfs, 
        function(z) {
          if (nrow(z) > 5e5) z <- sample_n(z, 5e5)
          select(z, one_of(model_predictors, "Stereotactic_RS"))
    })
  
  dfs <- 
    mclapply(dfs, 
             function(.x) complete(mice(.x, 2, seed = 10), action = "long")
             )
    
  dfs <- 
    dfs %>% 
    set_names(the_nms) %>% 
    map(mk_others) %>% 
    map(~ filter(.x, .imp == 1) %>% 
          select(-.id))
}

impute_l_mr <- function(the_dfs, the_nms, model_predictors) {
  dfs <- 
    map(the_dfs, 
        function(z) {
          select(z, one_of(model_predictors, "Stereotactic_RS"))
    })

  dfs %>%
    map(~ missRanger(.x, 
                     10, 
                     seed = 10, 
                     verbose = 1, 
                     pmm.k = 3,
                     sample.fraction = 0.05, 
                     num.trees = 100)
        ) %>%
    set_names(the_nms) %>%
    map(mk_others) %>% 
    map(~ mutate_all(.x, as.character))
}


mdl_prd_arg <- function(...) partial(..., model_predictors = model_predictors)

mk_glmnet_params <- function(glmnet_datum) {
  list(x = glmnet_datum[,-1],
       y = glmnet_datum[,1],
       family = "binomial",
       nfolds = 10#,
       #weights = rep_len(0.1, nrow(glmnet_datum[["xs"]]))
       #, parallel = TRUE
  )
}

impute_2 <- function(df_l, nms, prds, mc = FALSE) {

  if (mc) registerDoMC(12)
  
  imputed_data <- 
    impute_l_mr(df_l, 
                nms,
                model_predictors = prds) 
  
  variance <- 
    map_dfc(imputed_data, ~ map_int(.x, compose(length, unique, na.omit))) %>% 
    mutate(col_position = 1:n())

  map2(variance[,-grep("^col_position$", names(variance))], 
       imputed_data, 
       function(x, y) {
         acceptable_variance <- variance[["col_position"]][x > 1]
         y[,acceptable_variance]
       }) %>% 
  map(imputed_data, xs_y)
  
}

check_imputation <- function(i_matrix) {
  apply(i_matrix, 2, compose(length, which, is.na))
}

check_imputation_l <- function(i_matrices, .all = T) {
  if (.all) rtrn <- function(z) all(flatten_int(z) == 0)
  else rtrn <- function(z) z
  rtrn(map(i_matrices, check_imputation))
}


rm_low_variance <- function(dfs) {
  message("Removing low-variance columns")
  map_dfc(dfs, ~ map_int(.x, compose(length, unique, na.omit))) %>% 
    map(~ which(.x < 2)) %>% 
    map2(dfs, function(lv, i_df) select(i_df, -lv))
}


#srv_df <- function(nc_df) {
#  nc_df[["srv_months"]] <- nc_df[["dx_lastcontact_death_months"]]
#  nc_df[["censor"]] <- nc_df[["puf_vital_status_v"]] != "Alive"
#  nc_df[["srv"]] <- 
#    Surv(nc_df[["srv_months"]], 
#         event = nc_df[["censor"]], 
#         type = "right")
#  nc_df
  #filter(nc_df, year_of_diagnosis != "2015")
#}

ncdb_int_srv <- function(prd, dta) {
  dta  <- srv_df(relevel_predictors(dta))
  message(paste("Making a model for", prd))
  message(paste(sum(dta[[prd]], na.rm = T), prd, "treatments were found"))
  print(table(dta[,c(prd, "year_of_diagnosis")]))
  frm <- reformulate(paste(prd, "*year_of_diagnosis", sep = ""), "srv")
  coxph(frm, data = dta) 
}

tidy_int_srv <- function(mdl) {
  mdl %>% 
    broom::tidy() %>% 
    filter(grepl(":", term))
}

ncdb_srv <- function(prds, dta) {
  dta <- 
    relevel_predictors(dta) #%>% filter(!(age %in% c("18-30", "30-40")))
  dta <- dta[!(dta[["age"]] %in% c("18-30", "30-40")),]
  dta[["age"]] <- relevel(factor(as.character(dta[["age"]])), ref = "40-50")
  #  filter(!is.na(age) & !is.na(location) & !is.na(facility_type) & !is.na(urban) & !is.na(income_12) & 
  #           race != "Unknown")
  # rm_these <- 
  #   dta %>% 
  #   group_by_at(prds) %>% 
  #   summarise(the_count = n()) %>% 
  #   filter(the_count < 10) %>% 
  #   select(-the_count) %>% 
  #   gather(k, v) %>% 
  #   distinct()
  # browser()
  message(paste("Making a model for", paste(prds, collapse = ", ")))
  s <- function(...) safely(...)
  uni_cph <- function(prd) s(coxph)(reformulate(prd, "srv"), data = dta)[["result"]]
  mul_cph  <- function()    s(coxph)(reformulate(prds, "srv"), data = dta)[["result"]]
  list(univariable = map_dfr(prds, compose(broom::tidy, uni_cph)),
       multivariable = broom::tidy(mul_cph())) %>% 
  map2_dfr(c("univariable", "multivariable"), 
         ~ mutate(.x, model = .y)) %>% 
    modify_at(c("estimate", "conf.low", "conf.high"), exp)
}

relevel_predictors <- function(a_df) {
#  a_df[["race"]] <- relevel(factor(a_df[["race"]]), ref = "White")
#  a_df[["insurance"]] <- 
#    relevel(factor(a_df[["insurance"]]), ref = "Private Insurance")
#  a_df[["facility_type"]] <- 
#    relevel(factor(a_df[["facility_type"]]), ref = "Comprehensive Community Cancer Program")
#  a_df[["urban"]] <- relevel(factor(a_df[["urban"]]), ref = "Metro")
#  a_df[["year_of_diagnosis"]] <- relevel(factor(a_df[["year_of_diagnosis"]]), ref = "2010")
  a_df
}



pander_srv <- function(mdl_df) {
  mdl_df %>% 
    select(-estimate, -conf.low, -conf.high, -statistic, -std.error) %>% 
    select(site, brain_mets, model, term, hr, p.value) %>% 
    pander()
}

tidy_tidied_srv <- function(a_df) {
  dig <- dig_fns()
  a_df %>% 
    modify_at(c("estimate", "conf.low", "conf.high"), compose(dig$two, exp)) %>% 
    modify_at("p.value", dig$four) %>% 
    mutate(hr = paste(estimate, " (", conf.low, "-", conf.high, ")", sep = "")) 
}

clean_term <- function(trm) {
  reduce2(
    c("age", 
      "facility_type",
      "location",
      "race",
      "year_of_diagnosis",
      "urban", 
      "income_12", 
      "insurance", 
      "West", 
      "East", 
      " Cancer Program", 
      "Private Insurance", 
      "Other Government"
      ), 
    c("Age: ", 
      "Facility: ",
      "Region: ",
      "Race: ",
      "Year: ",
      "Location: ",
      "Income: ",
      "Insurance: ", 
      "W.", 
      "E.", 
      "", 
      "Private", 
      "Gov. NOS"
      ), 
    function(a, b, d) gsub(b, d, a), 
    .init = trm
  )
}

show_srs_or <- function(a_df) {
  d <- dig_fns()
  a_df %>% 
    modify_at(c("estimate", "conf.low", "conf.high", "p.value"), as.numeric) %>% 
    modify_at(c("estimate", "conf.low", "conf.high"), d$two) %>% 
    modify_at("p.value", d$three) %>% 
    modify_at("term", clean_term) %>% 
    select(term, estimate, conf.low, conf.high, p.value) %>% 
    pander(split.table = Inf)
}

show_srs_hr <- function(a_df, mdl = "univariable") {
  d <- dig_fns()
  a_df %>% 
    filter(model %in% mdl) %>% 
    select(-model) %>% 
    modify_at(c("estimate", "conf.low", "conf.high", "p.value"), as.numeric) %>% 
    modify_at(c("estimate", "conf.low", "conf.high"), d$two) %>% 
    modify_at("p.value", d$three) %>% 
    modify_at("term", clean_term) %>% 
    filter(!grepl("Term|2015", term)) %>% 
    select(term, estimate, conf.low, conf.high, p.value) %>% 
    pander(split.table = Inf)
}

timewrite_rds <- function(z, fp) {
  fp <- paste( fp, gsub("[^a-zA-Z0-9]", "-", Sys.time()), ".rds", sep = "")
  write_rds(z, path = fp)
}

mk_choro_df <- function(era_mdls) {
  list(
    `West North Central` = c("IA", "KS", "MN", "MO", "ND", "NE", "SD"),
    `South Atlantic` = c("DC", "DE", "FL", "GA", "MD", "NC", "SC", "VA", "WV"),
    `Pacific` = c("AK", "CA", "HI", "OR", "WA"),
    `New England` = c("CT", "MA", "ME", "NH", "RI", "VT"),
    `Mountain` = c("AZ", "CO", "ID", "MT", "NM", "NV", "UT", "WY"),
    `East North Central` = c("IL", "IN", "MI", "OH", "WI"),
    `East South Central` = c("AL", "KY", "MS", "TN"),
    `Middle Atlantic` = c("NJ", "NY", "PA"),
    `West South Central` = c("AR", "LA", "OK", "TX")
  ) %>% 
    imap_dfr(~ data.frame(level = .y, state.abb = .x, stringsAsFactors = FALSE)) %>% 
    left_join(data.frame(state = state.name, state.abb, stringsAsFactors = FALSE)) %>% 
    full_join(era_mdls[["models"]][[which(era_mdls$the_grouping== "location")]]) %>% 
    select(-state.abb) %>% 
    modify_at("state", ~ factor(tolower(.x))) %>% 
    filter(!(state %in% c("alaska", "hawaii"))) 
}

states_plot <- function(states_df) {
  ggplot(states_df, aes(map_id = state)) + 
    # map points to the fifty_states shape data
    geom_map(aes(fill = estimate), map = fifty_states) + 
    expand_limits(x = fifty_states$long, y = fifty_states$lat) +
    coord_map() +
    facet_wrap("year") + 
    scale_x_continuous(breaks = NULL) + 
    scale_y_continuous(breaks = NULL) +
    #scale_fill_viridis_c(end = 0.9) + 
    scale_fill_continuous(breaks = 13:7 / 10, type = "viridis", end = 0.9) + 
    labs(x = "", y = "", fill = "SRS Odds Ratio") +
    theme(legend.position = c(0.85, 0.3), 
          panel.background = element_blank(), 
          strip.text = element_text(color = "black"),
          legend.title = element_text(size = 10)
    ) + 
    guides(fill = guide_colorbar(barheight = 3))
}

era_plots <- function(era_df, continuous = FALSE) {
  era_df <- era_df %>% filter(!is.na(level))
  if (continuous) {
    p <- 
      ggplot(era_df, aes(x = level, y = estimate)) + 
        geom_point(position = position_dodge(1)) + 
        geom_errorbar(aes(ymax = conf.high, ymin = conf.low)) + 
        geom_hline(yintercept = 1, color = "grey20", linetype = 2) + 
        labs(x = "", y = "Odds Ratio (95% CI)")
  } else {
    p <- 
      ggplot(era_df, aes(x = year, y = estimate, color = level)) + 
        geom_point(position = position_dodge(1)) + 
        geom_errorbar(aes(ymax = conf.high, ymin = conf.low)) + 
        geom_hline(yintercept = 1, color = "grey20", linetype = 2) + 
        facet_wrap(~ level) + 
        # theme(legend.position = "bottom", 
        #       strip.text.x = element_blank(), 
        #       axis.text.x = element_text(angle = 35, hjust = 1)
        #       ) + 
        labs(x = "", y = "Odds Ratio (95% CI)")
  }
  p
}

munge_era_models <- function(the_mdls, continuous = FALSE) {
  if (continuous) {
    c_fn <- function(z) rename(z, the_grouping = k, level = v)
  } else {
    c_fn <- function(z) {
      rename(z, the_grouping = k, level = v, year = term) %>% 
      modify_at("year", ~ as.numeric(gsub("year_of_diagnosis", "", .x)))
    }
  }
  the_mdls %>% 
    select(-data) %>% 
    unnest() %>% 
    filter(term != "(Intercept)") %>% 
    select(-statistic, -std.error) %>% 
    c_fn() %>% 
    modify_at(c("estimate", "conf.low", "conf.high"), compose(exp, as.numeric)) %>% 
    group_by(the_grouping) %>% 
    nest(.key = "models")
}

mk_era_mdl <- function(the_df){
  mdl <- speedglm(Stereotactic_RS ~ year_of_diagnosis, data = the_df)
  bind_cols(broom::tidy(mdl), broom::confint_tidy(mdl)) %>% 
    mutate_all(as.character)
}

mk_time_dfs <- function(predictors, srs_mdls, continuous = FALSE) {
  if (continuous) yr_fn <- function(a) as.numeric(a)
  else yr_fn <- function(a) a
  dta <- 
    srs_mdls[["data"]][[1]] %>% 
    modify_at("year_of_diagnosis", yr_fn) %>% 
    modify_at("facility_type", ~ gsub(" Cancer Program", "", .x))
  map_dfr(predictors, 
          function(prd) {
            group_by_at(dta, prd) %>% 
              nest() %>% 
              gather(k, v, -data)
          }) %>% 
    filter(map_int(data, nrow) > 10) %>% 
    filter(map_dbl(data, ~ sum(.x[["Stereotactic_RS"]], na.rm = T) / 5) > 10)
}

table_vs_year <- function(the_dfs) {
  f <- function(the_df) {
    table(the_df[["puf_vital_status_v"]], the_df[["year_of_diagnosis"]]) %>% 
      as.data.frame()
  }
  
  the_dfs %>%
    mutate(vital_status = map(data, f)) %>% 
    select(brain_mets, vital_status) %>% 
    unnest()
}

table_srs_yr <- function(the_df) {
  table(the_df[["year_of_diagnosis"]], the_df[["Stereotactic_RS_ar"]]) %>% 
  as.data.frame()
}

apply_criteria <- function(the_df) {
  filter(
    the_df, 
      treatment_status == "Treatment_Given" & 
        site %in% "LungNSC" & 
        brain_mets == "Yes"
    )
}

mk_models_over_time <- function(nested_df, continuous = FALSE) {
  mk_time_dfs(predictors, nested_df, continuous = continuous) %>% 
    mutate(mdls = map(data, mk_era_mdl)) %>% 
    munge_era_models(continuous = continuous)
}

double_tidy <- function(a_df, nm) {
  dig <- dig_fns()
  a_df %>% 
    modify_at(c("estimate", "conf.low", "conf.high"), compose(dig$two, as.numeric)) %>% 
    mutate(!!nm := paste(estimate, " (", conf.low, "-", conf.high, " p:", p.value, ")", sep = "")) %>%
    select(term, !!nm)
}

show_all_srs_or <- function(mv, uv) {
  map2(list(mv, uv), c("Multivariable", "Univariable"), double_tidy) %>% 
    reduce(full_join) %>% 
    modify_at("term", clean_term)
}

show_all_srs_hr <- function(srvs) {
  map(c("univariable", "multivariable"), 
       function(mdl) {
         d <- dig_fns()
         srvs %>% 
           filter(model %in% mdl) %>% 
           select(-model) %>% 
           modify_at(c("estimate", "conf.low", "conf.high", "p.value"), as.numeric) %>% 
           modify_at(c("estimate", "conf.low", "conf.high"), d$two) %>% 
           modify_at("p.value", d$three) %>% 
           modify_at("term", clean_term) %>% 
           filter(!grepl("Term|2015", term)) %>% 
           select(term, estimate, conf.low, conf.high, p.value) %>% 
           double_tidy(mdl)
       }) %>% 
    reduce(left_join)
}


