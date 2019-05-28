

summarise_by_year <- function(a_df, facet_gp, otcm_gp) {
  
  # count_nm <- paste(quo_name(enquo(a_df)), "count", sep = "_")
  # prop_nm <- paste(quo_name(enquo(a_df)), "prop", sep = "_")
  count_nm <- "treated_count"
  prop_nm <- "treated_prop"
  
  numerator_df <-  
    group_by_(a_df, facet_gp, otcm_gp, "year_of_diagnosis") %>% 
    summarise(treated_count = n()) %>% 
    filter(treated_count > 100)
  denominator_df <- 
    group_by_(a_df, facet_gp, "year_of_diagnosis") %>% 
    summarise(denominator = n())
  
  sumerian <- left_join(numerator_df, denominator_df) 
  
  sumerian[["treated_prop"]] <- 
    sumerian[["treated_count"]] / sumerian[["denominator"]]
  
  sumerian
}

plot_measure_in_this_df <- 
  function(this_df, which_cncr, the_measure_to_plot, predictors, outcomes) {
    map(predictors, function(the_predictor) {
      map(outcomes, function(the_outcome) {
        the_data <- 
          treated %>% 
          filter(site == which_cncr) %>% 
          summarise_by_year(the_predictor, the_outcome)
        if (!exists("summarised", .GlobalEnv))
          summarised <<- list()
        if (!(the_predictor %in% names(summarised))) 
          summarised[[the_predictor]] <<- list()
        
        summarised[[the_predictor]][[the_outcome]] <<- the_data
        message(paste("Plotting", 
                      which_cncr, 
                      the_measure_to_plot, 
                      the_predictor, 
                      the_outcome)
        )
        ggplot(the_data, 
               aes_string(x = "year_of_diagnosis", 
                          y = the_measure_to_plot, 
                          color = the_predictor)) + 
          #plot_by_radiation()
          geom_point() + 
          stat_smooth(geom = "line", method = "loess", n = 6, se = FALSE, alpha = 0.3) + 
          facet_wrap(reformulate(the_outcome), scales = "free_y", ncol = 2) + 
          labs(x = "", y = "Number treated", 
               title = paste(which_cncr, the_predictor, the_outcome, collapse = ", ")) + 
          scale_y_continuous(labels = scales::comma)
      }) %>% set_names(outcomes) 
    }) %>% set_names(predictors)
}


plot_nms <- function(figures) {
  map(names(figures), 
      function(f_nm) 
        map(names(figures[[1]]), 
            function(prd_nm)
              map_chr(names(figures[[1]][[1]]), 
                      function(otcm_nm)
                        paste(f_nm, prd_nm, otcm_nm, sep = "_"))) %>% 
        flatten_chr()) %>% 
    flatten_chr()
}

munge_figures <- function(site_specific_figures) {
  for (i in seq_along(figures)) {
    new_nms <- plot_nms(figures[[i]])
    figures[[i]] <- flatten(flatten(figures[[i]])) %>% set_names(new_nms)
  }
  figures
}

flatten_figures <- #function(figures) {
#  map(figures, 
  function(the_fs) {
    for (i in seq_along(the_fs)) {
      new_nms <- plot_nms(the_fs[[i]])
      the_fs[[i]] <- flatten(flatten(the_fs[[i]])) #%>% set_names(new_nms)
    }
    the_fs
  }#)
#}

write_rds_freal <- function(obj, nm, return_obj = TRUE) {
  if (exists("for_real")) 
    if (for_real)
      write_rds(obj, nm)
  if (return_obj) obj
}

write_csv_freal <- function(obj, nm, return_obj = TRUE) {
  if (exists("for_real")) 
    if (for_real)
      write_csv(obj, nm)
  if (return_obj) obj
}

one_hot <- function(a_df, to_oh) {
  lvls <- unique(a_df[[to_oh]])
  new_nms <- 
    reduce2(c("\\ ", 
             "[nN]one",
             "[uU]nknown"),
            c("_", 
             paste("None_", to_oh, sep = ""),
             paste("Unknown_", to_oh, sep = "")),
             function(x, y, z) gsub(y, z, x),
             .init = lvls) %>% 
         make.names(unique = T)
  
  hot_cols <- 
    map(lvls, function(lvl) ifelse(a_df[[to_oh]] == lvl, 1, 0)) %>% 
    set_names(new_nms) %>% 
    as.data.frame(stringsAsFactors = FALSE)
  
  bind_cols(a_df, hot_cols)
}

#' @export
one_hot_ncdb <- function(ncdb_df, to_one_hot) {
  prvs_nms <- list()
  for (i in seq_along(to_one_hot)) {
    prvs_nms[[to_one_hot[i]]] <- unique(ncdb_df[[to_one_hot[i]]])
    ncdb_df <- one_hot(ncdb_df, to_one_hot[i])
  }
  list(prvs_nms, ncdb_df)
}


tidy_biglm <- function(z) {
  as.data.frame(summary(z)[["mat"]]) %>%  
    modify_at(c("Coef", "(95%", "CI)"), ~ sprintf("%.2f", exp(.x))) %>% 
    select(estimate   = `Coef`, 
           lower_conf = `(95%`, 
           upper_conf = `CI)`, 
           p) %>%
    rownames_to_column(var = "term") %>%
    select(term, everything())
}



exclude_this_fn <- function() {
tidy_glm <- function(mdl) {
  if (is(mdl, "biglm")) return(tidy_biglm(mdl))
  if (all(mdl == "")) return(mdl)
  bind_cols(
    broom::tidy(mdl), 
    broom::confint_tidy(mdl)
  ) %>% 
    modify_at(c("estimate", "conf.low", "conf.high"), 
              ~ sprintf("%.5f", exp(as.numeric(as.character(.x))))) %>% 
    (function(the_df) {
      if (is(mdl, "glm")) {
        modify_at(the_df, "p.value", ~ sprintf("%.4f", .x))
      } else if (is(mdl, "speedglm")) {
        modify_at(the_df, 
          "p.value", 
          function(the_p) {
            parts <- strsplit(as.character(the_p), "e")
            sprintf("%.3f", as.numeric(format(the_p, scientific = FALSE)))
          }) 
      } else {
        browser()
        #stop("Can't tidy this type of model!")
      }
    })%>% 
    select(one_of(c("term", "estimate", "conf.low", "conf.high", "p.value"))) %>% 
    filter(term != "(Intercept)") %>% 
    (function(tidied_glm)
      if (any(grepl("\\*|:", tidied_glm[["term"]])))
        filter(tidied_glm, grepl("\\*|:"))
      else tidied_glm
    )
}
}

bad_model <- function(value) {
  function(..., term = value) {
    data.frame(term = term,
               estimate = value, 
               conf.low = value, 
               conf.high = value, 
               p.value = value, 
               stringsAsFactors = FALSE)
  }
}

set_select <- function(l, nms, new_vr) {
  new_vr <- quo_name(enquo(new_vr))
  stopifnot(length(l) == length(nms))
  set_names(l, nms) %>% 
  imap_dfr(~ mutate(.x, !!new_vr := .y)) %>% 
    select(!!new_vr, everything())
}

make_glm_models <- function(a_df, predictors, outcomes) {
  #browser()
  map(predictors, function(predictor) {
    map(outcomes, function(outcome) {
      tryCatch({
        if (any(table(a_df[,c(predictor, outcome)]) <= 11)) 
          bad_model("suppressed")(term = outcome)
        speedglm::speedglm(reformulate(predictor, outcome), a_df) %>% 
          tidy_glm()
        },
        # error = tryCatch({
        #   biglm::bigglm(glm(reformulate(predictor, outcome), 
        #                     a_df, 
        #                     family = binomial()
        #                     )
        #                 ) %>% 
        #     tidy_glm()
        # }, 
        error = bad_model("bad") 
        #)
        )
    }) %>% 
    set_select(outcomes, outcome)
  }) %>% 
    set_names(predictors) %>% 
    imap_dfr(function(otcms_df, prd) {
      otcms_df[["term"]] <- 
        gsub(paste("^", prd, sep = ""), "", otcms_df[["term"]])
      mutate(otcms_df, predictor = prd) %>% 
      select(outcome, predictor, everything())
    })
}



m_fn <- function(a_l, otcm, filter_expr) {
  filter_expr <- enquo(filter_expr)
  map(a_l, 
      function(l_element) l_element %>% 
        filter(outcome == otcm & !!filter_expr) %>% 
        modify_at(c("estimate", "conf.low", "conf.high", "p.value"), as.numeric) %>% 
        ggplot(aes(x = term, y = estimate, color = brain_mets)) + 
        geom_errorbar(aes(ymin = conf.low, ymax = conf.high), 
                      position = position_dodge(width = 0.75)) + 
        geom_point(position = position_dodge(width = 0.75)) + 
        facet_wrap(~ predictor, scales = "free") + 
        theme(axis.text.x = element_text(angle = 45, hjust = 1))
  )
}

make_plots <- function(a_l, plot_outcomes, an_expr) {
  an_expr <- enquo(an_expr)
  
  map(plot_outcomes, ~ m_fn(a_l, .x, !!an_expr)
      ) %>% 
    set_names(plot_outcomes)
}

tidy_reference <- function(outcome, predictor, term) {
  data.frame(outcome    = outcome,
             predictor  = predictor,
             term       = term,
             estimate   = 1,
             conf.low   = 1,
             conf.high  = 1,
             p.value    = 1,
             brain_mets = c("Present", "Absent"), 
             stringsAsFactors = FALSE)
}


get_models_things <- function(models, thing) {
  map(models, ~ unique(.x[[thing]])) %>% 
    flatten_chr() %>% 
    unique()
}

reference_row <- function(otcm, prdctr) {
  data.frame(outcome = otcm, 
             predictor = prdctr,
             term = "2010", 
             estimate = "1.00", 
             conf.low = "1.00", 
             conf.high = "1.00", 
             p.value = "1.000", 
             stringsAsFactors = FALSE)
}

bind_rows_mutate <- function(df_l, newcol_name) {
  if (is_null(names(df_l))) stop("This function requires a named list")
  imap_dfr(df_l, 
          ~ mutate(.x, 
                   !!newcol_name := .y))
}



