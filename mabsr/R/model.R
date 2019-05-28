
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




