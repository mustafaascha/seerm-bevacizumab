
#' @import glmnet
#' @import dplyr
#' @import purrr
impute <- function(a_df, predictors, outcome = "bmab01", ...) {
  library(mice)
  get_missing_props <- function(a_df) {
    map_dbl(a_df, ~ length(which(is.na(.x))) / length(.x))
  }
  a_df <- 
    a_df[,c(predictors, outcome, "patient_id")] %>% 
    mutate_all(function(z) ifelse(z %in% c("Missing", "Unknown"), NA, z))

  outcome_col <- 
    length(c(predictors, outcome)) %>% 
    modify_if(is.factor, as.character)

  mice_out <- mice::mice(a_df, ...)
  mice_df <- mice::complete(mice_out)
  matrix_df <- 
      drop_na(mice_df)

  missing_prop <- 
    list(missing_before = get_missing_props(a_df), 
         missing_after = get_missing_props(mice_df))

  list(imputed = matrix_df,
       missing = missing_prop)
}

cv_glm_lasso <- function(xs, y, ...) {
  #message("Organize so that selected predictors are automatically selected for mv models")
  #browser()A
  ls <- 10 ^ seq(4, -4, by = -.1) 
  y <- data.matrix(y)
  xs <- data.matrix(xs)
  mdl <- glmnet(x = xs, y = y, family = "binomial")
  c_val <- cv.glmnet(x = xs, y = y, lambda = ls, nfolds = 2e2, ...) 
  list(mdl, c_val)
}

cv_cox_lasso <- function(cox_df, prds, primary_prd = "bmab01", ...) {
  #message("Organize so that selected predictors are automatically selected for mv models")
  ls <- 10 ^ seq(4, -4, by = -.1) 
  cox_df <- filter(cox_df, srv_time != 0)
  y <- as.matrix(cbind(time = cox_df[["srv_time"]], status = (cox_df[["alive"]] == "Dead")))
  xs <- data.matrix(map_df(select(cox_df, one_of(prds)), as.factor))
  #browser()
  ps <- 
    list(
      x = xs, 
      y = y, 
      family = "cox", 
      penalty.factor = colnames(xs) != primary_prd, 
      ...
    )
  do.call("cv.glmnet", ps)
}

rm_p <- function(ps, p) ps[-grep()]

do_glm <- function(a_df, predictors, outcome = "bmab01") {
  prds <- predictors
  message(paste("Making glm for", paste(predictors, collapse = ", ")))
  #library(speedglm)
  #tryCatch({
   #the_mdl <- speedglm(reformulate(predictors, outcome), data = a_df) 
   #if (length(prds) > 1) {
   #   the_mdl <- speedglm(reformulate(prds, outcome), data = a_df) 
   #} else { 
      the_mdl <- glm(reformulate(prds, outcome), data = a_df, family = "binomial") 
   #}
   tidy_glm(the_mdl)
  #}, 
  #error = empty_model_row
  #)
}

do_cox <- function(a_df, predictors, wph = FALSE, mdls = FALSE) {
  library(survival)
  library(coxphw)
  the_formula <- 
    as.formula(paste('Surv(srv_time, event = alive == "Dead")~', 
               paste(predictors, collapse = "+")))
  message(paste("Making cox model for", paste(predictors, collapse = ", ")))

  tryCatch({ 
    # use weighted estimation for non-proportional hazards
    #if (wph & ("bmab01" %in% predictors)) {
    if (wph) {
      a_df <- arrange(a_df, srv_time)
      #ctrl <- coxphw.control(add.constant = 0.1, gconv = 0.001)
      the_mdl <- 
        coxphw(the_formula, 
                     data = a_df, 
                     #control = ctrl,
                     robust = FALSE, 
                     trunc.weights = 0.95,
                     sorted = TRUE
                    )
      return_v <- tidy_mdl(the_mdl)
      the_mdl[c("linear_predictors", "caseweights", "w.matrix", "y", "x")] <- NULL
    } else {
      the_mdl <- coxph(the_formula, data = a_df)              
      return_v <- tidy_mdl(the_mdl) 
    }
    #only keep the interaction term
    if (any(grepl(":", return_v[["term"]]))) {               
      return_v <- filter(return_v, grepl(":", term))        
    }
    if (mdls) return(list(the_mdl, return_v)) else return(return_v)
    }, 
    error =  function(e) {
      message(paste("there was a problem processing", 
              paste(predictors, collapse = ","), 
              ":", e))
      #browser()
      empty_model_row(paste(predictors, collapse = "_"))
    })                                                       
}

do_prop_mdls <- function(params) {
  # USE LIFT?
  #mtch <- 
  #  ifelse("propensity_score" %in% names(params[["data"]]), "matched", "unmatched")
  mdl <- 
    do_cox(
      params[["data"]], 
      params[["prds"]], 
      wph = params[["wph"]], 
      mdls = params[["mdls"]]
    )

  #if (params[["wph"]]) lbl <- "AHR" else lbl <- "HR"
  # Avoid indexing by integers--prefer strings
  model_df <- 
    mdl[[2]] %>% 
      twocol_mdl(lbl = "HR") %>% 
      mutate(matched = attributes(params[["data"]])[["matched"]]) %>%
      modify_at("term", function(x) gsub("seer_br_metsSEER_Positive", "Brain Mets", x)) %>%
      modify_at("term", function(x) gsub("bmab01", "Bevacizumab", x)) %>% 
      mutate(weighted = ifelse(params[["wph"]], "weighted", "not weighted"))

  if (any(grepl("\\:", model_df$term))) model_df <- filter(model_df, grepl("\\:", term))
  
  if (params[["mdls"]]) return(list(mdl, model_df)) else return(model_df)
}

#' @import dplyr
gp_count <- function(df, ...){
  gps <- quos(...)
  to_return <- df %>% group_by(!!!gps) %>% summarise(cnt = n())
  to_return[["cnt"]][to_return$cnt <= 11 & to_return$cnt > 0] <- NA
  to_return
}


select_coefs <- function(net_cv) UseMethod("select_coefs", net_cv)

select_coefs.cv.glmnet <- function(net_cv) {
  coef(net_cv, s = "lambda.min") %>%   
    as.matrix() %>%
    as.data.frame() %>%
    rownames_to_column() %>%
    set_names(c("v", "s")) %>% 
    filter(s != 0) %>%
    select(v) %>%
    unlist() %>% 
    unname() %>% 
    grep("\\(", ., value = T, invert = T)
}


select_coefs.coxnet <- function(net_cv) {
  #browser()
  coef(net_cv, s = net_cv$lambda.min) %>%   
    as.matrix() %>%
    as.data.frame() %>%
    rownames_to_column() %>%
    set_names(c("v", "s")) %>% 
    filter(s != 0) %>%
    select(v) %>%
    unlist() %>% 
    unname() %>% 
    grep("\\(", ., value = T, invert = T)
}

f_prds <- function(frm) {
  strsplit(as.character(frm)[[3]], split = "\\ \\+\\ ") %>% unlist()
}

count_bmab_prds <- function(frm, unmatched, matched) {
  map(f_prds(frm), function(prd) {                           
    map(list(unmatched = unmatched, matched = matched), 
        ~ table(.x[["bmab01"]], .x[[prd]]))                       
  })  %>% 
  set_names(f_prds(frm))
}


evaluate_propensity_scoring <- function(frm, unmatched, matched) {
  list(
    ps_plots  = map(list(unmatched, matched), ggplot),
    ps_counts = count_bmab_prds(frm, unmatched, matched)
  )
}

add_matched <- function(obj, v) {
 attr(obj, "matched") <- v
 obj
}

bm_pos <- function(a_df) {
  a_df[with(a_df, brain_mets == "SEER_Positive"),]
}

postprocess_match <- function(a_df, frm, nm) {
  a_df[["propensity_score"]] <- predict(glm(frm, data = a_df))
  attr(a_df, "matched") <- nm
  a_df %>% modify_at("srv_time_mon", as.numeric)
}




srv2mbs <- function(a_df, cancers_df) {
  a_df <- 
    cancers_df %>% 
    select(patient_id, dx_year, dx_month, mab_start) %>% 
    left_join(a_df, by = "patient_id")
  if (!("srv_time_mab" %in% names(a_df))) {
    a_df <- srv_df(a_df, mb = TRUE)
  } 
  a_df[["srv_time"]] <- a_df[["srv_time_mab"]]
  a_df[["srv_time_mon"]] <- a_df[["srv_time_mab"]]
  a_df
}



