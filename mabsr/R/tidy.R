
empty_model_row <- function(z) {
  if (!missing("z")) term <- z else term <- ""
  data.frame(term = z, 
             estimate = "", conf.low = "", 
             conf.high = "", p.value = "", 
             stringsAsFactors = FALSE)
}

too_few_in_cats <- function(prd) {
  if (is.character(prd)) return(any(table(prd) <= 5))
  else return(FALSE)
}

twocol_mdl <- function(mdl_df, lbl) {
  mdl_df %>% 
    mutate(!!lbl := paste(estimate, " (", 
                          conf.low, "-", 
                          conf.high, ", p: ", 
                          p.value, ")", 
                          sep = "")) %>% 
    select(term, !!lbl)
}


cphw_pretidy <- function(the_mdl) {
    wm_coefs <- coef(the_mdl)
    wm_ci <- confint(the_mdl)
    #browser()
    to_return <- 
      data.frame(term = names(wm_coefs), 
                 estimate = (wm_coefs), 
                 conf.low = (wm_ci[,1]),
                 conf.high = (wm_ci[,2]), 
                 p.value = the_mdl[["prob"]], 
                 stringsAsFactors = FALSE)
    rownames(to_return) <- NULL
    to_return
}

add_binary_metrics <- function(the_mdl, tidied_mdl) {
  predictions <- predict(the_mdl, type = "response")
  outcome_values <- the_mdl[["y"]] 
  tidied_mdl[["h_l"]] <- 
   ResourceSelection::hoslem.test(
                                  outcome_values, predictions, g = 5 
                                 )[["p.value"]] 
  tidied_mdl[["c_s"]] <- pROC::auc(pROC::roc(outcome_values ~ predictions))
  tidied_mdl
}

tidy_mdl <- function(the_mdl) {
  if ("coxphw" %in% class(the_mdl)) {
    pretidied_mdl <- cphw_pretidy(the_mdl)
  } else {
    pretidied_mdl <- 
      quietly(bind_cols)(
                        broom::confint_tidy(the_mdl),
                        broom::tidy(the_mdl)
                        )[["result"]]
  }

  tidied_mdl <- 
      pretidied_mdl %>% 
      modify_at(c("estimate", "conf.low", "conf.high", "p.value"), 
                compose(as.numeric, as.character)) %>% 
      modify_at(c("estimate", "conf.low", "conf.high"),
               function(z) sprintf("%.2f", exp(z))) %>% 
      modify_at("p.value", partial(sprintf, fmt = "%.4f")) %>%
      dplyr::select(term, estimate, conf.low, conf.high, p.value) %>%
      filter(term != "(Intercept)") 
                            
  if ("glm" %in% class(the_mdl)) {
    tidied_mdl <- 
      add_binary_metrics(the_mdl = the_mdl, tidied_mdl = tidied_mdl)
  } else if ("coxph" %in% class(the_mdl)) {
    # global side effect: check proportional hazards assumption
    #if(!exists("cox_zphs")) assign("cox_zphs", list(), envir = .GlobalEnv)
    #cox_zphs[[as.character(deparse(the_mdl$formula))]] <<- 
    #  cox.zph(the_mdl)  
  } 

  tidied_mdl
}


