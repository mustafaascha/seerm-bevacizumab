


classification_metrics <- function(dflist, dfnm) {
  library(pROC)
  library(ResourceSelection)
  map(dflist, function(mdl) {
    predictions <- predict(mdl[[1]], type = "response")
    the_outcome <- mdl[[2]]
    outcome_values <- mdl[[1]][["data"]][[the_outcome]]
    test <- hoslem.test(outcome_values, predictions, g = 5)
    h_l <- hoslem.test(outcome_values, predictions, g = 5)[["p.value"]]
    c_s <- auc(roc(outcome_values ~ predictions))
    data.frame(c_stat = c_s, 
               h_l_stat =  h_l, 
               glm = dfnm, 
               outcome = the_outcome, 
               stringsAsFactors = FALSE
    )
  })
}

#' Rubin rules
#'
#' @param ps_glm A model constructed using propensity-score matching/adjustment 
#' @param u_df An unmatched data frame
#' @param m_df A matched data frame
#'
#' @return
#' @export
#'
#' @examples
rubin_rules <- function(u_df, m_df, frmla) {
  ps_glm <- 
    glm(frmla, 
        data = select(m_df, -distance, -weights), 
        family = "binomial")
  u_df$pscore <- predict(ps_glm, newdata = u_df)
  m_df$pscore <- predict(ps_glm, newdata = m_df)
  rule_1 <- function(df) {
    with(df, 
         abs(100 * (mean(pscore[Insurance == 0], na.rm=T) - 
                      mean(pscore[Insurance == 1], na.rm=T)) /
               sd(pscore, na.rm = T)))
  }
  rule_2 <- function(df) {
    with(df, var(pscore[Insurance == 0], na.rm = T) / 
           var(pscore[Insurance == 1], na.rm = T))
  }
  list(
    rule_1 = list(
      unmatched = rule_1(u_df),
      matched = rule_1(m_df)
      ),
    rule_2 = list(
      unmatched = rule_2(u_df),
      matched = rule_2(m_df)
    )
  )
}






