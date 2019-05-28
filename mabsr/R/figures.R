      
do_srvfit <- function(a_df, prd) {
  the_frm <- paste('Surv(srv_time, event = alive == "Dead")~', prd)
  do.call(survfit, 
          list(formula = as.formula(the_frm), 
               data = a_df))
}

do_plt <- function(a_fit, a_df, ...) {
  the_params <- 
    list(a_fit, 
         data = a_df, 
         break.time.by = 12, 
         #pval = TRUE,
         #pval.coord = c(0.8, 0.8),
         surv.median.line = "hv",
         conf.int = TRUE,
         # Add risk table
         risk.table = TRUE,
         tables.height = 0.2,
         tables.theme = theme_cleantable(),
         # Color palettes. Use custom color: c("#E7B800", "#2E9FDF"),
         # or brewer color (e.g.: "Dark2"), or ggsci color (e.g.: "jco")
         palette = "jco",
         ggtheme = theme_bw(),# Change ggplot2 theme, 
         ...
         )
  do.call(ggsurvplot, the_params)
}

bm_plt <- function(the_df) {
  the_df$B <- 
    ifelse(the_df$bmab01 == 1, "Bevacizumab", "None")
  the_fit <- do_srvfit(the_df, "B")
  do_plt(the_fit, the_df)
}

sbm_plt <- function(the_df) {
  the_df$SBM <- 
    ifelse(the_df$seer_br_mets == "SEER_Positive", "Present", "Absent")
  the_fit <- do_srvfit(the_df, "SBM")
  do_plt(the_fit, the_df)
}


switch_plt <- function(the_df, the_vr) {
  switch(the_vr, 
    bmab01 = {
      the_df$B <- 
        ifelse(the_df$bmab01 == 1, "Bevacizumab", "None")
        the_fit <- do_srvfit(the_df, "B")
      },
    seer_br_mets = {
      the_df$SBM <- 
        ifelse(the_df$seer_br_mets == "SEER_Positive", "Present", "Absent")
        the_fit <- do_srvfit(the_df, "SBM")
      }
    )
  do_plt(the_fit, the_df)
}

ggsurv_ <- function(frm, the_df, ...) {

  params <- 
    list(
      fit              = surv_fit(frm, data = the_df), 
      #data           = the_df, 
      break.time.by    = 12, 
      conf.int         = TRUE,
      risk.table       = "percentage",
      tables.height    = 0.2,
      palette          = "grey",
      tables.theme     = theme_cleantable(),
      ggtheme          = theme_bw(), 
      ylab             = "Survival Probability",
      legend.labs      = c("No Bev", "Bev"),
      linetype       = c("solid", "dashed"),
      surv.median.line = "v",
      censor = FALSE,
      ...
      )
  #browser()
  do.call("ggsurvplot", params) + 
    labs(color = "Treatment", fill = "Treatment", x = "Time")
}

ggsurvadj_ <- function(frm, the_df, ...) {
  params <- 
    list(
      fit            = do.call("coxph", list(formula = frm, data = the_df)), 
      data           = the_df, 
      variable       = "bmab01",
      xticks.by      = 12, 
      conf.int       = TRUE, 
      ggtheme        = theme_bw(), 
      legend.title   = "Treatment", 
      legend.labs    = c("No Bev", "Bev"),
      ylab           = "Survival Probability", 
      xlim           = c(0, 60), 
      legend         = "none", 
      method         = "average", 
      censor         = FALSE,
      linetype       = c("solid", "dashed"),
      ...
    )

  (do.call("ggadjustedcurves", params) %>%
    set_palette("grey")) + 
    #set_palette("jco")) + 
    labs(color = "Treatment", x = "Time")
}




