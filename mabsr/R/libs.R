load_libs <- function() {
  libraries <- 
    c(
      "speedglm",
      "biglm",
      "glmnet",
      "MatchIt",
      "parallel",
      "doMC",
      "tidyverse",
      #"fiftystater",
      "rlang",
      "tableone",
      "zeallot",
      "mice",
      "lmerTest",
      "missRanger",
      "survival",
      "lubridate",
      "survminer",
      "lars"#, "biglars", "leaps"
    )
  for (l in libraries) {
    if (!(l %in% installed.packages())) install.packages(l)
    library(l, character.only = T)
  }
}
