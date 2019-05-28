
library(tidyverse)

devtools::load_all("augur")

ptids <- 
  read_csv("cache/cancers.csv.gz", 
           col_types = cols_only(which_cancer = col_character(), 
                                 patient_id = col_character())
  ) %>%
  select(patient_id) %>% 
  unlist()

#nch_bev <- nch_w_this_ptid("seerm", ptids)

filetype <- 
  c(
   # "nch",  
   # "dme",
    "CCflag" 
   # "medpar", 
   # "outsaf", 
   # "pdesaf"
    )

these_files_w_these_ptids(filetype, "seerm", ptids) %>% 
  write_csv(path =
    paste(
      paste("cache/bevpts/bev-patients-", filetype, sep = ""), 
      ".csv.gz", sep = "")
  )


