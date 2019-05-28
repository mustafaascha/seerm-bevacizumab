
library(tidyverse)

devtools::load_all("augur")

px_codes <- 
  paste("01", c(21:25, 31, 51, 59), sep = "")

out_icd_dx <- 
  extract_slice_codes("seerm", "outsaf", 
                     augur::infiles[["outpat"]], 
                     px_codes, paste("prcdr_cd", 1:13, sep = "_")
                     ) %>% 
  bind_rows()

write_csv(out_icd_dx, "cache/neurosurg/icd-nsrg-out.csv.gz")







