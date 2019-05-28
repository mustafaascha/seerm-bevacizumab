
library(tidyverse)

devtools::load_all("augur")

px_codes <- 
  paste("01", c(21:25, 31, 51, 59), sep = "")

mpr_icd_dx <- 
  extract_slice_codes("seerm", "medpar", 
                      augur::infiles[["medpar"]], 
                      px_codes, 
                      "SRGCDE") %>% 
  bind_rows
write_csv(mpr_icd_dx, "cache/neurosurg/icd-nsrg-mpr.csv.gz")








