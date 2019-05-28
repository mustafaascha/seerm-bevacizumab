
library(tidyverse)

devtools::load_all("augur")

px_codes <- readLines("documentation/bev-codes.txt")

dme_cpt_mab <- 
  extract_slice_codes("seerm", "dme", augur::infiles[["dme"]], px_codes, "hcpcs") %>% 
  bind_rows
write_csv(dme_cpt_mab, "cache/mabs/cpt-mb-dme.csv.gz")







