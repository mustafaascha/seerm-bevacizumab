
library(tidyverse)

devtools::load_all("augur")

px_codes <- readLines("documentation/bev-codes.txt")

nch_cpt_mab <- 
  extract_slice_codes("seerm", "nch", 
                      augur::infiles[["nch"]], 
                      px_codes, "hcpcs") %>% 
  bind_rows
write_csv(nch_cpt_mab, "cache/mabs/cpt-mb-nch.csv.gz")







