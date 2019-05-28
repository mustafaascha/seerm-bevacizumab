
library(tidyverse)

devtools::load_all("augur")

px_codes <- readLines("documentation/bev-codes.txt")

out_cpt_img <- 
  extract_slice_codes("seerm", "outsaf", 
                      augur::infiles[["outpat"]], 
                      px_codes, "hcpcs") %>% 
  bind_rows()
write_csv(out_cpt_img, "cache/mabs/cpt-mb-out.csv.gz")







