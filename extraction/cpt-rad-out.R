
library(tidyverse)

devtools::load_all("augur")

code_table <- with(augur::dissertation_codes, 
                   table(Treatment.Type, Code.Classification))

px_codes <- 
  augur::dissertation_codes %>% filter(Treatment.Type == "Radiation Therapy")

out_cpt_rad <- 
  extract_slice_codes("seerm", "outsaf", augur::infiles[["outpat"]], px_codes, "hcpcs") %>% 
  bind_rows
write_csv(out_cpt_rad, "cache/rad/cpt-rad-out.csv.gz")







