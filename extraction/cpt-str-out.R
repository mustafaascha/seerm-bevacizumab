
library(tidyverse)

devtools::load_all("augur")

code_table <- with(augur::dissertation_codes, 
                   table(Treatment.Type, Code.Classification))

px_codes <- 
  augur::dissertation_codes %>% filter(Treatment.Type == "Stereotactic radiosurgery")

out_cpt_stereo <- 
  extract_slice_codes("seerm", "outsaf", augur::infiles[["outpat"]], px_codes, "hcpcs") %>% 
  bind_rows
write_csv(out_cpt_stereo, "cache/stereo/cpt-str-out.csv.gz")







