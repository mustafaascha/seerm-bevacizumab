
library(tidyverse)

devtools::load_all("augur")

code_table <- with(augur::dissertation_codes, 
                   table(Treatment.Type, Code.Classification))

px_codes <- 
  augur::dissertation_codes %>% 
  filter(Treatment.Type == "Neurosurgical resection" & 
         (Code.Classification == "CPT"  | Code.Classification == "HCPCS"))

outp <- 
  extract_slice_codes("seerm", "outsaf", augur::infiles[["outpat"]], px_codes, "hcpcs") %>% 
  bind_rows
write_csv(outp, "cache/neurosurg/cpt-nsrg-out.csv.gz")







