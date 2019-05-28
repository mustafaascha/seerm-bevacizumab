
library(tidyverse)

devtools::load_all("augur")

code_table <- with(augur::dissertation_codes, 
                   table(Treatment.Type, Code.Classification))

px_codes <- 
  augur::dissertation_codes %>% 
  filter(Treatment.Type == "Neurosurgical resection" & 
         (Code.Classification == "CPT"  | Code.Classification == "HCPCS"))

nch <- 
  extract_slice_codes("seerm", "nch", augur::infiles[["nch"]], px_codes, "hcpcs") %>% 
  bind_rows
write_csv(nch, "cache/neurosurg/cpt-nsrg-nch.csv.gz")





