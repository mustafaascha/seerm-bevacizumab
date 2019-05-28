
library(tidyverse)

devtools::load_all("mabsr")

chemo <- 
  read_csv("cache/chemotherapy.csv.gz") %>% 
  left_join(mabsr::drug_names) %>% 
  select(-drug) %>% 
  rename(drug = drugs_short) %>% 
  modify_at("drug", ~ ifelse(grepl("platin", .x), "platin", .x))

if (!exists("cancers")) {
  cancers <- 
    read_csv("cache/cancers-analytic.csv.gz", col_types = c_cols()) %>% 
    reg_df() %>% 
    srv_df()
}

chemo <- munge_merged(chemo, cancers)

write_rds(chemo, "cache/chemo-analytic.rds")


