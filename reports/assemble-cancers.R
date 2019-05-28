
cancers <- 
  read_csv("cache/cancers-analytic.csv.gz") %>% 
  exclude() %>% 
  left_join(read_chemo("cache/chemotherapy.csv.gz")) %>% 
  select(-d_ssg00) %>% 
  mutate(d_ssg00 = d_ajcc_s_v) %>% 
  modify_at(
     c("n_carboplatin", "n_cisplatin", 
     "n_dexamethasone", "n_paclitaxel", "n_pemetrexed",
     "Stereo_RS", "NSRG", "RAD"
     ),
    ~ ifelse(is.na(.x), "None", .x)
  )
source("reports/add-mabs.R")
source("reports/add-stereo.R")
write_csv(cancers, "cache/cancers-bev.csv.gz")
