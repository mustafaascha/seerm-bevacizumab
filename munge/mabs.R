

library(tidyverse)

if (!exists("cancers")) {
 cancers <- 
    read_csv("cache/cancers_postrecode.csv.gz")
}

these <- 
  c("patient_id", 
    "from_dtm", "from_dtd", "from_dty",
    "thru_dtm", "thru_dtd", "thru_dty",
    "hcpcs")

these_cols <- 
  cols_only(
    patient_id = col_character(),
    from_dtm = col_character(),
    from_dtd = col_character(),
    from_dty = col_character(),
    thru_dtm = col_character(),
    thru_dtd = col_character(),
    thru_dty = col_character(),
    hcpcs = col_character()
  )

mab_codes <- 
  data.frame(
    hcpcs = 
      c("S0116",  "J9035", "C9257", "C9214", "Q2024"),
    mab = rep_len("Bevacizumab", 5), 
    stringsAsFactors = FALSE
  )

# should specify one of as col_types
mabs <- 
  map_dfr(list.files("cache/mabs", full.names = TRUE), 
      function(fl) read_csv(fl, col_types = these_cols)
        #select(read_csv(fl, col_types = cols(.default = "c")), 
               #one_of(these))
  ) %>% 
  left_join(mab_codes) %>% 
  select(-hcpcs) %>% 
  distinct() %>% 
  mutate(
    from_dmy = paste(from_dtd, from_dtm, from_dty, sep = "-"),
    thru_dmy = paste(thru_dtd, thru_dtm, thru_dty, sep = "-")
  ) %>% 
  modify_at(c("from_dmy", "thru_dmy"), lubridate::dmy) %>% 
  mutate(drug_dur = as.numeric(thru_dmy - from_dmy)) %>% 
  group_by(patient_id, mab) %>% 
  summarise(
    mab_strt = min(from_dmy, na.rm = T),
    mab_time = sum(drug_dur, na.rm = T), 
    n_mab_courses = n()
  ) %>% 
  modify_at("mab_time", ~ ifelse(.x <= 0, 1, .x))
  
cancers <- 
  left_join(cancers, mabs) %>% 
  modify_at("mab", ~ ifelse(is.na(.x), "None", .x))



















write_csv(cancers, "cache/cancers_w_mabs.csv.gz")


