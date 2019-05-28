
#cancers <- add_claims(cancers, read_matches("cache/stereo"), "stereo")

if (!exists("cancers")) {
 cancers <- 
    read_csv("cache/cancers_postrecode.csv.gz")
}

these <- 
  c("patient_id", 
    "from_dtm", "from_dtd", "from_dty",
    "thru_dtm", "thru_dtd", "thru_dty",
    "hcpcs")

#dflt_c <- 
#  cols_only(
#    patient_id = col_character(),
#    from_dtm = col_character(),
#    from_dtd = col_character(),
#    from_dty = col_character(),
#    thru_dtm = col_character(),
#    thru_dtd = col_character(),
#    thru_dty = col_character(),
#    hcpcs = col_character()
#  )
#
#the_cols <- 
#  list(dflt_c, 
#  cols_only(
#    patient_id = col_character(),
#    from_dtm = col_character(),
#    from_dtd = col_character(),
#    from_dty = col_character(),
#    thru_dtm = col_character(),
#    thru_dtd = col_character(),
#    thru_dty = col_character(),
#    hcpcs = col_character()
#  ), 
#  dflt_c)

str_cols <- 
  c("patient_id", "from_dtm",
  "from_dtd", "from_dty",
  "thru_dtm", "thru_dtd",
  "thru_dty", "hcpcs")

# should specify one of as col_types
strs <- 
  map(list.files("cache/stereo", full.names = TRUE), read_csv) %>% 
  map_dfr(function(z) {
    names(z) <- tolower(names(z))
    z <- select(z, one_of(str_cols))
    z[] <- lapply(z, as.character)
    z
  }) %>% 
  distinct() %>% 
  mutate(
    from_dmy = paste(from_dtd, from_dtm, from_dty, sep = "-"),
    thru_dmy = paste(thru_dtd, thru_dtm, thru_dty, sep = "-")
  ) %>% 
  modify_at(c("from_dmy", "thru_dmy"), lubridate::dmy) %>% 
  mutate(str_dur = as.numeric(thru_dmy - from_dmy)) %>% 
  group_by(patient_id) %>% 
  summarise(
    str_start = min(from_dmy, na.rm = T),
    str_durat = sum(str_dur, na.rm = T), 
    n_str_crs = n()
  ) %>% 
  modify_at("str_durat", ~ ifelse(.x <= 0, 1, .x)) 
  
cancers <- left_join(cancers, strs, by = "patient_id") 



