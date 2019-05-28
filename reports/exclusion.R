
library(tableone)
library(tidyverse)
devtools::load_all("augur")
devtools::load_all("frequencies")
devtools::load_all("mabsr")

cancers <- 
  read_csv("cache/cancers_before_exclusion.csv.gz", progress = FALSE) %>% 
  filter(dx_year >= 2010 & which_cancer == "lung")

paper_products <- list()

paper_products[["not_nscl"]] <- 
  filter(cancers,
         (hist03v %in% mabsr::histo_exclusion)) %>% 
  (function(z) data.frame(table(z$siterwho, z$hist03v)))
paper_products[["stage_iiib_sbm"]] <- 
  data.frame(table(cancers$d_ajcc_s_v, cancers$hist03v, cancers$csmetsdxbr_pub_v))

cancers <- 
  filter(cancers,
         siterwho == 22030 & 
         !(hist03v %in% mabsr::histo_exclusion) & 
         ! (d_ajcc_s_v == "Stage IIIB" & csmetsdxbr_pub_v == "Yes"))

expected_whocodes <- c(22030, 25010, 26000)
expected_cancers <- 
  map(expected_whocodes, ~ table(cancers$siterwho == .x))
names(expected_cancers) <- paste("x", expected_whocodes, sep = "")

#  we're working with lung, skin, and breast cancers
#  is this necssaru pr correct?                 
paper_products[["not_right_cancer"]] <- expected_cancers

#fix NAs so they're correctly represented
cancers[,grep("csmetsdx", names(cancers))] <- 
  lapply(cancers[,grep("csmetsdx", names(cancers))], 
        function(x) ifelse(x == "Unknown" | x == "Not applicable", NA, x)
        )

#cancers$dx_year <- factor(cancers$dx_year)
cancers$dx_year_c <- cut(as.numeric(cancers$dx_year), c(1975, 2007:2012))

# exclusions table ==========================================
# nobody had missing or unknown sex                          
# nobody had missing or unknown date of death                
#   see: table(cancers$ser_dody, cancers$stat_rec)

vfe <- vars_for_exclusion <-
    c("hmocnt", "age_dx", "dx_year", "typefup",
      "med_stcd", "med_stcd", "med_stcd", 
      "med_stcd", "race",  "srv_time_mon",
     #"hmocnt1991_2015", "age_dx_e",
      "vrfydth", "payerdx_v",
      "payerdx_v"
      #, "all_same_cancer"
      )
vte <- values_to_exclude <-
    c("Used HMO", "Under 65 years", "Before 2008",
      "Autopsy/Death Certificate",
      "ESRD only", "Aged with ESRD", "Disabled",
      "Disabled with ESRD", "Unknown", "9999",
     #"Used HMO", "Age < 65", 
      "No", "Insurance status unknown",
      "Not insured"
      #, "Unaccounted primaries present"
      )
exclusion_df_vars <- 
  c("patient_id", "which_cancer", 
    names(cancers)[grep("_e\\.?[0-5]?$", names(cancers))])
exclusion_df <- cancers[,exclusion_df_vars]

medicare_eligibility_vars <- 
  c("hmocnt_e", "age_dx_e",
    "med_stcd_v_e", "med_stcd_v_e.1", "med_stcd_v_e.2",
    "med_stcd_v_e.3", "payerdx_v_e", "payerdx_v_e.1")

medicare_exclusions <- 
  exclusion_df[,c("patient_id", "which_cancer", medicare_eligibility_vars)]
medicare_exclusions[["any_medicare"]] <- 
  apply(medicare_exclusions[,3:10], 1, any)
medicare_exclusions_summary <- 
  medicare_exclusions %>% 
  group_by(which_cancer, any_medicare) %>% 
  summarise(ex_med = n())

general_exclusions <- 
  exclusion_df[,!(names(exclusion_df) %in% medicare_eligibility_vars)]
general_exclusions[["any_general"]] <- 
  apply(general_exclusions[,3:7], 1, any)
general_exclusions_summary <- 
  general_exclusions %>% 
  group_by(which_cancer, any_general) %>% 
  summarise(ex_gen = n()) 

m_bind <- medicare_exclusions_summary %>% 
    filter(any_medicare) %>% select(-any_medicare)
g_bind <- general_exclusions_summary %>% 
    filter(any_general) %>% select(-any_general)

paper_products[["selected_ns"]] <- 
  reduce(list(m_bind, g_bind), left_join)


exclusion_map <- 
  exclusion_df %>% group_by(which_cancer) %>% 
# distinct(patient_id, .keep_all = TRUE) %>%
  select(-patient_id) %>%  nest
exclusion_map$n <- map_int(exclusion_map$data, ~ nrow(.x))

exclusion_sums <- 
  map_dfr(exclusion_map$data, 
          function(exclusion_df) {
          map_df(exclusion_df, 
                 function(exclusion_column){ 
                   sum(exclusion_column, na.rm = TRUE)
                 })
          })

#transpose and return to dataframe with variable column
#backup_em <- exclusion_map
exclusion_map <- bind_cols(exclusion_map, exclusion_sums) %>% select(-data)
exclusion_map <- data.frame(t(exclusion_map), stringsAsFactors = FALSE)
names(exclusion_map) <- exclusion_map[1,]
exclusion_map <- exclusion_map[-1,]
exclusion_map$Var <- rownames(exclusion_map)
rownames(exclusion_map) <- NULL
exclusion_map$Vars <- c("n", vfe)
exclusion_map$Values <- c("", vte)

paper_products[["exclusions"]] <- exclusion_map

write_rds(paper_products, "results/exclusion.rds")



