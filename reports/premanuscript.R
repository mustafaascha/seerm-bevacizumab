
library(tidyverse)
library(tableone)
library(zeallot)
library(devtools)
devtools::load_all("augur")
devtools::load_all("frequencies")
devtools::load_all("mabsr")

paper_products <- list()

if(!exists("cancers")) {
  source("reports/load_exclude.R")
}

#-------------------------------------------------

cancers <- filter(cancers, dx_year >= 2010)


paper_products[["df_counts"]] <- 
             cancers %>% 
               group_by(which_cancer) %>% 
               summarise(count = n()) %>% 
               mutate(df = "2010-2013") %>% 
  spread(df, count)
     
# winnowing histology

cancers <- cancers %>% group_by(which_cancer) %>% nest()

cancers[["data"]] <- 
  map2(cancers[["data"]], cancers[["which_cancer"]], df_top_n)

cancers <- unnest(cancers)

# treatment evidence =============================

tx <- 
  map2_dfc(c("nsrg_cnt", "rad_cnt", "stereo_cnt"), 
           c("Neurosurgery", "Radiation", "SRS"), 
           function(vctr, nm) {
             vctr <- 
               ifelse(is.na(cancers[[vctr]]), paste("No", nm), nm)
             tr <- data.frame(vctr, stringsAsFactors = FALSE)
             names(tr) <- nm
             tr
           })
          
cancers <- 
  bind_cols(cancers, tx) %>% 
  filter(which_cancer == "lung") %>% 
  mutate(bevacizumab = ifelse(mab == "Bevacizumab", "Bevacizumab", "None")) %>% 
  filter(!grepl("Autopsy", radiatn_v))

#annum counts=========================================

paper_products[["bev"]] <- 
  gp_count(cancers, 
           which_cancer, 
           seer_br_mets, 
           bevacizumab, 
           histo)

paper_products[["bev_race_sex"]] <- 
  gp_count(cancers, 
           which_cancer, 
           seer_br_mets, 
           race_v, 
           s_sex_v, 
           bevacizumab, 
           histo)

paper_products[["cccc_retreat"]] <- 
  gp_count(cancers, 
           which_cancer, 
           seer_br_mets, 
           SRS,
           s_sex_v, 
           bevacizumab
           )

paper_products[["bev_age_race_sex"]] <- 
  gp_count(cancers, 
           which_cancer, 
           seer_br_mets, 
           age_cut,
           race_v,
           s_sex_v, 
           bevacizumab
           )

paper_products[["with_rad"]] <- 
  gp_count(cancers, 
           which_cancer, 
           seer_br_mets, 
           Radiation,
           s_sex_v, 
           bevacizumab
           )


paper_products[["with_nsrg"]] <- 
  gp_count(cancers, 
           which_cancer, 
           seer_br_mets, 
           age_cut,
           Neurosurgery,
           s_sex_v, 
           bevacizumab
           )

#table ones===============================================

#to add: 
# income, education, charlson, 
# Consider filtering for first primary...?

varnames <- 
  c("age_dx", 
    "race",
    "histo", 
    "s_sex_v",
    "beh03v",
    "d_ssg00",
    "csmetsdxb_pub_v", 
    "csmetsdxliv_pub_v",   
    "csmetsdxlung_pub_v",
    "Stereo_RS", 
    "NSRG",
    "RAD", 
    "Bevacizumab"
)

to_numeric <- c("age_dx", "cs_size", "eod10_pn", "cs_mets")
cancers[,to_numeric] <- lapply(cancers[,to_numeric], as.numeric)

the_strata <- c("seer_br_mets", "default", "bevacizumab")

cancers <- reg_df(cancers)

paper_products[["table_ones"]] <- 
  map(the_strata, ~ tbl_one(cancers, varnames, .x))

names(paper_products[["table_ones"]]) <- the_strata

paper_products[["table_ones"]][["seer_br_mets"]]  <- 
  paper_products[["table_ones"]][["seer_br_mets"]][,1:2] 

paper_products[["table_ones"]][["bevacizumab"]]  <- 
  paper_products[["table_ones"]][["bevacizumab"]][,1:2]

paper_products[["table_ones"]][["mabs_br"]] <- 
  tbl_one(cancers, varnames, c("seer_br_mets", "bevacizumab"))[,1:4]

paper_products[["histo_key"]] <-
  with(cancers, table(hist03v, hist03v_v)) %>%
  data.frame %>%
  filter(Freq != 0 & hist03v_v %in% toupper(unique(cancers$histo))) %>%
  arrange(desc(Freq))

paper_products <-
  modify_if(paper_products, is.data.frame,
            ~ as.data.frame(ungroup(.x)), stringsAsFactors = FALSE)

write_rds(paper_products, "table-ones.rds")



