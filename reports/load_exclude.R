
library(tidyverse)
library(zeallot)

#options(warn = 2)
#options(error = browser)

paper_products <- list()

devtools::load_all("augur")
devtools::load_all("frequencies")
devtools::load_all("mabsr")

qr <- quietly(read_csv)

cancers <- 
  qr("cache/cancers.csv.gz", progress = FALSE)[["result"]] %>% 
  filter(dx_year >= 2010)

cancers <- 
  filter(cancers,
         siterwho == 22030 & 
         !(hist03v %in% mabsr::histo_exclusion))

to_numeric <- c("age_dx", "cs_size", "eod10_pn", "cs_mets")                  
cancers[,to_numeric] <- lapply(cancers[,to_numeric], as.numeric)  
rm(to_numeric)

cgp <- cancers %>% group_by(which_cancer, dx_year) %>% tidyr::nest()
cgp$data <- map2(cgp$data, cgp$which_cancer, df_top_n)
cancers <- unnest(cgp)
rm(cgp) 

cancers$cs_size[as.numeric(cancers$cs_size) >= 988] <- NA

cancers$age_cut <- 
  cut(as.numeric(cancers$age_dx), c(65, 70, 75, 80, 85, 115), include.lowest = TRUE,
    labels = c("65 to 69", "70 to 74", "75 to 79", "80 to 84", "85+"))


cancers$histo <- stringr::str_to_title(tolower(cancers$histo))

cancers$the_strata <- cancers$histo
  #with(cancers, ifelse(which_cancer == "skin", breslow, histo))


cancers$the_strata <- stringr::str_to_title(cancers$the_strata)
cancers$the_strata <- 
  ifelse(cancers$the_strata == "Unknown", "Other", 
         cancers$the_strata)

old_values <- 
  c("Unknown", 
    "Mucinous Adenocarcinoma",
    "Adenoid Cystic & Cribriform Ca.",
    "Lobular And Other Ductal Ca.", 
    "Not 2010\\+ Breast" 
    #"Face", "Upper Limb", "Lower Limb"
    )
new_values <- 
  c("Other", 
    "Other", 
    "Other", 
    "Duct Carcinoma",
    "Other"
    #"Face or Limb", "Face or Limb", "Face or Limb"
   )

cancers$the_strata <- 
  reduce2(old_values, new_values, 
          function(init, a, b) gsub(a, b, init),
          .init = cancers$the_strata)

cancers$race_v <- cancers$rac_recy_v

hispanic_origin_values <- 
  c("Cuban", 
    "Dominican Republic", 
    "Mexican", 
    "NHIA Surname Match Only", 
    "Other specified Spanish/Hispanic Origin including Europe",
    "Puerto Rican", 
    "South or Central American excluding Brazil", 
    "Spanish/Hispanic/Latino, NOS")

cancers[["race_v"]] <- 
  ifelse(cancers[["nhiade_v"]] %in% hispanic_origin_values & 
         cancers[["race_v"]] == "White", 
         "White Hispanic", 
         gsub("White", "White Non-Hispanic", cancers[["race_v"]]))
         

cancers <- 
  reduce2(list(c("American Indian", "Asian/Pacific Islander"), 
              c("American Indian", "Black", "Asian/Pacific Islander")), 
         c("breast", "skin"), 
         function(df, oths, cncr){
          df[["race_v"]] <- 
            ifelse(df[["which_cancer"]] == cncr & df[["race_v"]] %in% oths, 
                   "Other", df[["race_v"]])
            df},
          .init = cancers)


write_csv(cancers, "cache/cancers-analytic.csv.gz")
