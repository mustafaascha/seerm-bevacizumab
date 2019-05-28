

library(tidyverse)
library(devtools)

filetypes <- 
  c(
    dme    = "dme", 
    outsaf = "outsaf", 
    #medpar = "medpar",
    pdesaf = "pdesaf",
    nch    = "nch"
  )

lookups <- 
  map(list.files("../documentation", full.names = TRUE, pattern = "csv"), 
      read_csv)
names(lookups) <- gsub("\\.csv$", "", list.files("../documentation", pattern = "csv"))
names(lookups) <- gsub("-", "_", names(lookups))

ndc_hcpcs_xwalk <- 
  readxl::read_xlsx("../documentation/ndc-hcpcs-xwalk.xlsx", 
                    sheet = "09-05-2018 NDC-HCPCS XWalk")

ndc_hcpcs_xwalk <- 
  ndc_hcpcs_xwalk[,
    c("NDC", 
      "HCPCS", 
      "NDC Label",
      "HCPCS Description"
      )
  ]

names(ndc_hcpcs_xwalk) <- c("ndc", "hcpcs", "ndc_nm", "hcpcs_nm")

ndc_hcpcs_xwalk[["ndc"]] <- gsub("\\-", "", ndc_hcpcs_xwalk[["ndc"]])

ndc_hcpcs_chemo <- 
  ndc_hcpcs_xwalk[
    ndc_hcpcs_xwalk[["ndc"]] %in% lookups[["chemo_therapy_drugs"]][["NDC2"]]
  ,]


drug_names <- 
  read_csv("data-raw/drug-names.csv") %>% 
  modify_at("drugs_short", tolower)


# https://seer.cancer.gov/oncologytoolbox/canmed/ndconc/
canmed <- list()
canmed[["ndc"]] <- read_csv("data-raw/ndconc_results.csv")
canmed[["hcpcs"]] <- read_csv("data-raw/hcpcs_results.csv")


usethis::use_data(
  canmed, 
  ndc_hcpcs_xwalk, 
  ndc_hcpcs_chemo, 
  filetypes, 
  lookups, 
  drug_names,
  overwrite = TRUE
)


rm_these <- 
  c(
    "infusion therapy", 
    "photochemotherapy", 
    "planning chemotherapy", 
    "prescription oral chemotherapeutic", 
    "unclassified biologics", 
    "cns infusion", 
    "unlisted chemotherapy", 
    "Other", 
    "stage iii chemo", 
    "antineoplastic drug", 
    "cns infusion therapy", 
    "phototherapy", 
    "unlisted therapeutic, prophylactic or diagnostic intravenous or intra-arterial i", 
    "stage iii adjuvant", 
    "unlisted therapeutic, prophylactic, or diagnostic intravenous or intra-arterial", 
    "er/pr positive breast"
  )

usethis::use_data(rm_these, overwrite = TRUE)

