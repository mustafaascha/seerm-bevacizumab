

library(tidyverse); library(zeallot); devtools::load_all("augur")
cancers <- read_csv("cache/cancers_joined-vars.csv.gz", progress = FALSE)

#joining for two things: code counts and date nearest to primary diagnosis

dx_matches <- claims_dates_df("cache/diagnoses")
img_df <- claims_dates_df("cache/dx-imaging")

date_differences <- days_between_claims(dx_matches, img_df, "cpt_img")

cancers <- left_join(cancers, date_differences)

rad_df <- claims_dates_df("cache/rad")
date_differences <- days_between_claims(dx_matches, rad_df, "dx_rad")

cancers <- left_join(cancers, date_differences)

#=============================================

write_csv(cancers, "cache/cancers_prerecode.csv.gz")


