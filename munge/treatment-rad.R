



library(tidyverse)
library(zeallot)
devtools::load_all("augur") 
cancers <- read_csv("cache/cancers_nsrg.csv.gz")
#joining for two things: code counts and date nearest to primary diagnosis

cancers <- add_claims(cancers, read_matches("cache/rad"), "rad")

#for (i in seq_along(claims_dfs)) {
#  cancers <- add_claims(cancers, claims_dfs[[i]], names(claims_dfs)[i])
#}

#cancers <- reduce2(claims_dfs, names(claims_dfs), add_claims, .init = cancers)

write_csv(cancers, "cache/cancers_rad.csv.gz")



