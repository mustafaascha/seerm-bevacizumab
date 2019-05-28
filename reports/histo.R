
library(tidyverse) 
library(zeallot)
library(Hmisc)

papes <- read_rds("paper_products.rds")

devtools::load_all("augur")
devtools::load_all("frequencies")

# For description of histology codes in methods =============================

histo_key <- 
  papes$histo_key %>% 
  rename(histo = hist03v_v) %>% 
  select(-Freq) %>% 
  modify_at("histo", 
            function(z) str_to_title(gsub("\\,\\ NOS$", "", z))) %>% 
  modify_at("histo", 
            function(z) tolower(gsub("Mal\\.\\ Mel\\.\\ In\\ Junct\\.\\ Nevus", 
                                     "malignant Melanoma in junctional nevus", z))) %>% 
  modify_at("histo", 
            function(z) gsub("ca\\.$", "carcinoma", z)) %>% 
  group_by(which_cancer, histo) %>% 
  nest()

#I don't need to do this biggest/smallest stuff...I need to remove it
histo_key[["codes"]] <- 
  map_chr(histo_key[["data"]], 
          function(df) {
            df[["hist03v"]] <- as.numeric(as.character(df[["hist03v"]]))
            df <- arrange(df, hist03v)
            biggest <- max(df[["hist03v"]], na.rm = TRUE)
            smallest <- min(df[["hist03v"]], na.rm = TRUE)
            if (biggest == (smallest + nrow(df) - 1) & nrow(df) > 1) {
              return(paste(smallest, biggest, sep = "-"))
            } else {
              rle_hyp(df[["hist03v"]])
            }
          })  

# Function for in-text reference
hst_c <- function(cnc, hst) {
  if (is.numeric(hst)) { 
    which_hst <- hst
  } else if (is.character(hst)) {
    which_hst <- grep(hst, histo_key$histo[histo_key$which_cancer == cnc])
  }
  hst <- histo_key$histo[histo_key$which_cancer == cnc][which_hst]
  with(histo_key, 
    list(h = hst, c = codes[which_cancer == cnc & histo == hst]))
}














