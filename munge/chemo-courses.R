
max_rows <- Inf

library(tidyverse)
library(foreach)
library(doMC)
devtools::load_all("mabsr")

# table out the code types to see what we can investigate
the_filetypes <- c("dme", "medpar", "nch", "outsaf", "pdesaf")

filepaths <- 
  map(paste("bev-patients-", the_filetypes, "*", sep = ""), 
      ~ list.files("cache", pattern = .x, full.names = T)) %>% 
  map(function(z) Filter(function(a) !grepl("bk|chemo", a), z)) %>% 
  set_names(the_filetypes)

#chemotherapy <- 
#  map_dfr(filepaths, 
#          function(fps)  map_dfr(fps, ~ join(read_bev(.x, n_max = 1e4))))

#registerDoMC(6)

chemotherapy <- 
  map(filepaths, 
        function(fps) {
          #mclapply(fps, 
          lapply(fps, 
            function(fp) {
              tryCatch(
                distinct(join(read_bev(fp, n_max = Inf))), 
                error = function(e) {
                  data.frame(patient_id = "", from = "", thru = "", hcpcs = "", ndc = "")
                }
              )
          })
        }) %>% 
  map(compose(distinct, bind_rows))

#chemotherapy <- list()
#
#for (i in seq_along(filepaths)) {
#  chemotherapy[[names(filepaths)[i]]] <- 
#    mclapply(filepaths[[i]], function(fp) distinct(join(read_bev(fp, n_max = Inf)))) %>% 
#    compose(distinct, bind_rows)()
#}

chemo <- map_if(chemotherapy, is_packed, unpack_chemo)

pde <- chemo[["pdesaf"]]
chemo[["pdesaf"]] <- NULL

chemo <- 
  bind_clean(chemo) %>% 
  left_join(mabsr::drug_names) %>% 
  select(-drug) %>% 
  rename(drug = drugs_short)

write_csv(chemo, "cache/chemotherapy.csv.gz")
write_csv(chemo, "cache/chemo-pdesaf.csv.gz")

#for (i in seq_along(types_and_joins)) {
#  chemotherapy[[names(t_n_j)[i]]] <- 
#    inner_join(read(names(t_n_j)[i]), t_n_j[[i]])
#}
#products <- 
#  read_tsv("documentation/products.txt", 
#           cols = cols_only(PRODUCTNDC = col_))
#package <- read_tsv("documentation/package.txt")


