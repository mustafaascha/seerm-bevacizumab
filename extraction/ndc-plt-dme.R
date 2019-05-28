
library(tidyverse)

devtools::load_all("augur")

codes_site <- "https://crn.cancer.gov/resources/ctcodes-drugs.csv.gz"

platin_agents <- 
  read_csv(codes_site) %>% 
  filter(grepl("latin", tolower(GENERIC_NAME)))

write_csv(platin_agents, 
          paste("cache/platinated-codes", 
                Sys.Date(), 
                ".csv.gz", 
                sep = ""))

plt_codes <- platin_agents[["NDC2"]]

dme_ndc_plt <- 
  extract_slice_codes("seerm", 
                      "dme", 
                      augur::infiles[["dme"]], 
                      plt_codes, 
                      "prod_srvc_id"
  ) %>% 
  bind_rows()

write_csv(dme_ndc_plt, "cache/mabs/dme-ndc-plt.csv.gz")







