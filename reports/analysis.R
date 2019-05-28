devtools::load_all("mabsr")
load_libs()
switches(
  do_reload       = FALSE
  , do_glms       = F, do_coxs       = T
  , do_match      = F, do_mtch_mdls  = F
  , srs_analysis  = F, do_tbl_ones   = F
  , write_results = F
  , .force = T
  )
set.seed(1001)

if (do_reload | !exists("cancers")) { source("reports/assemble-cancers.R") 
} else { cancers <- read_csv("cache/cancers-bev.csv.gz") }

if (!exists("results")) results <- list()
if (!exists("unmatched_data")) { source("reports/unmatched-data.R") } 

if (do_glms)  { source("reports/glms.R") }
if (do_match) { source("reports/match.R") } 
if (do_coxs)  { 
  source("reports/cox-models.R") 
  source("reports/figures.R")
}

drug_prds <- c(paste("n_", prds$drugs,sep = ""), paste("bmab01*n_", prds$drugs, sep = ""))

if (do_mtch_mdls) { 
  source("reports/prop-model-params.R")
  results[["propensity_hrs"]]    <- map(prop_mdl_params,    ~ do_prop_mdls(.x))
  results[["bm_propensity_hrs"]] <- map(bm_prop_mdl_params, ~ do_prop_mdls(.x))  
} 

if (srs_analysis){ source("reports/srs.R") }
if (do_tbl_ones) { source("reports/table-ones.R") } 

results <- modify_if(results, is.data.frame, compose(as.data.frame, ungroup))

# save environment? 
write_rds(results$to, "results/table-ones.rds")
save(results, file = "results/models-mabs.rdata")
write_rds(results, "results/models-mabs.rds", compress = "gz")
writeLines("", "cache/analysis")

