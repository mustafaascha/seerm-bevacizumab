
devtools::load_all("mabsr")
load_libs()

prds <- list()

prds[["drugs"]] <- 
  c(
    #"abatacept",
    "carboplatin",
    "cisplatin",
    #"denosumab",
    "dexamethasone",
    #"docetaxel",
    #"gemcitabine",
    "paclitaxel",
    "pemetrexed"
    #"vinorelbine"
  )

prds[["candidates"]] <- 
  c("age_cut",           #"cs_size", 
    "bone_mets",         "liver_mets",  "brain_mets", 
    "d_ssg00",
    "urbrur",
    "race",            
    "histo", "s_sex_v",     #"grade_v"
    #"NSRG",   
    #"RAD",    
    "surgery",
    "Stereo_RS", 
    paste("n", prds[["drugs"]], sep = "_")
  )

prds[["not_cox"]] <- 
  paste(c("bone_mets","brain_mets", "liver_mets", "surgery", "cs_size"), 
        collapse = "|"
        )
prds[["cox"]] <- 
  c("bmab01", 
    prds[["candidates"]][-grep(prds[["not_cox"]], prds[["candidates"]])])

prds[["impute"]] <- 
  c("seer_br_mets", 
    prds[["cox"]][-1], 
    "srv_time_mon", 
    "srv_time_mon_flag_v", 
    "stat_rec_v")

prds[["plts"]] <- 
  Filter(function(z) !grepl("\\*|NSRG|RAD|Stereo", z), 
         c("bmab01", "seer_br_mets", prds[["cox"]]))

prds[["to"]] <- 
  c("age_dx", 
    "race",
    "histo",
    "histo_to",
    "s_sex_v", 
    "beh03v",
    "d_ajcc_s_v",
    "csmetsdxb_pub_v",
    "csmetsdxliv_pub_v",
    "csmetsdxlung_pub_v",
    "Stereo_RS", 
    "NSRG", 
    "RAD",
    "Bevacizumab",
    paste("n_", prds$drugs, sep = ""))

prds[["tn"]] <- c("age_dx", "cs_size", "eod10_pn", "cs_mets")

prds[["strata"]] <- 
  list(
    br = "seer_br_mets", 
    overall = "default", 
    bev = "mab", 
    br_bev = c("seer_br_mets", "mab")
  )

prop <- list()
prop[["frm"]] <- reformulate(paste(prds[["cox"]][-1]),"bmab01")
prop[["on"]] <- list("bmab01", "bmab01*seer_br_mets", "seer_br_mets")



setwd("mabsr")
use_data(prds, prop, overwrite = T)
setwd("..")

