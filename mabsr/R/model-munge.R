
too_few_in_cats <- function(prd) {
  if (is.character(prd)) return(any(table(prd) <= 5))
  else return(FALSE)
}

#' @import dplyr
reg_df <- function(a_df, keep_stage = FALSE) {
  to_return <- 
    a_df %>% 
    modify_at("no_surg_v", ~ ifelse(is.na(.x), "Unknown", .x)) %>% 
    mutate( 
           Bevacizumab = 
             ifelse(mab == "Bevacizumab", "Bevacizumab", "None"), 
           Stereo_RS = ifelse(stereo_cnt > 2, "SRS", "None"), 
           NSRG = ifelse(nsrg_cnt > 2, "NSRG", "None"), 
           RAD  = ifelse(rad_cnt > 2,  "RAD", "None"),
           bone_mets = ifelse(csmetsdxb_pub_v == "Yes", "Bone Mets", 
                          ifelse(csmetsdxb_pub_v == "Unknown", NA, "Absent")), 
           liver_mets = ifelse(csmetsdxliv_pub_v == "Yes", "Liver Mets", 
                          ifelse(csmetsdxliv_pub_v == "Unknown", NA, "Absent")), 
           brain_mets = ifelse(seer_br_mets == "SEER_Positive", "Yes", 
                          ifelse(seer_br_mets == "SEER_Negative", "Absent", NA)),
           race = ifelse(race_v %in% c("White Non-Hispanic"), "White Non-Hispanic" , "Other"), 
           histo = ifelse(hist03v_v != "ADENOCARCINOMA, NOS", "aaother", hist03v_v), 
           surgery = ifelse(no_surg_v == "Surgery performed", "Surgery", "Not documented")
           ) %>% 
    mutate(bmab01 = as.numeric(Bevacizumab == "Bevacizumab")) %>% 
    mutate(histo_to = fct_lump(hist03v_v, prop = 0.1)) %>% 
    filter(which_cancer == "lung") %>% 
    modify_at("urbrur", ~ fct_collapse(.x, 
                            A_Not_Metro = c("Urban", "Less Urban", "Rural"),
                            B_Metro = c("Big Metro", "Metro"))
    ) %>% 
    modify_at("histo", ~ as.character(fct_recode(.x, adeno = "ADENOCARCINOMA, NOS"))
    ) %>% 
    modify_at("age_cut", function(ac) {
      ifelse(grepl("8|75", as.character(ac)), "75+", as.character(ac))
    }) %>% 
    modify_at("age_cut", ~ factor(.x, levels = c("65 to 69", "70 to 74", "75+"))) %>% 
    modify_at("race", function(rv) {
      fct_recode(rv, A_WNH = "White Non-Hispanic", B_Other = "Other")
    }) %>% 
    modify_at(c("age_cut", "urbrur", "srv_time_mon_flag", "race"), as.character) %>% 
    modify_at(c("NSRG", "RAD", "Stereo_RS"), ~ ifelse(is.na(.x), "None", .x)) %>% 
    modify_at(mabsr::prds$tn, as.numeric) %>%
    modify_at("d_ssg00", ~ gsub("[AB]", "", .x)) %>% 
    filter(d_ssg00 != "Not applicable") %>% 
    modify_at(
      "d_ssg00", 
      function(x) { 
        if (!keep_stage) ifelse(x == "Stage I", "0I", gsub("Stage", "", x))
        else x
    })
    #modify_at("d_ssg00", ~ relevel(factor(.x), ref = "Stage I"))
  #histos_to_keep <- table(to_return$hist03v_v) > 20
  #to_return[["histo"]] <- 
  #  ifelse(to_return$hist03v_v %in% names(histos_to_keep)[histos_to_keep], 
  #         to_return$hist03v_v, "Other")
  to_return
}

#' @import dplyr
#' @import survival
srv_df <- function(a_df, mb = FALSE) {
  if (mb) {
    mb_fn <- function(z) mutate(z,
      srv_time_mab = make_time(srv_time_mon, mab_start, dx_year, dx_month)
    )
  } else mb_fn <- function(a) a
  
  nonzero_surv <- 
    "Complete dates are available and there are more than 0 days of survival"
  nonzero_cens <- 
    "Incomplete dates are available and there cannot be zero days of follow-up"
  a_df <- 
    filter(a_df,  
           srv_time_mon_flag_v %in% c(nonzero_surv, nonzero_cens)) %>% 
    mutate(srv_time = as.numeric(gsub("^0+", "", srv_time_mon)),
           alive = stat_rec_v, 
           bmab10 = 1 - bmab01) %>% 
    mb_fn() %>% 
    modify_at(paste("d", mabsr::prds$drugs, sep = "_"), 
              ~ relevel(factor(.x), ref = "None"))
  a_df
}

make_time <- function(s, m, dy, dm) {
  dt <- ymd(paste(dy, "-", dm, "-15", sep = ""))
  tt_m <- as.numeric(difftime(dt, ymd(m), units = "days")) / 30
  tt_m <- ifelse(is.na(tt_m), 0, tt_m)
  as.numeric(gsub("^0+", "", s)) - tt_m
}

sw_mab <- function(a_df) {
  mutate(
    a_df, 
    bmab01 = case_when(
      bmab01 == 0 ~ "a_none", 
      bmab01 == 1 ~ "bevacizumab",
      bmab01 == "a_none" ~ "0", 
      bmab01 == "bevacizumab" ~ "1"
    )
  )
}





