
file_sig <- function(filetype) {
  reduce2(
    c(
"bev-patients-", "saf|pa[tr]", "cache\\/", "\\..*"
), 
    c("", "", "", ""), 
    function(e, f, g) gsub(f, g, e),
    .init = filetype
  )
}

read_bev <- function(fp, ...) {
  message(paste("Reading one of the", file_sig(fp), "files"))
  bev_df <- 
    quietly(read_csv)(fp, col_types = choose_col(file_sig(fp)), ...)[["result"]]
  structure(bev_df, class = c(file_sig(fp), class(bev_df)))
} 
 
read <- function(filetype, ...) {
  map_dfr(
    list.files("cache", pattern = filetype, full.names = T), 
    function(fp) read_bev(fp, ...)
  )
}

from_thru <- function(a_df) {
  a_df %>% 
  mutate(from = paste(from_dtd, from_dtm, from_dty, sep = "-")) %>% 
  mutate(thru = paste(thru_dtd, thru_dtm, thru_dty, sep = "-"))
}

choose_col <- function(fp) {
  switch(fp, 
    dme = dme_bev_cols(), 
    out = out_bev_cols(), 
    med = med_bev_cols(), 
    pde = pde_bev_cols(), 
    nch = nch_bev_cols()
  )
}

hcpcs_lookup <- function() {
    rename(mabsr::lookups$chemo_therapy_procs, hcpcs = PX2) %>% 
    mutate(hcpcs_nm = paste(DESCRIPTION, 
                            CHEMOCAT, 
                            CHEMO_TYPE, 
                            TREATMENT_TYPE, 
                            CHEMOTHERAPY,
                            sep = "***")
    ) %>%  
    select(hcpcs, hcpcs_nm) %>% 
    bind_rows(select(mabsr::ndc_hcpcs_chemo, hcpcs, hcpcs_nm))
}

inner_join_hcpcs <- function(z) {
  if (("PX2" %in% names(z))) z[["hcpcs"]] <- z[["PX2"]]
  inner_join(z, hcpcs_lookup())
}

left_join_hcpcs <- function(z) {
  if (("PX2" %in% names(z))) z[["hcpcs"]] <- z[["PX2"]]
  left_join(z, hcpcs_lookup())
}

inner_join_ndc <- function(z) {
  inner_join(z, mabsr::ndc_hcpcs_chemo)
}

left_join_ndc <- function(z) {
  left_join(z, mabsr::ndc_hcpcs_chemo)
}

select_rename <- function(a_df) {
  if (!("ndc_nm" %in% names(a_df))) a_df[["ndc_nm"]] <- NA
  select(a_df, 
         patient_id, from, thru, hcpcs = hcpcs_nm, ndc = ndc_nm)
}

join <- function(.df, ...) UseMethod("join", .df)

join.dme <- function(.df, ...) {
  #from_thru(.df) %>% left_join_hcpcs() %>% left_join_ndc()
  compose(select_rename, from_thru, left_join_ndc, left_join_hcpcs)(.df)
}

join.out <- function(.df, ...) {
  compose(select_rename, from_thru, inner_join_hcpcs)(.df)
}

join.med <- function(.df, ...) {
  inner_join_ndc(.df)
}

join.pde <- function(.df, ...) {
  .df %>% 
  mutate(from = paste(srvc_day, srvc_mon, srvc_yr, sep = "-"), 
         thru = NA, 
         hcpcs = NA
         ) %>% 
  select(patient_id, from, thru, hcpcs, gnn)
}

join.default <- function(.df, ...) {
  compose(select_rename, from_thru, inner_join_hcpcs)(.df)
}


mk_other_chemo <- function(z) {
  censor <- 
    group_by(z, drug) %>% 
    summarise(count_ = n()) %>% 
    filter(count_ < 10) %>% 
    select(drug) %>% 
    unlist() %>% 
    unname()
  z[z[["drug"]] %in% censor,"drug"] <- "Other"
  z
}

munge_merged <- function(chemo_df, cancers_df) {
 left_join(chemo_df, cancers_df) %>% 
    #filter(seq_num %in% c("00", "01")) %>% 
    modify_at(c("from", "thru"), lubridate::dmy) %>% 
    mutate(dx_date = lubridate::dmy(paste("01", dx_month, dx_year, sep = "-"))
    ) %>% 
    mutate(from_since_dx = as.numeric(difftime(from, dx_date, units = "days")),
           thru_since_dx = as.numeric(difftime(thru, dx_date, units = "days"))
    ) %>% 
    mutate(
      RAD = ifelse(rad_cnt >= 1, "RAD", "No RAD"),
      NSRG = ifelse(nsrg_cnt >= 1, "NSRG", "No NSRG"),
      Stereo_RS = ifelse(stereo_cnt >= 1, "SRS", "No SRS")
    ) %>% 
    select(
      -from, 
      -thru, 
      -dx_date, 
      -dx_year, 
      -dx_month, 
      -srv_time_mon, 
      -srv_time_mon_flag_v, 
      -csmetsdxb_pub_v, 
      -csmetsdxliv_pub_v, 
      -stat_rec_v, 
      -nsrg_cnt, 
      -rad_cnt,
      -stereo_cnt
      #, -Bevacizumab
    ) %>% 
    rename(
      urban = urbrur,
      stage = d_ssg00,
      m_rad = RAD, 
      m_nsrg = NSRG, 
      m_srs = Stereo_RS,
      no_surg = no_surg_v, 
      m_brain_mets = counts_dx_matches,
      bmab = mab, 
      age = age_cut, 
      sex = s_sex_v
    ) %>% 
    #filter(from_since_dx >= 0) %>% 
    filter(!grepl("^SMALL", hist03v_v)) %>% 
    select(-hist03v_v) %>% 
    modify_at("patient_id", openssl::md5)
}


fix_term_names <- function(the_df) {
  modify_at(
    the_df, 
    "term",
    ~ reduce2(
      c("age_cut",  "bone_mets", "liver_mets", "d_ssg00",
        "urbrur", "s_sex_v", "NSRGNSRG", "RADRAD", "surgerySurgery",
        "Stereo_RSSRS", "raceB_", "raceC_Other", "histoOther",
        "brain_metsYes"
       ),
      c("age: ", "", "", "Stage: ",
        "", "", "Neurosurgery", "Radiotherapy", "Any Onco-Surgery",
        "Stereo RS", "Race: ", "Race: Other", "Non-Adenocarcinoma",
        "Brain Mets"
       ),
      function(a, b, d) gsub(b, d, a),
      .init = .x
      )
  )
}


fix_to_names <- function(the_table) {
  rownames(the_table) <- 
    reduce2(
    c("age_dx", "race", "[ABC]_", "d_ssg00", "csmetsdxb_pub_v", "csmetsdxliv_pub_v",
      "csmetsdxlung_pub_v",
      "beh03v", "s_sex_v = ", "histo = Other", "Stereo_RS", "NSRG = NSRG", "RAD = RAD"), 
    c("Age at Dx", "Race", "", "Stage", "Bone Mets", "Liver Mets",
      "Lung Mets",
      "Behavior", "", "Non-Adenocarcinoma", "Stereo RS", "Neurosurgery", "Radiotherapy"), 
    function(a, b, d) gsub(b, d, a),
    .init = rownames(the_table)
    )

  colnames(the_table) <- 
    reduce2(
      c("SEER_Positive:Bevacizumab", "SEER_Negative:Bevacizumab"), 
      c("SEER Positive", "SEER Negative"), 
      function(a, b, d) gsub(b, d, a),
      .init = colnames(the_table)
    )
  the_table

}




table_bev_odds <- function(r) {
  imap_dfr(r, 
           function(dta, nm) {
             dta %>% 
               modify_at(
                 c("estimate", "conf.low", "conf.high"), 
                 ~ sprintf("%.2f", as.numeric(.x))
               ) %>% 
               mutate(
                 model = ifelse(nm == "mv_glm", "Multivariable", "Univariable"),
                 OR = paste(estimate, " (", conf.low, "-", conf.high, ", p: ", p.value, ")", sep = "")
               )
           }) %>% 
    select(-estimate, -conf.low, -conf.high, -p.value) %>% 
    #fix_term_names() %>% 
    spread(model, OR)
  
}

glm_ggplot <- function(z) {
  ggplot(z, aes(x = reorder(term, estimate), y = estimate, color = mdl)) +
    geom_hline(yintercept = 1, color = "grey60") +
    geom_point(position = position_dodge(1)) +
    geom_errorbar(aes(ymin = conf.low, ymax = conf.high), position = position_dodge()) +
    coord_flip() +
    facet_wrap("pop") +
    labs(x = "Odds Ratio: Bevacizumab Treatment (95% CI)",
         y = "",
         color = "Model")  +
    scale_color_hue(labels = c("Multi.", "Uni.")) +
    theme_bw() +
    theme(axis.text.x = element_text(color = "black"),
          axis.text.y = element_text(color = "black"))
}




read_chemo <- function(fp) {
  # ten most frequent drugs used with Bev
  keep_these <- 
    c(
      "dexamethasone",
      "carboplatin",
      "pemetrexed", 
      "paclitaxel", 
      "docetaxel",  
      "abatacept",  
      "gemcitabine",
      "denosumab",  
      "cisplatin",  
      "vinorelbine"
      )

  read_csv(fp) %>% 
    select(-is_chemo) %>%
    distinct() %>%
    modify_at(
      "drug",
      ~ forcats::fct_recode(.x,
        testosterone = "testosterone cypionate and estradiol cypion",
        interferon_alpha1 = "interferon alfacon-1",
        interferon_alpha2 = "interferon, alfa-2a"
      ) %>%
      gsub("[^a-zA-Z0-9]+", "_", .) %>%
      substr(1, 40) %>% 
      gsub("_$", "", .) %>% 
      (function(a) ifelse(grepl("estosteron", a), "testosterone", a))
      ) %>% 
    filter(drug %in% keep_these) %>% 
    modify_at(c("from", "thru"), lubridate::dmy) %>%
    mutate(drug_dur = as.numeric(thru - from)) %>% 
    group_by(patient_id, drug) %>% 
    summarise(
      start = as.character(min(from, na.rm = T)),
      d = sum(drug_dur, na.rm = T),
      n = n()
    ) %>% 
    ungroup() %>% 
    modify_at(c("d", "n"), ~ ifelse(.x <= 0, 1, .x)) %>% 
    #(function(a) browser()) %>% 
    gather(k, v, -patient_id, -drug) %>% 
    mutate(d_nm = paste(k, drug, sep = "_")) %>% 
    select(-drug, -k) %>% 
    spread(d_nm, v) %>% 
    modify_at(
      paste("d", mabsr::prds$drugs, sep = "_"), 
      ~ case_when(is.na(.x) ~ "None", .x < 30 ~ "1-30 days", .x > 30 ~ "> 30 days")
    ) %>% 
    modify_at(
      paste("n", mabsr::prds$drugs, sep = "_"), 
      ~ case_when(is.na(.x) ~ "None", .x >= 1 ~ "Prescribed")
    )
}



