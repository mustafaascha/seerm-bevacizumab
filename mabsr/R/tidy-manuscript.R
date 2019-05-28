reduce_terms <- function(the_terms) {
  reduce2(
    c("n_|prescribed", 
      "seer_br_metsseer_positive",
      "age_cut", 
      "d_ssg00", 
      "urbrurnot_metro",
      "race[bc]_", 
      "stereo_rssrs",
      "nsrgnsrg", 
      "histoother", 
      "s_sex_vmale",
      "radrad", 
      "bmab01",
      "liver_metsliver mets",
      "bone_metsbone mets",
      "surgerysurgery",
      "Summary "
    ), 
    c("",
      "Had brain mets",
      "Age: ", 
      "Summary Stage: ", 
      "Not Metro Area",
      "Race: ", 
      "SRS", 
      "Neurosurgery",
      "Non-Adenocarcinoma",
      "Male", 
      "Radiotherapy", 
      "Bevacizumab",
      "Liver Mets",
      "Bone Mets",
      "Surgery", 
      ""
    ), 
    function(a, b, d) gsub(b, d, a),
    .init = tolower(the_terms)
  )
}


clean_terms <- function(mdl_df) {
  modify_at(mdl_df, "term", compose(str_to_title, reduce_terms))
}

munge_mdl_df <- function(mdl_df) {
 mdl_df %>% 
    modify_at(c("estimate", "conf.low", "conf.high"), exp) %>% 
    modify_at(c("estimate", "conf.low", "conf.high"), ~ sprintf("%.2f", .x)) %>% 
    #modify_at("p.value", ~ format(.x, digits = 4, scientific = F)) %>% 
    modify_at("p.value", ~ sprintf("%.4f", .x)) %>% #~ format(.x, digits = 4, scientific = F)) %>% 
    clean_terms() %>% 
    select(term, estimate, conf.low, conf.high, p.value) 
    
}

indiv_plots <- function(p_df) {
  p <- function(b) b %>% 
    modify_at("term", ~ gsub("bevacizumab:", "", .x)) %>%
    modify_at("matched", str_to_title) %>% 
    ggplot(aes(y = hr, x = reorder(term, hr), lty = matched)) + 
    geom_hline(yintercept = 1, color = "grey60", lty = 2) + 
    geom_point(position = position_dodge(0.75)) + 
    geom_errorbar(aes(ymin = lci, ymax = hci), 
                  position = position_dodge(0.75)) + 
    #facet_wrap("pop", nrow = 2, scales = "free") + 
    coord_flip() + 
    theme(axis.text.y = element_text(color = "black"), 
          legend.position = c(0.75, 0.15), 
          legend.background = element_rect(fill = alpha('black', 0.1))) + #"bottom") + 
    labs(y = "HR (95% CI)", x = "Predictor", lty = "Matching")
  p_df %>%  
    group_by(pop, int, weighted) %>% 
    nest() %>% 
    (function(a) {
      set_names(map(a$data, p), 
                gsub(" ", "_", paste(a$pop, a$int, a$weighted, sep = "_")))
    })
}

unnest_cox <- function(cox_models) {
  cox_models %>% map_dfr(2)
}

rm_ <- function(a, rm) gsub(rm, "", a)

separate_estimate <- function(tidy_mdl_df) {
  tidy_mdl_df %>% 
  mutate(hr  = gsub("\\ \\(.*", "", HR),
         lci = gsub(".*\\(|-.*", "", HR),
         hci = gsub(".*-|,.*", "", HR),
         phr = gsub(".*p:\\ |\\)", "", HR)) %>% 
  select(-HR) %>% 
  modify_at(c("hr", "lci", "hci", "phr"), as.numeric)
}

int_plot <- function(p_df) {
  p_df %>% 
    ggplot(aes(y = hr, x = reorder(term, hr), color = matched)) + 
    geom_hline(yintercept = 1, color = "grey60", lty = 2) + 
    geom_point(position = position_dodge(0.75)) + 
    geom_errorbar(aes(ymin = lci, ymax = hci), 
                  position = position_dodge(0.75)) + 
    #facet_wrap("pop", nrow = 2, scales = "free") + 
    coord_flip() + 
    theme(axis.text.y = element_text(color = "black"))
}

plot_glms <- function(results) {
  with(results, 
       list(
         overall_mv = mv_glm, 
         overall_uni = uni_glm, 
         bm_mv = bm_mv_glm, 
         bm_uni = bm_uni_glm)
  ) %>% 
    imap_dfr(~ mutate(.x, mdl = .y)) %>% 
    separate("mdl", c("pop", "mdl"), "_", fill = "left") %>% 
    modify_at("pop", ~ ifelse(.x == "bm", "Brain Mets", "Overall")) %>% 
    clean_terms() %>% 
    modify_at(c("estimate", "conf.low", "conf.high", "p.value"), as.numeric) %>% 
    filter(!grepl("Situ", term))
}

plot_coxs <- function(results, f = T) {
  results[c("mv_cox", "bm_mv_cox")] <- 
    map(results[c("mv_cox", "bm_mv_cox")], 
        ~ rename(.x[[2]], hr = estimate, lci = conf.low, hci = conf.high, phr = p.value) %>% 
          mutate(matched = "unmatched", weighted = "unweighted") %>% 
          modify_at(c("hr", "lci", "hci", "phr"), as.numeric))
  if (f) fn <- filter_cox else fn <- function(z) z
    
  with(results, 
       list(overall_uni = compose(separate_estimate, unnest_cox)(propensity_hrs),
            overall_mv  = mv_cox,
            bm_uni      = compose(separate_estimate, unnest_cox)(bm_propensity_hrs),
            bm_mv       = bm_mv_cox
       )
  ) %>% 
    imap_dfr(~ mutate(.x, mdl = .y)) %>% 
    separate("mdl", c("pop", "mdl"), "_", fill = "left") %>% 
    clean_terms() %>% 
    modify_at("pop", ~ ifelse(.x == "bm", "Brain Mets", "Overall")) %>% 
    fn()
}

filter_cox <- function(a_df) { 
  a_df %>% 
  filter(weighted != "weighted" & 
           matched != "matched" &
           !(pop == "overall" & term == "Bevacizumab") &  
           !(grepl("Bevacizumab:brain Mets", term)) & 
           term != "Brain Mets" & 
           !(grepl("Situ", term))
           #!(matched == "matched" & term != "Bevacizumab") & 
         #  !grepl(":", term)
         ) %>% 
    bind_rows(filter(a_df, pop == "overall" & weighted == "weighted" & term == "Bevacizumab"))
}

clean_table_one <- function(to) {
  to <- to[,-c(5:6)]
  to <- to[!grepl("Malignant|Bevacizumab|NSRG", rownames(to)),]
  rownames(to) <- 
    reduce2(
      c(
        "age_dx"
        , "race = B_Other"
        , "s_sex_v = "
        , "d_ssg00"
        , "csmetsdxb_pub_v"
        , "csmetsdxliv_pub_v"
        , "csmetsdxlung_pub_v"
        , "^n_"
        , "carboplatin"
        , "cisplatin"
        , "dexamethasone"
        , "paclitaxel"
        , "pemetrexed"
        , "histo = Other"
        , " = Prescribed"
        ), 
      c(
        "Age"
        , "Non-White"
        , "" #
        , "Stage"
        , "Bone Metastases"
        , "Liver Metastases"
        , "Lung Metastases"
        , ""#
        , "Carboplatin"
        , "Cisplatin"
        , "Dexamethasone"
        , "Paclitaxel"
        , "Pemetrexed"
        , "Non-Adenocarcinoma"
        , "" #
      ), 
      function(a, b, d) gsub(b, d, a),
      .init = rownames(to)
    )
  colnames(to) <- 
    reduce2(
      c("SEER_Negative", 
        "SEER_Positive", 
        "Bevacizumab", 
        "None"), 
      c("No BM", 
        "BM", 
        "Bev.", 
        "No Bev."),
      function(a, b, d) gsub(b, d, a),
      .init = colnames(to)
    )
  to
}

arrange_srv <- function(a, b, d) {
  gridExtra::grid.arrange(
    a, b, d,
    layout_matrix = 
      rbind(
        c(1, 1, 1), 
        c(1, 1, 1),
        c(1, 1, 1),
        c(1, 1, 1),
        c(2, 2, 2),
        c(2, 2, 2),
        c(2, 2, 2),
        c(3, 3, 3)
      )
  )
}

annotation_gen <- function(the_tbl) {
  the_tbl <- the_tbl[the_tbl[["matched"]] != "matched",]
  
  function(plt, mdl, overall = F, hr_term = "HR") {
    if (overall) the_col <- the_tbl[["Overall"]] else the_col <- hr_table[["Brain_Mets"]]
    if (mdl == "mv") {
      lb <- "Multivariable" 
      vl <- the_col[the_tbl[["mdl"]] == "mv"]
    } else { 
      lb <- "Univariable"
      vl <- the_col[the_tbl[["mdl"]] == "uni"]
    }
    the_label <- paste(lb, hr_term, ":\n", vl)
    plt + annotate(geom = "label", label = the_label, x = 44, y = 0.8)
  }
}
