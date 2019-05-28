

histo_key <- function(cncr_df) {
  with(cncr_df, data.frame(table(hist03v, hist03v_v))) %>% 
  filter(Freq != 0) %>% # & hist03v %in% toupper(unique(cncr_df[["histo"]]))) %>%
  arrange(desc(Freq))
}



