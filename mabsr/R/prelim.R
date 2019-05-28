

munge_cc <- function(the_df, rm_these) {
  the_df %>% 
  group_by_at(vars(-from_since_dx, -thru_since_dx)) %>% 
    summarise(from_since_dx = min(from_since_dx, na.rm = T), 
              thru_since_dx = max(thru_since_dx, na.rm = T)) %>% 
    ungroup() %>% 
    mutate(from_thru = paste(from_since_dx, thru_since_dx, sep = "**")) %>% 
    filter(!(drug %in% rm_these | is.na(drug)) & from_since_dx < (365 * 3)) %>% 
    mutate(alph = ifelse(bmab == "Bevacizumab", 1, 0.1))
}

censor <- function(the_df) {
  enough_observations <- 
    the_df %>% 
    group_by(drug, bmab) %>% 
    summarise(count_ = n()) %>% 
    filter(count_ > 10) %>% 
    select(drug) %>% 
    unlist() %>% 
    unname()
  the_df %>% 
    filter(drug %in% enough_observations)
}

drug_freq_order <- function(the_df) {
  the_df[["drug"]] <- 
    reduce2(c("[^0-9a-zA-Z]+", "_$"), 
            c("_", ""), 
            function(a, b, d) gsub(b, d, a), 
            .init = the_df[["drug"]]
    )
  the_df[["drug"]] <- str_to_title(the_df[["drug"]])
  orders <- sort(table(the_df[["drug"]]))
  the_df[["drug"]] <- factor(the_df[["drug"]], levels = rev(names(orders)))
  the_df
}

spread_cc <- function(the_df) {
  select(the_df, -from_since_dx, -thru_since_dx) %>%
  spread(drug, from_thru)
}

drug_plot <- function(the_df, lim, sz = 3) { 
  the_df_sum <-
    the_df %>% 
    group_by(drug, bmab) %>% 
    summarise(bd_count = n()) %>% 
    filter(bd_count > 10) %>% 
    modify_at("bd_count", 
              function(x) {
                the_x <- trimws(x,"both")
                if (all(grepl("[0-9]", unlist(strsplit(the_x, split = ""))))) {
                  prettyNum(round(as.numeric(the_x),2), big.mark = ",")
                } else {
                  x
                }
              }) 
  
  bv_the_df  <- the_df_sum %>% filter(bmab == "Bevacizumab")
  the_df     <- filter(the_df, drug %in% bv_the_df[["drug"]])
  the_df_sum <- filter(the_df_sum, drug %in% bv_the_df[["drug"]])
  
  g <- 
    ggplot(the_df, 
           aes(
             y = from_since_dx, 
             x = drug,
             color = drug)
    ) + 
    guides(color = "none") + 
    coord_flip() + 
    labs(y = "Days after primary diagnosis", x = "")
  
  g + 
    geom_boxplot(coef = 3, outlier.alpha = 0) +
    scale_y_continuous(breaks = seq(0, lim + 50, 365), limits = c(0, lim + 50)) + 
    geom_text(data = the_df_sum, 
              aes(x = drug, 
                  y = lim, 
                  label = bd_count), 
              color = "black", 
              size = sz
    ) + 
    theme(
      axis.text.x = element_text(color = "black"),
      axis.text.y = element_text(color = "black")
    ) 
}
