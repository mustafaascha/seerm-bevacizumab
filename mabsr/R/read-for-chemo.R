

read <- function(fp, ...) {
  read_csv(fp, col_types = c_cols(), ...)
}

bind_clean <- function(the_dfs) {
  the_df <- reduce(the_dfs, bind_rows)
  the_df[["drug"]] <- 
    with(the_df, 
      ifelse(!is.na(c_desc), c_desc, 
        ifelse(!is.na(ndc), ndc, 
          ifelse(!is.na(is_chemo) & nchar(is_chemo > 1), 
            is_chemo, 
            NA
          )
        )
      )
    )
  the_df[["is_chemo"]][nchar(the_df[["is_chemo"]]) > 1] <- NA
  the_df[,c("c_desc", "ndc")] <- NULL
  the_df %>% 
    select(patient_id, from, thru, drug, is_chemo) %>% 
    drop_na()
}

unpack_chemo <- function(the_df) {
  separate(the_df, 
    "hcpcs", 
    c("c_desc", "c_cat", "c_type", "tx_type", "is_chemo"), 
    sep = "\\*\\*\\*", 
    fill = "left"
  )
}

is_packed <- function(z) {
  any(grepl("\\*\\*\\*", z[["hcpcs"]]))
}








