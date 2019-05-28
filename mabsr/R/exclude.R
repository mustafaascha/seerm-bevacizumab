
exclude <- function(a_df) {
  filter(
    a_df, 
    !grepl("^Small Cell Carcinoma, Nos$", histo) & 
      !grepl("^8", age_cut) & 
      dx_year >= 2010 & 
      seq_num %in% c("00", "01") & 
      !grepl("situ", d_ssg00)
  )
}



