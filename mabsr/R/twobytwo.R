


tidy_tbt <- function(tbt) {
  data.frame(
    cbind(
      tbt[["measures"]], 
      tbt[["p.value"]]
    ),
    stringsAsFactors = FALSE
  )
}



two_by_twos <- function(a_df) {

  two_by_two <- 
    partial(                               
      Epi::twoby2,                         
      outcome = a_df[["bmab01"]],
      print = FALSE                        
    )

  grep(
    "bmab01|id|propensity_score|brain_mets|srv.*|stat", 
    names(a_df), 
    value = T, 
    invert = T
   )  %>% 
   (function(vrs) {
     map(
      vrs,
      function(vr_nm) {
 
       predictor <- a_df[[vr_nm]]

       lu_p <- length(unique(predictor))
 
       if (lu_p > 10 | lu_p == 1) {
       } else if (lu_p > 2) {
         map(
           unique(predictor), 
           function(u_p) {
             factor(predictor == u_p, levels = 0:1, labels = c("None", u_p)
             ) %>% 
             two_by_two()
         })
       } else { two_by_two(predictor) }
     }) %>% 
     set_names(vrs)

  })
}
