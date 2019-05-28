




srv_frm <- function(prd) {
  if (length(prd) > 1) prd <- paste(prd, collapse = "+")
  as.formula(paste('Surv(srv_time, event = alive == "Dead") ~ ', prd))
}

rm_pr <- function(vr) mabsr::prds$candidates[-grep(vr, mabsr::prds$candidates)]





