



ex_rmd <- function(rmd) {
  e <- grep('^```$', rmd)
  s <- grep('```\\{', rmd)
  stopifnot(length(e) == length(s))
  r <- list()
  for (i in seq_along(e)) {
    r[i] <- list(rmd[seq(s[i] + 1, e[i] - 1)])
  }
  r
}
