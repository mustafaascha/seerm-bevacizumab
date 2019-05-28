
dig <- function() {
  map(2:4, function(a) function(z) sprintf(paste("%.", a, "f", sep = ""), z)
  ) %>% 
  set_names(c("two", "three", "four"))
}






