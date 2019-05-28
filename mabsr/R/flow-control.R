


switches <- function(..., .force = FALSE) {
  sws <- dots_list(...)
  for (i in seq_along(sws)) {
    if (!exists(names(sws)[i], envir = .GlobalEnv) | .force) {
      assign(names(sws)[i], unlist(sws[i]), envir = .GlobalEnv)
    }
  }
}

only_if_real <- function(an_expr, filepath, ...) {
  library(rlang)
  else_statement <- 
    substitute({
      thing_nms <- map_chr(exprs(...), as.character)
      things <- read_rds(filepath)
      for (i in seq_along(thing_nms)) {
        the_nm <- thing_nms[i]
        assign(the_nm, things[[the_nm]], envir = .GlobalEnv)
      }
    })
  if (exists("for_real")) {
    if (for_real) {
      eval(an_expr)
      write_rds(..., filepath)
    } else {
      eval(else_statement)
    }
  } else {
    eval(else_statement)
  }
}

