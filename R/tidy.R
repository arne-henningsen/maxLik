
require_tibble_package <- function () {
  if (! requireNamespace("tibble", quietly = TRUE)) {
    stop("The `tibble` package must be installed to use tidy() or glance() methods")
  }
}


tidy.maxLik <- function (x,  ...) {
  require_tibble_package()
  
  s <- summary.maxLik(x)
  ret <- tibble::as_tibble(.estimate(s), rownames = "term")
  colnames(ret) <- c("term", "estimate", "std.error", "statistic", "p.value")
  
  ret
}

tidy.maxim <- function(x, ...) {
   coefs <- .estimate(x)
   term <- names(coefs)
   if(is.null(term))
      term <- seq(along=coefs)
   tibble::tibble(term,
                  estimate=coefs,
                  gradient=.gradient(x)
                  )
}
