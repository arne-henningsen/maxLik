glance.maxLik <- function (x, ...) {
  require_tibble_package()

  ll <- logLik(x)
  nobs <- tryCatch(nObs(x), error = function(e) NA)
                           # nobs = NA in case of error
  ret <- tibble::tibble(
           df     = attr(ll, "df"),
           logLik = as.numeric(ll),
           AIC    = AIC(x),
           nobs   = nobs
         )
  
  ret
}
