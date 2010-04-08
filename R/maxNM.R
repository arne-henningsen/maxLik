maxNM <- function(fn, grad=NULL, hess=NULL,
                  start, fixed = NULL,
                  print.level=0,
                  iterlim=500,
                  constraints=NULL,
                  tol=1e-8, reltol=tol,
                  parscale=rep(1, length=length(start)),
                  alpha=1, beta=0.5, gamma=2,
                  ...) {
   ## contraints    constraints to be passed to 'constrOptim'
   ## ... : further arguments to fn()
   ##
   ## Note: grad and hess are for compatibility only, SANN uses only fn values

   result <- maxOptim( fn = fn, grad = grad, hess = hess,
      start = start, method = "Nelder-Mead", fixed = fixed,
      print.level = print.level, iterlim = iterlim, constraints = constraints,
      tol = tol, reltol = reltol, parscale = parscale,
      alpha = alpha, beta = beta, gamma = gamma,
      ... )

   return(result)
}
