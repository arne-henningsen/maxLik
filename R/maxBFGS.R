maxBFGS <- function(fn, grad=NULL, hess=NULL,
                    start,
                    print.level=0,
                    iterlim=200,
                    constraints=NULL,
                    tol=1e-8, reltol=tol,
                    parscale=rep(1, length=length(start)),
                    ## sumt parameters
                    ...) {
   ## contraints    constraints to be passed to 'constrOptim'
   ## ...           further arguments to fn() and grad()

   result <- maxOptim( fn = fn, grad = grad, hess = hess,
      start = start, method = "BFGS", print.level = print.level,
      iterlim = iterlim, constraints = constraints,
      tol = tol, reltol = reltol, parscale = parscale,
      ... )

   return(result)
}
