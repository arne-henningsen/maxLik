maxSANN <- function(fn, grad=NULL, hess=NULL,
                    start,
                    print.level=0,
                    iterlim=10000,
                    constraints = NULL,
                    tol=1e-8, reltol=tol, cand = NULL,
                    temp=10, tmax=10, parscale=rep(1, length=length(start)),
                    random.seed = 123, ... ) {
   ## ... : further arguments to fn()
   ##
   ## Note: grad and hess are for compatibility only, SANN uses only fn values

   # save seed of the random number generator
   if( exists( ".Random.seed" ) ) {
      savedSeed <- .Random.seed
   }

   # set seed for the random number generator (used by 'optim( method="SANN" )')
   set.seed( random.seed )

   # restore seed of the random number generator on exit
   # (end of function or error)
   if( exists( "savedSeed" ) ) {
      on.exit( assign( ".Random.seed", savedSeed, envir = sys.frame() ) )
   } else {
      on.exit( rm( .Random.seed, envir = sys.frame() ) )
   }

   result <- maxOptim( fn = fn, grad = grad, hess = hess,
      start = start, method = "SANN", print.level = print.level,
      iterlim = iterlim, constraints = constraints,
      tol = tol, reltol = reltol, parscale = parscale,
      temp = temp, tmax = tmax, random.seed = random.seed, cand = cand,
      ... )

   return(result)
}

