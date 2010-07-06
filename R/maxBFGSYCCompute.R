maxBFGSYCCompute <- function(fn,
                         start, print.level=0,
                  tol=1e-6, reltol=sqrt(.Machine$double.eps),
                  gradtol=1e-6, steptol=1e-10,
                  iterlim=150,
                             finalHessian=TRUE,
                  fixed=NULL,
                  ...) {
   ##                       iterlim = 2000,
   ##                       tol = 1E-06,
   ## This function is originally developed by Yves Croissant (and placed in 'mlogit' package).
   ## Fitted for 'maxLik' by Ott Toomet
   ## 
   ## BFGS maximisation, implemented by Yves Croissant
   ## Parameters:
   ## fn          - the function to be minimized.  Returns either scalar or
   ##               vector value with possible attributes 
   ##               constPar and newVal
   ##               fn must return the value with attribute 'gradient'
   ##               (and also attribute 'hessian' if it should be returned)
   ##               fn must have an argument sumObs
   ## start       - initial parameter vector (eventually w/names)
   ## steptol     - minimum step size
   ## ...         - extra arguments for fn()
   ## The stopping criteria
   ## tol         - maximum allowed absolute difference between sequential values
   ## reltol      - maximum allowed reltive difference (stops if < reltol*(abs(fn) + reltol)
   ## gradtol     - maximum allowed norm of gradient vector
   ## iterlim     - maximum # of iterations
   ## finalHessian  include final Hessian?  As computing final hessian does not carry any extra penalty for NR method, this option is
   ##               mostly for compatibility reasons with other maxXXX functions.
   ##               TRUE/something else  include
   ##               FALSE                do not include
   ## fixed       - a logical vector -- which parameters are taken as fixed.
   ##               Other paramters are treated as variable (free).
   ##
   ## RESULTS:
   ## a list of class "maxim":
   ## maximum     function value at maximum
   ## estimate    the parameter value at maximum
   ## gradient        gradient
   ## hessian         Hessian
   ## code        integer code of success:
   ##             1 - gradient close to zero
   ##             2 - successive values within tolerance limit
   ##             3 - could not find a higher point (step error)
   ##             4 - iteration limit exceeded
   ##             100 - initial value out of range
   ## message     character message describing the code
   ## last.step   only present if code == 3 (step error).  A list with following components:
   ##             theta0    - parameter value which led to the error
   ##             f0        - function value at these parameter values
   ##             climb     - the difference between theta0 and the new approximated parameter value (theta1)
   ##             fixed     - logical vector, which parameters are constant (fixed, inactive, non-free)
   ## fixed       logical vector, which parameters were treated as constant (fixed, inactive, non-free)
   ## iterations  number of iterations
   ## type        "BFGS-YC maximisation"
   ## 

   ##
   maxim.type <- "BFGS-YC maximization"
  param <- start
   nimed <- names(start)
   nParam <- length(param)
   ##
  chi2 <- 1E+10
  iter <- 0
  # eval a first time the function, the gradient and the hessian
  x <- sumKeepAttr( fn( param, fixed = fixed, sumObs = FALSE,
      returnHessian = FALSE, ... ) )
            # sum of log-likelihood value but not sum of gradients
   if (print.level > 0)
    cat(paste("Initial value of the function :", as.numeric(x), "\n"))
   if(is.na(x)) {
      result <- list(code=100, message=maximMessage("100"),
                     iterations=0,
                     type=maxim.type)
      class(result) <- "maxim"
      return(result)
   }
   if(is.infinite(x) & (x > 0)) {
                                        # we stop at +Inf but not at -Inf
      result <- list(code=5, message=maximMessage("5"),
                     iterations=0,
                     type=maxim.type)
      class(result) <- "maxim"
      return(result)
   }
   ##
   ## gradient by individual observations, used for BHHH approximation of initial Hessian.
   ## If not supplied by observations, we use the summed gradient.
   gri <- attr( x, "gradient" )
   gr <- sumGradients( gri, nParam = length( param ) )
   if(print.level > 2) {
      cat("Initial gradient value:\n")
      print(gr)
   }
   if(any(is.na(gr[!fixed]))) {
      stop("NA in the initial gradient")
   }
   if(any(is.infinite(gr[!fixed]))) {
      stop("Infinite initial gradient")
   }
   if(length(gr) != nParam) {
      stop( "length of gradient (", length(gr),
         ") not equal to the no. of parameters (", nParam, ")" )
   }
   ## initial approximation for inverse Hessian
   if(observationGradient(gri, length(param))) {
      invHess <- solve(crossprod(gri[,!fixed]))
                           # initial approximation of inverse Hessian (as in BHHH), if possible
   }
   else
       invHess <- 1e-5*diag(1, nrow=length(gr))
                           # if not possible (Is this OK?)
   if( print.level > 1) {
      cat( "----- Initial parameters: -----\n")
      cat( "fcn value:",
      as.vector(x), "\n")
      a <- cbind(start, gr, as.integer(!fixed))
      dimnames(a) <- list(nimed, c("parameter", "initial gradient",
                                          "free"))
      print(a)
   }
   samm <- NULL
                           # structure for too low 'step' value
  direction <- rep(0, nParam)
   repeat {
      if( iter >= iterlim) {
         code <- 4; break
      }
    step <- 1
     direction[!fixed] <- -as.vector(invHess %*% gr[!fixed])
      iter <- iter + 1
    oldx <- x
     oldgr <- gr
    oldparam <- param
    param[!fixed] <- oldparam[!fixed] - step * direction[!fixed]
    x <- sumKeepAttr( fn( param, fixed = fixed, sumObs = FALSE,
      returnHessian = FALSE, ... ) )
            # sum of log-likelihood value but not sum of gradients
    while((is.na(x) | x < oldx) & step > steptol) {
       step <- step/2
       if(print.level > 2) {
          cat("function values: old ", oldx, ", new ", x, ", difference ", x - oldx, " -> step ", step, "\n", sep="")
          if(print.level > 3) {
             resdet <- cbind(param = param, gradient = gr, direction=direction, active=!fixed)
             print(resdet)
          }
       }
       param[!fixed] <- oldparam[!fixed] - step * direction[!fixed]
       x <- sumKeepAttr( fn( param, fixed = fixed, sumObs = FALSE,
         returnHessian = FALSE, ... ) )
            # sum of log-likelihood value but not sum of gradients
    }
    if(step < steptol) {
                           # we did not find a better place to go...
       samm <- list(theta0=oldparam, f0=oldx, climb=direction)
    }
     gri <- attr( x, "gradient" )
                           # observation-wise gradient.  We only need it in order to compute the BHHH Hessian, if asked so.
     gr <- sumGradients( gri, nParam = length( param ) )
      incr <- step * direction
      y <- gr - oldgr
     invHess <- invHess +
        outer( incr[!fixed], incr[!fixed]) *
          (sum(y[!fixed] * incr[!fixed]) +
           as.vector( t(y[!fixed]) %*% invHess %*% y[!fixed])) / sum(incr[!fixed] * y[!fixed])^2 -
             (invHess %*% outer(y[!fixed], incr[!fixed])
              + outer(incr[!fixed], y[!fixed]) %*% invHess)/
                  sum(incr[!fixed] * y[!fixed])
      chi2 <- -  crossprod(direction[!fixed], oldgr[!fixed])
    if (print.level > 0){
       cat("--- iteration ", iter, ", step = ",step,
                      ", lnL = ", x,", chi2 = ",
           chi2,"\n",sep="")
       if (print.level > 1){
          resdet <- cbind(param = param, gradient = gr, direction=direction, active=!fixed)
          print(resdet)
          if(print.level > 3) {
             cat("Approximated Hessian:\n")
             print(invHess)
          }
          cat("--------------------------------------------\n")
       }
    }
      if( step < steptol) {
         code <- 3; break
      }
      if( sqrt( t(gr[!fixed])%*%gr[!fixed]) < gradtol) {
         code <-1; break
      }
      if(x - oldx < tol) {
         code <- 2; break
      }
      if(x - oldx < reltol*(x + reltol)) {
         code <- 2; break
      }
      if(is.infinite(x) & x > 0) {
         code <- 5; break
      }
  }
   if( print.level > 0) {
      cat( "--------------\n")
      cat( maximMessage( code), "\n")
      cat( iter, " iterations\n")
      cat( "estimate:", param, "\n")
      cat( "Function value:", x, "\n")
   }
   if( is.matrix( gr ) ) {
      if( dim( gr )[ 1 ] == 1 ) {
         gr <- gr[ 1, ]
      }
   }
   names(gr) <- names(param)
   # calculate (final) Hessian
   if(tolower(finalHessian) == "bhhh") {
      if(observationGradient(gri, length(param))) {
          hessian <- -t(gri) %*% gri
          attr(hessian, "type") <- "BHHH"
      } else {
         hessian <- NULL
         warning("For computing the final Hessian by 'BHHH' method, the log-likelihood or gradient must be supplied by observations")
      }
   }
   else if(finalHessian) {
      hessian <- attr( fn( param, fixed = fixed, returnHessian = TRUE, ... ) ,
         "hessian" )
   }
   else
       hessian <- NULL
   if( !is.null( hessian ) ) {
      rownames( hessian ) <- colnames( hessian ) <- nimed
   }
   ## remove attributes from final value of objective (likelihood) function
   attributes( x )$gradient <- NULL
   attributes( x )$hessian <- NULL
   ##
   result <-list(
                  maximum = unname( drop( x ) ),
                  estimate=param,
                  gradient=gr,
                 hessian=hessian,
                  code=code,
                  message=maximMessage( code),
                  last.step=samm,
                                        # only when could not find a
                                        # lower point
                  fixed=fixed,
                  iterations=iter,
                  type=maxim.type)
   if(observationGradient(gri, length(param))) {
       colnames( gri ) <- names( param )
       result$gradientObs <- gri
   }
   class(result) <- c("maxim", class(result))
   invisible(result)
}
