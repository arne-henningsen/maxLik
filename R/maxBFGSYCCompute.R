maxBFGSYCCompute <- function(fn, grad=NULL, hess=NULL,
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
   ##               If fn returns the value with attributes 'gradient'
   ##               and 'hessian', those are used instead of calculatin
   ##               new gradient and hessian values
   ## grad        - gradient function (numeric used if missing).  Must return either
   ##               * vector, length=nParam
   ##               * matrix, dim=c(nObs, 1).  Treated as vector
   ##               * matrix, dim=c(M, nParam), where M is arbitrary.  In this case the
   ##                 rows are simply summed (useful for maxBHHH).
   ## hess        - hessian function (numeric used if missing)
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
   ##             activePar - logical vector, which parameters are active (not constant)
   ## activePar   logical vector, which parameters were treated as free (resp fixed)
   ## iterations  number of iterations
   ## type        "BFGS-YC maximisation"
   ## 
   func <- function(theta, sumObs = TRUE, ...) {
      f <- fn(theta, ...)
      if( sumObs ) {
         fAttr <- attributes(f)
         f <- sum(f)
         mostattributes(f) <- fAttr
      }
      return( f )
   }
   gradient <- function(theta, sumObs = TRUE,
                        suppliedValue=NULL, ...) {
      ## suppliedValue: use gradient value supplied from elsewhere
      ##           (attribute to fn) and only do sanity checks
      if(is.null(gr <- suppliedValue)) {
         if(!is.null(grad)) {  # use user-supplied if present
            gr <- grad(theta, ...)
         } else {
            gr <- numericGradient(f = func, t0 = theta,
                                  activePar=activePar, sumObs = sumObs, ...)
                           # Note we need nObs rows x nParam cols
         }
      }
      ## Now check if the gradient is vector or matrix...
      if(!sumObs) {
         ## return (preferably) by observations, ensure it will be a matrix
         if(observationGradient(gr, length(theta))) {
            gr <- as.matrix(gr)
         }
      }
      else {
         ## We need just summed gradient
         if( !is.null(dim(gr))) {
            gr <- colSums(gr)
         } else {
            ## ... or vector if only one parameter
            if(length(gr) > nParam ) {
               gr <- sum(gr)
            }
         }
      }
      ## Why do we need this?
      if( is.null( dim( gr ) ) ) {
         gr[ !activePar ] <- NA
      } else {
         gr[ , !activePar ] <- NA
      }
      return(gr)
   }
   hessian <- function(theta, activePar=activePar,
                       suppliedValue=NULL, ...) {
      ## Hessian is only used for the final hessian (if asked)
      ## 
      ## suppliedValue: use hessian value supplied from elsewhere
      ##           (attribute to fn) and only do sanity checks
      ##
      ## Note: a call to hessian must follow a call to gradient using
      ## /exactly the same/ parameter values in this program. 
      ## This ensures compatibility with maxBHHH.  This warning does not
      ## apply for user programs.
      if(is.null(h <- suppliedValue)) {
         if(!is.null(hess)) {
            h <- as.matrix(hess(theta, ...))
         } else {
            h <- numericHessian( f = func, grad = gradient, t0 = theta,
                                activePar=activePar, ...)
         }
      }
      if((dim(h)[1] != nParam) | (dim(h)[2] != nParam)) {
         stop("Wrong hessian dimension.  Needed ", nParam, "x", nParam,
              " but supplied ", dim(h)[1], "x", dim(h)[2])
      }
      h[ !activePar, ] <- NA
      h[ , !activePar ] <- NA 
      return( h )
   }
   ##
   maxim.type <- "BFGS-YC maximization"
   argNames <- c( "fn", "grad", "hess", "start", "print.level",
      "tol", "reltol", "gradtol", "steptol", 
      "iterlim", "activePar", "fixed" )
   checkFuncArgs( fn, argNames, "fn", "maxNR" )
   if( !is.null( grad ) ) {
      checkFuncArgs( grad, argNames, "grad", "maxNR" )
   }
  param <- start
   nimed <- names(start)
   nParam <- length(param)
   ## establish the active parameters.  Internally, we just use 'activePar'
   activePar <- !fixed
   rm( fixed )
   ##
  chi2 <- 1E+10
  iter <- 0
  # eval a first time the function, the gradient and the hessian
  x <- func(param, ...)
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
   gri <- gradient(param, suppliedValue=attr(x, "gradient"), sumObs=FALSE, ...)
   if(observationGradient(gri, length(param))) 
       gr <- colSums(gri)
   else {
      gr <- gri
   }
   if(print.level > 2) {
      cat("Initial gradient value:\n")
      print(gr)
   }
   if(any(is.na(gr[activePar]))) {
      stop("NA in the initial gradient")
   }
   if(any(is.infinite(gr[activePar]))) {
      stop("Infinite initial gradient")
   }
   if(length(gr) != nParam) {
      stop( "length of gradient (", length(gr),
         ") not equal to the no. of parameters (", nParam, ")" )
   }
   ## initial approximation for inverse Hessian
   if(observationGradient(gri, length(param))) {
      invHess <- solve(crossprod(gri[,activePar]))
                           # initial approximation of inverse Hessian (as in BHHH), if possible
   }
   else
       invHess <- 1e-5*diag(1, nrow=length(gr))
                           # if not possible (Is this OK?)
   if( print.level > 1) {
      cat( "----- Initial parameters: -----\n")
      cat( "fcn value:",
      as.vector(x), "\n")
      a <- cbind(start, gr, as.integer(activePar))
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
     direction[activePar] <- -as.vector(invHess %*% gr[activePar])
      iter <- iter + 1
    oldx <- x
     oldgr <- gr
    oldparam <- param
    param[activePar] <- oldparam[activePar] - step * direction[activePar]
    x <- func(param, ...)
    while((is.na(x) | x < oldx) & step > steptol) {
       step <- step/2
       if(print.level > 2) {
          cat("function values: old ", oldx, ", new ", x, ", difference ", x - oldx, " -> step ", step, "\n", sep="")
          if(print.level > 3) {
             resdet <- cbind(param = param, gradient = gr, direction=direction, active=activePar)
             print(resdet)
          }
       }
       param[activePar] <- oldparam[activePar] - step * direction[activePar]
       x <- func(param, ...)
    }
    if(step < steptol) {
                           # we did not find a better place to go...
       samm <- list(theta0=oldparam, f0=oldx, climb=direction)
    }
     gri <- gradient(param, suppliedValue=attr(x, "gradient"), sumObs=FALSE, ...)
                           # observation-wise gradient.  We only need it in order to compute the BHHH Hessian, if asked so.
      if(observationGradient(gri, length(param)))
          gr <- colSums(gri)
      else
          gr <- gri
      incr <- step * direction
      y <- gr - oldgr
     invHess <- invHess +
        outer( incr[activePar], incr[activePar]) *
          (sum(y[activePar] * incr[activePar]) +
           as.vector( t(y[activePar]) %*% invHess %*% y[activePar])) / sum(incr[activePar] * y[activePar])^2 -
             (invHess %*% outer(y[activePar], incr[activePar])
              + outer(incr[activePar], y[activePar]) %*% invHess)/
                  sum(incr[activePar] * y[activePar])
      chi2 <- -  crossprod(direction[activePar], oldgr[activePar])
    if (print.level > 0){
       cat("--- iteration ", iter, ", step = ",step,
                      ", lnL = ", x,", chi2 = ",
           chi2,"\n",sep="")
       if (print.level > 1){
          resdet <- cbind(param = param, gradient = gr, direction=direction, active=activePar)
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
      if( sqrt( t(gr[activePar])%*%gr[activePar]) < gradtol) {
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
      if(observationGradient(gri, length(param)))
          hessian <- -t(gri) %*% gri
      else {
         hessian <- logLikHess( param, fnOrig = fn,  gradOrig = grad,
                                hessOrig = hess, ... )
         warning("For computing Hessian by 'BHHH' method, the log-likelihood or gradient must be supplied by observations")
      }
   }
   else if(finalHessian) {
      hessian <- logLikHess( param, fnOrig = fn,  gradOrig = grad,
                            hessOrig = hess, ... )
   }
   else
       hessian <- NULL
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
                  activePar=activePar,
                  iterations=iter,
                  type=maxim.type)
   if(observationGradient(gri, length(param)))
       result$gradientObs <- gri
   class(result) <- c("maxim", class(result))
   invisible(result)
}
