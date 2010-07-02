maxNRCompute <- function(fn, grad=NULL, hess=NULL,
                         start, print.level=0,
                  tol=1e-8, reltol=sqrt(.Machine$double.eps),
                  gradtol=1e-6, steptol=1e-10,
                  lambdatol=1e-6,
                  qrtol=1e-10,
                  iterlim=150,
                         finalHessian=TRUE,
                  fixed=NULL,
                  ...) {
   ## Newton-Raphson maximisation
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
   ## lambdatol   - max lowest eigenvalue when forcing pos. definite H
   ## qrtol       - tolerance for qr decomposition
   ## ...         - extra arguments for fn()
   ## The stopping criteria
   ## tol         - maximum allowed absolute difference between sequential values
   ## reltol      - maximum allowed reltive difference (stops if < reltol*(abs(fn) + reltol)
   ## gradtol     - maximum allowed norm of gradient vector
   ## 
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
   ## type        "Newton-Raphson maximisation"
   
   argNames <- c( "fn", "grad", "hess", "start", "print.level",
      "tol", "reltol", "gradtol", "steptol", "lambdatol", "qrtol",
      "iterlim", "activePar", "fixed" )
   checkFuncArgs( fn, argNames, "fn", "maxNR" )
   if( !is.null( grad ) ) {
      checkFuncArgs( grad, argNames, "grad", "maxNR" )
   }
   if( !is.null( hess ) ) {
      checkFuncArgs( hess, argNames, "hess", "maxNR" )
   }
   max.eigen <- function( M) {
      ## return maximal eigenvalue of (symmetric) matrix
      val <- eigen(M, symmetric=TRUE, only.values=TRUE)$values
      val[1]
      ## L - eigenvalues in decreasing order, [1] - biggest in abs value
   }
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
                       suppliedValue=NULL, gradAttr = FALSE, ...) {
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
            if( gradAttr ) {
               gradFunc <- function( theta, ... ) {
                  gr <- attr( func( theta, ... ), "gradient" )
                  if( !is.null( dim( gr ) ) ) {
                     gr <- colSums(gr)
                  } else {
                     if( length( gr ) > nParam ) {
                        gr <- sum( gr )
                     }
                  }
                  return( gr )
               }
            } else {
               gradFunc <- gradient
            }
            h <- numericHessian( f = func, grad = gradFunc, t0 = theta,
                                activePar=activePar, ...)
         }
      }
      if((dim(h)[1] != nParam) | (dim(h)[2] != nParam)) {
         stop("Wrong hessian dimension.  Needed ", nParam, "x", nParam,
              " but supplied ", dim(h)[1], "x", dim(h)[2])
      }
      ## Why do we need this?
      h[ !activePar, ] <- NA
      h[ , !activePar ] <- NA 
      return( h )
   }
   ## -------------------------------------------------
   maxim.type <- "Newton-Raphson maximisation"
   nimed <- names(start)
   nParam <- length(start)
   ## establish the active parameters.  Internally, we just use 'activePar'
   activePar <- !fixed
   rm( fixed )

   samm <- NULL
   I <- diag(rep(1, nParam))
                           # I is unit matrix
   start1 <- start
   iter <- 0
   f1 <- func(start1, ...)
   if(print.level > 2) {
      cat("Initial function value:", f1, "\n")
   }
   if(is.na( f1)) {
      result <- list(code=100, message=maximMessage("100"),
                     iterations=0,
                     type=maxim.type)
      class(result) <- "maxim"
      return(result)
   }
   if(is.infinite( f1) & (f1 > 0)) {
                                        # we stop at +Inf but not at -Inf
      result <- list(code=5, message=maximMessage("5"),
                     iterations=0,
                     type=maxim.type)
      class(result) <- "maxim"
      return(result)
   }
   G1 <- gradient(start1, suppliedValue=attr(f1, "gradient"), ...)
   if(print.level > 2) {
      cat("Initial gradient value:\n")
      print(G1)
   }
   if(any(is.na(G1[activePar]))) {
      stop("NA in the initial gradient")
   }
   if(any(is.infinite(G1[activePar]))) {
      stop("Infinite initial gradient")
   }
   if(length(G1) != nParam) {
      stop( "length of gradient (", length(G1),
         ") not equal to the no. of parameters (", nParam, ")" )
   }
   H1 <- hessian(start1, activePar=activePar,
                 suppliedValue=attr(f1, "hessian"),
                 gradAttr = !is.null( attr( f1, "gradient" ) ), ... )
   if(any(is.na(H1[activePar, activePar]))) {
      stop("NA in the initial Hessian")
   }
   if(any(is.infinite(H1))) {
      stop("Infinite initial Hessian")
   }
   if( print.level > 1) {
      cat( "----- Initial parameters: -----\n")
      cat( "fcn value:",
      as.vector(f1), "\n")
      a <- cbind(start, G1, as.integer(activePar))
      dimnames(a) <- list(nimed, c("parameter", "initial gradient",
                                          "free"))
      print(a)
      cat( "Condition number of the (active) hessian:",
          kappa( H1[activePar, activePar]), "\n")
      if( print.level > 3) {
         print( H1)
      }
   }
   repeat {
      if( iter >= iterlim) {
         code <- 4; break
      }
      iter <- iter + 1
      lambda <- 0
      start0 <- start1
      f0 <- f1
      G0 <- G1
      if(any(is.na(G0[activePar]))) {
         stop("NA in gradient (at the iteration start)")
      }
      H0 <- H1
      if(any(is.na(H0[activePar, activePar]))) {
         stop("NA in Hessian (at the iteration start)")
      }
      step <- 1
      H <- H0
      ## check whether hessian is positive definite
      while((me <- max.eigen( H[activePar,activePar,drop=FALSE])) >= -lambdatol |
         (qRank <- qr(H[activePar,activePar], tol=qrtol)$rank) < sum(activePar)) {
                                        # maximum eigenvalue -> negative definite
                                        # qr()$rank -> singularity
         lambda <- abs(me) + lambdatol + min(abs(diag(H)[activePar]))/1e7
                           # The third term corrects numeric singularity.  If diag(H) only contains large values,
                           # (H - (a small number)*I) == H because of finite precision
         H <- H - lambda*I
                                        # how to make it better?
      }
      amount <- vector("numeric", nParam)
      amount[activePar] <- qr.solve(H[activePar,activePar,drop=FALSE],
                                    G0[activePar], tol=qrtol)
      start1 <- start0 - step*amount
      f1 <- func(start1, ...)
      ## Are we requested to fix some of the parameters?
      constPar <- attr(f1, "constPar")
      if(!is.null(constPar)) {
         if(any(is.na(constPar))) {
            stop("NA in the list of constants")
         }
         activePar <- rep(TRUE, nParam)
         activePar[constPar] <- FALSE
      }
      ## Are we asked to write in a new value for some of the parameters?
      if(is.null(newVal <- attr(f1, "newVal"))) {
         ## no ...
         while( is.na( f1) | ( ( f1 < f0) & ( step >= steptol))) {
                                        # We end up in a NA or a higher value.
                                        # try smaller step
            step <- step/2
            if(print.level > 2) {
               cat("function value difference", f1 - f0, "-> step", step, "\n")
            }
            start1 <- start0 - step*amount
            f1 <- func(start1, ...)
            ## Find out the constant parameters -- these may be other than
            ## with full step
            constPar <- attr(f1, "constPar")
            if(!is.null(constPar)) {
               if(any(is.na(constPar))) {
                  stop("NA in the list of constants")
               }
               activePar[constPar] <- FALSE
               ## Any new values requested?
               if(!is.null(newVal <- attr(f1, "newVal"))) {
                  ## Yes.  Write them to parameters and go for
                  ## next iteration
                  start1[newVal$index] <- newVal$val
                  break;
               }
            }
         }
         if(step < steptol) {
            # we did not find a better place to go...
            start1 <- start0
            f1 <- f0
            samm <- list(theta0=start0, f0=f0, climb=amount)
         }
      } else {
         ## Yes, indeed.  New values given to some of the params.
         ## Note, this may result in a lower function value,
         ## hence we do not check f1 > f0
         start1[newVal$index] <- newVal$val
         print(start1)
      }
      G1 <- gradient(start1, suppliedValue=attr(f1, "gradient"), ...)
      if(any(is.na(G1[activePar]))) {
         cat("Iteration", iter, "\n")
         cat("Parameter:\n")
         print(start1)
         if(length(G1) < 30) {
            cat("Gradient:\n")
            print(G1)
         }
         stop("NA in gradient")
      }
      if(any(is.infinite(G1))) {
         code <- 6; break;
      }
      H1 <- hessian(start1, activePar=activePar,
                    suppliedValue=attr(f1, "hessian"), 
                    gradAttr = !is.null( attr( f1, "gradient" ) ), ... )
      if( print.level > 1) {
        cat( "-----Iteration", iter, "-----\n")
      }
      if(any(is.infinite(H1))) {
         code <- 7; break
      }
      if(print.level > 2) {
         cat( "lambda ", lambda, " step", step, " fcn value:",
            formatC(as.vector(f1), digits=8, format="f"),  "\n")
         a <- cbind(amount, start1, G1, as.integer(activePar))
         dimnames(a) <- list(names(start0), c("amount", "new param",
                                             "new gradient", "active"))
         print(a)
         if( print.level > 3) {
            cat("Hessian\n")
            print( H1)
         }
         if(!any(is.na(H1[activePar, activePar]))) {
            cat( "Condition number of the hessian:",
                kappa(H1[activePar,activePar,drop=FALSE]), "\n")
         }
      }
      if( step < steptol) {
         code <- 3; break
      }
      if( sqrt( t(G1[activePar])%*%G1[activePar]) < gradtol) {
         code <-1; break
      }
      if(is.null(newVal) & f1 - f0 < tol) {
         code <- 2; break
      }
      if(is.null(newVal) & f1 - f0 < reltol*(f1 + reltol)) {
         code <- 2; break
      }
      if(is.infinite(f1) & f1 > 0) {
         code <- 5; break
      }
   }
   if( print.level > 0) {
      cat( "--------------\n")
      cat( maximMessage( code), "\n")
      cat( iter, " iterations\n")
      cat( "estimate:", start1, "\n")
      cat( "Function value:", f1, "\n")
   }
   names(start1) <- nimed
   G1 <- gradient(start1, suppliedValue=attr(f1, "gradient"),
                  sumObs=FALSE, ...)
   if(observationGradient(G1, length(start1))) {
      gradientObs <- G1
      colnames( gradientObs ) <- nimed 
      G1 <- colSums(as.matrix(G1 ))
   }
   else {
      gradientObs <- NULL
   }
   names( G1 ) <- nimed
   ## calculate (final) Hessian
   if(tolower(finalHessian) == "bhhh" & !is.null(gradientObs)) {
      hessian <- -t(gradientObs) %*% gradientObs
      attr(hessian, "type") <- "BHHH"
   }
   else if(finalHessian != FALSE) {
      if(tolower(finalHessian) == "bhhh")
          warning("Final BHHH Hessian: gradient by observations not available.  Using Hessian matrix")
      hessian <- hessian( start1, activePar=activePar,
                         suppliedValue=attr(f1, "hessian"), 
                         gradAttr = !is.null( attr( f1, "gradient" ) ), ... )
   }
   else
       hessian <- NULL
   rownames( hessian ) <- colnames( hessian ) <- nimed
   ## remove attributes from final value of objective (likelihood) function
   attributes( f1 )$gradient <- NULL
   attributes( f1 )$hessian <- NULL
   ##
   result <-list(
                  maximum = unname( drop( f1 ) ),
                  estimate=start1,
                  gradient=drop(G1),
                 hessian=hessian,
                  code=code,
                  message=maximMessage( code),
                  last.step=samm,
                                        # only when could not find a
                                        # lower point
                  activePar=activePar,
                  iterations=iter,
                  type=maxim.type)
   if( exists( "gradientObs" ) ) {
      result$gradientObs <- gradientObs
   }

   class(result) <- c("maxim", class(result))
   invisible(result)
}

returnCode.maxim <- function(x, ...)
    x$code
