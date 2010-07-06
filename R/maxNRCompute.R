maxNRCompute <- function(fn,
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
   ##               fn must return the value with attributes 'gradient'
   ##               and 'hessian'
   ##               fn must have an argument sumObs
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
   ##             fixed     - logical vector, which parameters are constant (fixed, inactive, non-free)
   ## fixed       logical vector, which parameters were treated as constant (fixed, inactive, non-free)
   ## iterations  number of iterations
   ## type        "Newton-Raphson maximisation"
   
   max.eigen <- function( M) {
      ## return maximal eigenvalue of (symmetric) matrix
      val <- eigen(M, symmetric=TRUE, only.values=TRUE)$values
      val[1]
      ## L - eigenvalues in decreasing order, [1] - biggest in abs value
   }

   ## -------------------------------------------------
   maxim.type <- "Newton-Raphson maximisation"
   nimed <- names(start)
   nParam <- length(start)

   samm <- NULL
   I <- diag(rep(1, nParam))
                           # I is unit matrix
   start1 <- start
   iter <- 0
   f1 <- fn(start1, fixed = fixed, sumObs = TRUE, ...)
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
   G1 <- attr( f1, "gradient" )
   if(print.level > 2) {
      cat("Initial gradient value:\n")
      print(G1)
   }
   if(any(is.na(G1[!fixed]))) {
      stop("NA in the initial gradient")
   }
   if(any(is.infinite(G1[!fixed]))) {
      stop("Infinite initial gradient")
   }
   if(length(G1) != nParam) {
      stop( "length of gradient (", length(G1),
         ") not equal to the no. of parameters (", nParam, ")" )
   }
   H1 <- attr( f1, "hessian" )
   if(any(is.na(H1[!fixed, !fixed]))) {
      stop("NA in the initial Hessian")
   }
   if(any(is.infinite(H1))) {
      stop("Infinite initial Hessian")
   }
   if( print.level > 1) {
      cat( "----- Initial parameters: -----\n")
      cat( "fcn value:",
      as.vector(f1), "\n")
      a <- cbind(start, G1, as.integer(!fixed))
      dimnames(a) <- list(nimed, c("parameter", "initial gradient",
                                          "free"))
      print(a)
      cat( "Condition number of the (active) hessian:",
          kappa( H1[!fixed, !fixed]), "\n")
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
      if(any(is.na(G0[!fixed]))) {
         stop("NA in gradient (at the iteration start)")
      }
      H0 <- H1
      if(any(is.na(H0[!fixed, !fixed]))) {
         stop("NA in Hessian (at the iteration start)")
      }
      step <- 1
      H <- H0
      ## check whether hessian is positive definite
      while((me <- max.eigen( H[!fixed,!fixed,drop=FALSE])) >= -lambdatol |
         (qRank <- qr(H[!fixed,!fixed], tol=qrtol)$rank) < sum(!fixed)) {
                                        # maximum eigenvalue -> negative definite
                                        # qr()$rank -> singularity
         lambda <- abs(me) + lambdatol + min(abs(diag(H)[!fixed]))/1e7
                           # The third term corrects numeric singularity.  If diag(H) only contains large values,
                           # (H - (a small number)*I) == H because of finite precision
         H <- H - lambda*I
                                        # how to make it better?
      }
      amount <- vector("numeric", nParam)
      amount[!fixed] <- qr.solve(H[!fixed,!fixed,drop=FALSE],
                                    G0[!fixed], tol=qrtol)
      start1 <- start0 - step*amount
      f1 <- fn(start1, fixed = fixed, sumObs = TRUE, ...)
      ## Are we requested to fix some of the parameters?
      constPar <- attr(f1, "constPar")
      if(!is.null(constPar)) {
         if(any(is.na(constPar))) {
            stop("NA in the list of constants")
         }
         fixed <- rep(FALSE, nParam)
         fixed[constPar] <- TRUE
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
            f1 <- fn(start1, fixed = fixed, sumObs = TRUE, ...)
            ## Find out the constant parameters -- these may be other than
            ## with full step
            constPar <- attr(f1, "constPar")
            if(!is.null(constPar)) {
               if(any(is.na(constPar))) {
                  stop("NA in the list of constants")
               }
               fixed[constPar] <- TRUE
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
      G1 <- attr( f1, "gradient" )
      if(any(is.na(G1[!fixed]))) {
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
      H1 <- attr( f1, "hessian" )
      if( print.level > 1) {
        cat( "-----Iteration", iter, "-----\n")
      }
      if(any(is.infinite(H1))) {
         code <- 7; break
      }
      if(print.level > 2) {
         cat( "lambda ", lambda, " step", step, " fcn value:",
            formatC(as.vector(f1), digits=8, format="f"),  "\n")
         a <- cbind(amount, start1, G1, as.integer(!fixed))
         dimnames(a) <- list(names(start0), c("amount", "new param",
                                             "new gradient", "active"))
         print(a)
         if( print.level > 3) {
            cat("Hessian\n")
            print( H1)
         }
         if(!any(is.na(H1[!fixed, !fixed]))) {
            cat( "Condition number of the hessian:",
                kappa(H1[!fixed,!fixed,drop=FALSE]), "\n")
         }
      }
      if( step < steptol) {
         code <- 3; break
      }
      if( sqrt( t(G1[!fixed])%*%G1[!fixed]) < gradtol) {
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
   F1 <- fn( start1, fixed = fixed, sumObs = FALSE, ... )
   G1 <- attr( F1, "gradient" )
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
   if(tolower(finalHessian) == "bhhh") {
      if(!is.null(gradientObs)) {
         hessian <- - crossprod( gradientObs )
         attr(hessian, "type") <- "BHHH"
      } else {
         hessian <- NULL
         warning("For computing the final Hessian by 'BHHH' method, the log-likelihood or gradient must be supplied by observations")
      }
   } else if( finalHessian != FALSE ) {
      hessian <- attr( F1, "hessian" )
   } else {
       hessian <- NULL
   }
   if( !is.null( hessian ) ) {
      rownames( hessian ) <- colnames( hessian ) <- nimed
   }

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
                  fixed=fixed,
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
