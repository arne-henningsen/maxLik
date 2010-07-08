### SUMT function borrowed from 'clue' package
### 
### Adapted for linear constraints
sumt <- function(fn, grad=NULL, hess=NULL,
                 start,
                 maxRoutine, constraints, 
                 SUMTTol = sqrt(.Machine$double.eps),
                           # difference between estimates for successive outer iterations
                 SUMTPenaltyTol = sqrt(.Machine$double.eps),
                           # maximum allowed penalty
                 SUMTQ = 10,
                 SUMTRho0 = NULL,
                 print.level=0,
                 SUMTMaxIter=100,
                 ...) {
   ## constraints    list w/components eqA and eqB.  Maximization will
   ##                be performed wrt to the constraint
   ##                A %*% theta + B = 0
   ##                The user must ensure the matrices are in correct
   ##                form
   ## maxSUMTiter    how many SUMT iterations to perform max
   ##
   penalty <- function(theta) {
      p <- A %*% theta + B
      sum(p*p)
   }
   ## Penalty gradient and Hessian are used only if corresponding function
   ## for the likelihood function is provided
   gPenalty <- function(theta) {
      2*(t(theta) %*% t(A) %*% A - t(B) %*% A)
   }
   hessPenalty <- function(theta) {
      2*t(A) %*% A
   }

   ## strip possible arguments of maxRoutine and call the function thereafter
   callWithoutMaxArgs <- function(theta, fName, ...) {
      return( callWithoutArgs( theta, fName = fName,
         args = names(formals(maxRoutine)), ... ) )
   }

   SUMTMessage <- function(code) {
      message <- switch(code,
                        "1" = "penalty close to zero",
                        "2" = "successive function values within tolerance limit",
                        "4" = "Outer iteration limit exceeded (increase SUMTMaxIter ?).",
                        paste("Code", code))
      return(message)
   }
   ## the penalized objective function
   Phi <- function(theta, ...) {
      llVal <- callWithoutMaxArgs( theta, "logLikFunc", fnOrig = fn,
         gradOrig = grad, hessOrig = hess, ... )
      llVal <- llVal - rho * penalty( theta ) / length( llVal )
      g <- attributes( llVal )$gradient
      if( !is.null( g ) ) {
         if( is.matrix( g ) ) {
            g <- g - matrix(
               rep( rho * gPenalty( theta ) / nrow( g ), each = nrow( g ) ),
               nrow = nrow( g ), ncol = ncol( g ) )
         } else {
            g <- g - rho * gPenalty( theta )
         }
         attributes( llVal )$gradient <- g
      }
      h <- attributes( llVal )$hessian
      if( !is.null( h ) ) {
         attributes( llVal )$hessian <- h - rho * hessPenalty( theta )
      }
      return( llVal )
   }
   ## gradient of the penalized objective function
   if(!is.null(grad)) {
      gradPhi<- function(theta, ...) {
         g <- grad(theta, ...)
         if(is.matrix(g)) {
            g <- g - matrix(
               rep( rho * gPenalty( theta ) / nrow( g ), each = nrow( g ) ),
               nrow = nrow( g ), ncol = ncol( g ) )
         } else {
            g <- g - rho * gPenalty( theta )
         }
         return( g )
      }
   } else {
      gradPhi <- NULL
   }
   ## Hessian of the penalized objective function
   if(!is.null(hess)) {
      hessPhi <- function(theta, ...) {
         return( hess(theta, ...) - rho*hessPenalty(theta) )
      }
   } else {
      hessPhi <- NULL
   }

   ##
   A <- constraints$eqA
   B <- constraints$eqB
    ## <NOTE>
    ## For the penalized minimization, the Newton-type nlm() may be
    ## computationally infeasible (although it works much faster for
    ## fitting ultrametrics to the Phonemes data).
    ## De Soete recommends using Conjugate Gradients.
    ## We provide a simple choice: by default, optim(method = "CG") is
    ## used.  If method is non-null and not "nlm", we use optim() with
    ## this method.  In both cases, control gives the control parameters
    ## for optim().
    ## If method is "nlm", nlm() is used, in which case control is
    ## ignored.  Note that we call nlm() with checking analyticals
    ## turned off, as in some cases (e.g. when fitting ultrametrics) the
    ## penalty function is not even continuous ...
   ## 
    ## Note also that currently we do not check whether optimization was
    ## "successful" ...
    ## </NOTE>
   ##
   rho <- 0
   result <- maxRoutine(fn=Phi, grad=gradPhi, hess=hessPhi,
                   start=start,
                   print.level=max(print.level - 1, 0),
                   ...)
   theta <- coef(result)
   if(print.level > 0) {
      cat("SUMT initial: rho = ", rho,
          ", function = ", callWithoutMaxArgs( theta, "logLikFunc",
            fnOrig = fn, gradOrig = grad, hessOrig = hess, ... ),
          ", penalty = ", penalty(theta), "\n")
      cat("Estimate:")
      print(theta)
   }
   ## <TODO>
   ## Better upper/lower bounds for rho?
   if(is.null(SUMTRho0)) {
       rho <- max( callWithoutMaxArgs( theta, "logLikFunc", fnOrig = fn,
         gradOrig = grad, hessOrig = hess, ... ), 1e-3) /
         max(penalty(start), 1e-3)
   } else {
       rho <- SUMTRho0
   }
   ## </TODO>
   iter <- 1L
   repeat {
      thetaOld <- theta
      result <- maxRoutine(fn=Phi, grad=gradPhi, hess=hessPhi,
                      start=thetaOld,
                      print.level=max(print.level - 1, 0),
                      ...)
      theta <- coef(result)
      if(print.level > 0) {
         cat("SUMT iteration ", iter,
             ": rho = ", rho, ", function = ", callWithoutMaxArgs( theta,
             "logLikFunc", fnOrig = fn, gradOrig = grad, hessOrig = hess, ... ),
             ", barrier value = ", penalty(theta), "\n", sep="")
         cat("Estimate:")
         print(theta)
      }
      if(max(abs(thetaOld - theta)) < SUMTTol) {
         SUMTCode <- 2
         break
      }
      if(penalty(theta) < SUMTPenaltyTol) {
         SUMTCode <- 1
         break
      }
      if(iter >= SUMTMaxIter) {
         SUMTCode <- 4
         break
      }
      iter <- iter + 1L
      rho <- SUMTQ * rho
   }
   ## Now we replace the resulting gradient and Hessian with those,
   ## calculated on the original function
   llVal <- callWithoutMaxArgs( theta, "logLikFunc", fnOrig = fn,
      gradOrig = grad, hessOrig = hess, sumObs = FALSE, ... )
   gradient <- attr( llVal, "gradient" )
   if( is.null( gradient ) ) {
      gradient <- callWithoutMaxArgs( theta, "logLikGrad", fnOrig = fn,
         gradOrig = grad, hessOrig = hess, sumObs = FALSE, ... )
   }
   if( !is.null( dim( gradient ) ) ) {
      if( nrow( gradient ) > 1 ) {
         gradientObs <- gradient
      }
      gradient <- colSums( gradient )
   } else if( length( start ) == 1 && length( gradient ) > 1 ) {
      gradientObs <- matrix( gradient, ncol = 1 )
      gradient <- sum( gradient )
   }
   result$gradient <- gradient
   names( result$gradient ) <- names( result$estimate )

   result$hessian <- callWithoutMaxArgs( theta, "logLikHess", fnOrig = fn,
      gradOrig = grad, hessOrig = hess, ... )
   result$constraints <- list(type="SUMT",
                             barrier.value=penalty(theta),
                              code=SUMTCode,
                              message=SUMTMessage(SUMTCode),
                             outer.iterations=iter
                             )
   if( exists( "gradientObs" ) ) {
      result$gradientObs <- gradientObs
      colnames( result$gradientObs ) <- names( result$estimate )
   }

   if( result$constraints$barrier.value > 0.001 ) {
      warning( "problem in imposing equality constraints: the constraints",
         " are not satisfied (barrier value = ",
         result$constraints$barrier.value, "). Try setting 'SUMTTol' to 0" )
   }
   return(result)
}
