print.summary.maxLik <- function( x, ... ) {
   cat("--------------------------------------------\n")
   cat("Maximum Likelihood estimation\n")
   cat(maximType(x), ", ", nIter(x), " iterations\n", sep="")
   cat("Return code ", returnCode(x), ": ", returnMessage(x), "\n", sep="")
   if(!is.null(x$estimate)) {
      cat("Log-Likelihood:", x$loglik, "\n")
      cat(x$NActivePar, " free parameters\n")
      cat("Estimates:\n")
      printCoefmat(x$estimate)
   }
   if(!is.null(x$constraints)) {
      cat("\nWarning: constrained likelihood estimation.",
          "Inference is probably wrong\n")
      cat("Constrained optimization based on", x$constraints$type,
          "\n")
      cat(x$constraints$outer.iterations,
          " outer iterations, barrier value",
          x$constraints$barrier.value, "\n")
   }
   cat("--------------------------------------------\n")
}

summary.maxLik <- function(object, eigentol=1e-12,... ) {
   ## object      object of class "maxLik"
   ## 
   ## RESULTS:
   ## list of class "summary.maxLik" with following components:
   ## maximum    : function value at optimum
   ## estimate   : estimated parameter values at optimum
   ## gradient   :           gradient at optimum
   ## code       : code of convergence
   ## message    : message, description of the code
   ## iterations : number of iterations
   ## type       : type of optimisation
   ##
   result <- object$maxim
   nParam <- length(coef <- coef.maxLik(object))
   if(!is.null(object$activePar)) {
      activePar <- object$activePar
   } else {
      activePar <- rep(TRUE, nParam)
   }
   if(object$code < 100) {
     hess <- hessian(object)[activePar, activePar] 
     hessev <- abs(eigen(hess, symmetric=TRUE, only.values=TRUE)$values)
     if(min(hessev) > (eigentol*max(hessev))){      
#    if(min(abs(eigen(hessian(object)[activePar,activePar],
#             symmetric=TRUE, only.values=TRUE)$values)) > eigentol) {
       varcovar <- matrix(0, nParam, nParam)
       varcovar[activePar,activePar] <-
         solve(-hessian(object)[activePar,activePar])
       hdiag <- diag(varcovar)
       if(any(hdiag < 0)) {
         warning("Diagonal of variance-covariance matrix not positive!\n")
       }
       stdd <- sqrt(hdiag)
       t <- coef/stdd
       p <- 2*pnorm( -abs( t))
     } else {
       stdd <- 0
       t <- 0
       p <- 0
     }
     results <- cbind("Estimate"=coef, "Std. error"=stdd, "t value"=t, "Pr(> t)"=p)
     Hess <- NULL
   } else {
     results <- NULL
     Hess <- NULL
   }
   summary <- list(maximType=object$type,
                   iterations=object$iterations,
                   returnCode=object$code,
                   returnMessage=object$message,
                   loglik=object$maximum,
                   estimate=results,
                   hessian=Hess,
                   activePar=object$activePar,
                   NActivePar=sum(object$activePar),
                   constraints=object$constraints)
   class(summary) <- "summary.maxLik"
   summary
}
