### SUMT function borrowed from 'clue' package
### 
### Adapted for linear constraints
sumt <- function(fn, grad=NULL, hess=NULL,
                 start,
                 maxRoutine, constraints, 
                 SUMTTol = sqrt(.Machine$double.eps),
                 SUMTQ = 10,
                 SUMTRho0 = max(fn(start), 1e-4)/max(penalty(start), 1e-4),
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
   Phi <- function(theta)
       fn(theta) - rho * penalty(theta)
   if(!is.null(grad)) {
      gradPhi<- function(theta) 
         grad(theta) - rho*gPenalty(theta)
   }
   else
       gradPhi <- NULL
   if(!is.null(hess)) {
      hessPhi <- function(theta) 
         hess(theta) - rho*hessPenalty(theta)
   }
   else
       hessPhi <- NULL
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
   a <- maxRoutine(fn=Phi, grad=gradPhi, hess=hessPhi,
                   start=start,
                   print.level=print.level - 1,
                   ...)
   theta <- coef(a)
   if(print.level > 0) {
      cat("Initial: rho = ", rho, ", function = ", fn(theta),
          ", penalty = ", penalty(theta), "\n")
      cat("Estimate:")
      print(theta)
   }
   ## <TODO>
   ## Better upper/lower bounds for rho?
   rho <- SUMTRho0
   ## </TODO>
   iter <- 1L
   repeat {
      ## <TODO>
      ## Shouldnt't we also have maxiter, just in case ...?
      ## </TODO>
      thetaOld <- theta
      a <- maxRoutine(fn=Phi, grad=gradPhi, hess=hessPhi,
                      start=start,
                      print.level=print.level - 1,
                      ...)
      theta <- coef(a)
      if(print.level > 0) {
         cat("SUMT iteration ", iter,
             ": rho = ", rho, ", function = ", fn(theta),
             ", penalty = ", penalty(theta), "\n", sep="")
         cat("Estimate:")
         print(theta)
      }
      if(max(abs(thetaOld - theta)) < SUMTTol)
          break
      if(iter >= SUMTMaxIter)
          break
      iter <- iter + 1L
      rho <- SUMTQ * rho
   }
   return(a)
}
