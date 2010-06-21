maxBFGS <- function(fn, grad=NULL, hess=NULL,
                    start, fixed = NULL,
                    print.level=0,
                    iterlim=200,
                    constraints=NULL,
                    tol=1e-8, reltol=tol,
                    finalHessian=TRUE,
                    parscale=rep(1, length=length(start)),
                    ## sumt parameters
                    ...) {
   ## Wrapper of optim-based 'BFGS' optimization
   ## 
   ## contraints    constraints to be passed to 'constrOptim'
   ## finalHessian:   how (and if) to calculate the final Hessian:
   ##            FALSE   not calculate
   ##            TRUE    use analytic/numeric Hessian
   ##            bhhh/BHHH  use information equality approach
   ##
   ## ...           further arguments to fn() and grad()

   result <- maxOptim( fn = fn, grad = grad, hess = hess,
      start = start, method = "BFGS", fixed = fixed,
      print.level = print.level, iterlim = iterlim, constraints = constraints,
      tol = tol, reltol = reltol,
                      finalHessian=finalHessian,
                      parscale = parscale,
      ... )

   return(result)
}
