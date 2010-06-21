maxNM <- function(fn, grad=NULL, hess=NULL,
                  start, fixed = NULL,
                  print.level=0,
                  iterlim=500,
                  constraints=NULL,
                  tol=1e-8, reltol=tol,
                  finalHessian=TRUE,
                  parscale=rep(1, length=length(start)),
                  alpha=1, beta=0.5, gamma=2,
                  ...) {
   ## Wrapper of optim-based 'Nelder-Mead' optimization
   ## 
   ## contraints    constraints to be passed to 'constrOptim'
   ## hessian:   how (and if) to calculate the final Hessian:
   ##            FALSE   not calculate
   ##            TRUE    use analytic/numeric Hessian
   ##            bhhh/BHHH  use information equality approach
   ## ... :      further arguments to fn()
   ##
   ## Note: grad and hess are for compatibility only, SANN uses only fn values

   result <- maxOptim( fn = fn, grad = grad, hess = hess,
      start = start, method = "Nelder-Mead", fixed = fixed,
      print.level = print.level, iterlim = iterlim, constraints = constraints,
      tol = tol, reltol = reltol,
                      finalHessian=finalHessian,
                      parscale = parscale,
      alpha = alpha, beta = beta, gamma = gamma,
      ... )

   return(result)
}
