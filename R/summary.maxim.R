print.summary.maxim <- function( x,
                                max.rows=getOption("max.rows", 20),
                                max.cols=getOption("max.cols", 7),
                                ... ) {
   summary <- x
   cat("--------------------------------------------\n")
   cat(.maximType(summary), "\n")
   cat("Number of iterations:", .nIter(summary), "\n")
   cat("Return code:", .returnCode(summary), "\n")
   cat(.returnMessage(summary), "\n")
   unsuccStep <- .unsuccStep( summary )
   if(!is.null(unsuccStep)) {
      cat("Last (unsuccessful) step: function value", unsuccStep$value,
         "\n")
      print(unsuccStep$parameters)
   }
   estimate <- .estimate( summary )
   if(!is.null(estimate)) {
      cat("Function value:", .maximum( summary ), "\n")
      cat("Estimates:\n")
      printRowColLimits(estimate, max.rows, max.cols, ...)
      hess <- .hessian( summary )
      if(!is.null(hess)) {
         cat("Hessian:\n")
         printRowColLimits(hess, max.rows, max.cols, ...)
      }
   }
   constraints <- .constraints( summary )
   if(!is.null(constraints)) {
      cat("\nConstrained optimization based on", constraints$type,
          "\n")
      if(!is.null(constraints$code))
         cat("Return code:", constraints$code, "\n")
                           # note: this is missing for 'constrOptim'
      if(!is.null(constraints$message))
         cat(constraints$message, "\n")
                           # note: this is missing for 'constrOptim'
      cat(constraints$outer.iterations,
          " outer iterations, barrier value",
          constraints$barrier.value, "\n")
   }
   cat("--------------------------------------------\n")
}

summary.maxim <- function(object, hessian=FALSE, unsucc.step=FALSE,
   ... ) {
   ## The object of class "maxim" should include following components:
   ## maximum    : function value at optimum
   ## estimate   : matrix, estimated parameter values and gradient at optimum
   ## hessian    :           hessian
   ## code       : code of convergence
   ## message    : message, description of the code
   ## last.step  : information about last step, if unsuccessful
   ## iterations : number of iterations
   ## type       : type of optimisation
   ##
   code <- .returnCode( object )
   if(code == 3 & unsucc.step) {
      lastStep <- .lastStep( object )
      a <- cbind(lastStep$theta0, lastStep$theta1)
      dimnames(a) <- list(parameter=.parameterNames(object),
                        c("current par", "new par"))
      unsucc.step <- list(value=lastStep$f0,
                        parameters=a)
   } else {
      unsucc.step <- NULL
   }
   estimate <- cbind("estimate"=.estimate( object ),
                     "gradient"=.gradient( object ))
   if(hessian) {
      H <- .hessian( object )
   }
   else {
      H <- NULL
   }
   summary <- list(maximum=.maximum( object ),
                   type=.maximType( object ),
                  iterations=.nIter( object ),
                  code=code,
                  message=.returnMessage( object ),
                   unsucc.step=unsucc.step,
                  estimate=estimate,
                   hessian=H,
                   constraints=.constraints( object ))
   class(summary) <- c("summary.maxim", class(summary))
   summary
}
