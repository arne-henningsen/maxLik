
### shoud move checkMaxControl to a separate file but how to do it?

checkMaxControl <- function(object) {
   ## check validity of MaxControl objects
   if(!inherits(object, "MaxControl")) {
      stop("'MaxControl' object required.  Currently '",
           class(object), "'")
   }
   ##
   errors <- character(0)
   ##
   for(s in slotNames(object)) {
      if(length(slot(object, s)) != 1) {
         errors <- c(errors,
                     paste("'", s, "' must be of length 1, not ",
                           length(slot(object, s)), sep=""))
      }
   }
   ##
   if(slot(object, "tol") < 0) {
      errors <- c(errors, paste("'tol' must be non-negative, not ",
                                slot(object, "tol"), sep=""))
   }
   if(slot(object, "reltol") < 0) {
      errors <- c(errors, paste("'reltol' must be non-negative, not",
                                slot(object, "reltol")))
   }
   if(slot(object, "gradtol") < 0) {
      errors <- c(errors, paste("'gradtol' must be non-negative, not",
                                slot(object, "gradtol")))
   }
   if(slot(object, "steptol") < 0) {
      errors <- c(errors, paste("'steptol' must be non-negative, not",
                                slot(object, "steptol")))
   }
   if(slot(object, "lambdatol") < 0) {
      errors <- c(errors, paste("'lambdatol' must be non-negative, not",
                                slot(object, "lambdatol")))
   }
   if(!pmatch(slot(object, "qac"), c("stephalving", "marquardt"))) {
      errors <- c(errors, paste("'qac' must be 'stephalving' or 'marquadt', not",
                                slot(object, "qac")))
   }
   if(slot(object, "qrtol") < 0) {
      errors <- c(errors, paste("'qrtol' must be non-negative, not",
                                slot(object, "qrtol")))
   }
   if(slot(object, "lambda0") < 0) {
      errors <- c(errors, paste("'lambda0' must be non-negative, not",
                                slot(object, "lambda0")))
   }
   if(slot(object, "lambdaStep") <= 1) {
      errors <- c(errors, paste("'lambdaStep' must be > 1, not",
                                slot(object, "lambdaStep")))
   }
   if(slot(object, "maxLambda") < 0) {
      errors <- c(errors, paste("'maxLambda' must be non-negative, not",
                                slot(object, "maxLambda")))
   }
   if(slot(object, "iterlim") < 0) {
      errors <- c(errors, paste("'iterlim' must be non-negative, not",
                                slot(object, "iterlim")))
   }
   if(length(errors) > 0)
      return(errors)
   return(TRUE)
}

### MaxControls contains all control parameters for max* family
setClass("MaxControl",
         slots=representation(
             tol="numeric",
             reltol="numeric",
             gradtol="numeric",
             steptol="numeric",
                           #
             lambdatol="numeric",
             ## Qadratic Approximation Control
             qac="character",
             qrtol="numeric",
             lambda0="numeric",
             lambdaStep="numeric",
             maxLambda="numeric",
         ## Optim Nelder-Mead:
         alpha="numeric",
         beta="numeric",
         gamma="numeric",
         ##
             iterlim="integer",
             ##
             printLevel="integer"),
         prototype=prototype(
             tol=1e-8,
             reltol=sqrt(.Machine$double.eps),
             gradtol=1e-6,
             steptol=1e-10,
                           #
             lambdatol=1e-6,
                           #
             qac="stephalving",
             qrtol=1e-10,
             lambda0=1e-2,
             lambdaStep=2,
             maxLambda=1e12,
         ## Optim Nelder-Mead
         alpha=1,
         beta=0.5,
         gamma=2,
                           #
             iterlim=150L,
             printLevel=0L),
         validity=checkMaxControl)
