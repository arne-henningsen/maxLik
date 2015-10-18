
## Function overwrite parameters of an existing MaxControl object using
## parameters supplied in a single list.
## We do not make it to a method: the signature would be indistinguishable
## from add(maxControl, ...) where ... is a single list
addControlList <- function(x, y) {
   ## add list y to the control
   setSlot <- function(openName,
                       slotName=openName,
                       convert=function(x) x
                       ) {
      if(!is.null(y[[openName]])) {
         i <- tail(which(names(y) %in% openName), 1)
                           # pick the last occurrence: allow user to overwrite defaults
         slot(x, slotName) <- convert(y[[i]])
      }
      assign("x", x, envir=parent.frame())
                           # save modified x into parent frame
   }
   if(!inherits(x, "MaxControl")) {
      stop("'x' must be of class 'MaxControl'")
   }
   if(is.null(y)) {
      return(x)
   }
   if(!inherits(y, "list")) {
      stop("Control arguments to 'maxControl' must be supplied in the form of a list")
   }
   knownNames <- union(openParam(x), slotNames(x))
   if(any(uNames <- !(names(y) %in% knownNames))) {
      cat("Unknown control options:\n")
      print(names(y)[uNames])
      stop("Unknown options not accepted")
   }
   ##
   setSlot("tol")
   setSlot("reltol")
   setSlot("gradtol")
   setSlot("lambdatol")
   setSlot("qrtol")
   ## QAC
   setSlot("QAC", "qac")
   setSlot("qac")
   setSlot("Marquardt_lambda0")
   setSlot("Marquardt_lambdaStep")
   setSlot("Marquardt_maxLambda")
   setSlot("marquardt_lambda0")
   setSlot("marquardt_lambdaStep")
   setSlot("marquardt_maxLambda")
   ## NM
   setSlot("NM_alpha", "nm_alpha")
   setSlot("NM_beta", "nm_beta")
   setSlot("NM_gamma", "nm_gamma")
   setSlot("nm_alpha", "nm_alpha")
   setSlot("nm_beta", "nm_beta")
   setSlot("nm_gamma", "nm_gamma")
   setSlot("alpha", "nm_alpha")
   setSlot("beta", "nm_beta")
   setSlot("gamma", "nm_gamma")
   ## SANN
   setSlot("SANN_cand", "sann_cand")
   setSlot("SANN_temp", "sann_temp")
   setSlot("SANN_tmax", "sann_tmax", as.integer)
   setSlot("SANN_randomSeed", "sann_randomSeed", as.integer)
   setSlot("sann_cand")
   setSlot("sann_temp")
   setSlot("sann_tmax", convert=as.integer)
   setSlot("sann_randomSeed", "sann_randomSeed", as.integer)
   setSlot("cand", "sann_cand")
   setSlot("temp", "sann_temp")
   setSlot("tmax", "sann_tmax", as.integer)
   setSlot("random.seed", "sann_randomSeed", as.integer)
   ##
   setSlot("iterlim", convert=as.integer)
   setSlot("print.level", "printLevel", as.integer)
   setSlot("printLevel", convert=as.integer)
   ##
   validObject(x)
   return(x)
}

### Method for 'MaxControl' objects: add the second argument, list
setMethod("maxControl", signature("MaxControl"), addControlList)
