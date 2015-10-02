
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
         slot(x, slotName) <- convert(y[[openName]])
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
   setSlot("qac")
   setSlot("lambda0")
   setSlot("lambdaStep")
   setSlot("maxLambda")
   ## NM
   setSlot("alpha", "NM_alpha")
   setSlot("beta", "NM_beta")
   setSlot("gamma", "NM_gamma")
   ## SANN
   setSlot("cand", "SANN_cand")
   setSlot("temp", "SANN_temp")
   setSlot("tmax", "SANN_tmax", as.integer)
   setSlot("random.seed", "SANN_randomSeed", as.integer)
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
