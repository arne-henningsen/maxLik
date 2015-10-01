
## Method to overwrite parameters of an existing MaxControl object
addControlDddot <- function(x, ...) {
   ## add ... to the control
   ## only fetch openParam from ...
   setSlot <- function(openName,
                       slotName=openName,
                       convert=function(x) x
                       ) {
      if(!is.null(dddot[[openName]])) {
         slot(x, slotName) <- convert(dddot[[openName]])
      }
      assign("x", x, envir=parent.frame())
                           # save modified x into parent frame
   }
   dddot <- list(...)
   dddot <- dddot[names(dddot) %in% openParam(x)]
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
   setSlot("alpha")
   setSlot("beta")
   setSlot("gamma")
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

setMethod("maxControl", "MaxControl", addControlDddot)
