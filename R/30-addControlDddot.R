
## Method to overwrite parameters of an existing MaxControl object
addControlDddot <- function(x, ...) {
   ## add ... to the control
   ## only fetch openParam from ...
   dddot <- list(...)
   dddot <- dddot[names(dddot) %in% openParam(x)]
   if("iterlim" %in% names(dddot)) {
      slot(x, "iterlim") <- as.integer(dddot$iterlim)
      dddot$iterlim <- NULL
   }
   if("print.level" %in% names(dddot)) {
      slot(x, "printLevel") <- as.integer(dddot$print.level)
      dddot$print.level <- NULL
   }
   if("printLevel" %in% names(dddot)) {
      slot(x, "printLevel") <- as.integer(dddot$printLevel)
      dddot$printLevel <- NULL
   }
   for(p in names(dddot)) {
      slot(x, p) <- dddot[[p]]
   }
   validObject(x)
   return(x)
}

setMethod("maxControl", "MaxControl", addControlDddot)
