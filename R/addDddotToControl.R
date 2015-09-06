addDddotToControl <- function(control, ...) {
   ## add ... to the control
   ## only fetch openParam from ...
   dddot <- list(...)
   dddot <- dddot[names(dddot) %in% openParam(control)]
   if("iterlim" %in% names(dddot)) {
      slot(control, "iterlim") <- as.integer(dddot$iterlim)
      dddot$iterlim <- NULL
   }
   if("print.level" %in% names(dddot)) {
      slot(control, "printLevel") <- as.integer(dddot$print.level)
      dddot$print.level <- NULL
   }
   if("printLevel" %in% names(dddot)) {
      slot(control, "printLevel") <- as.integer(dddot$printLevel)
      dddot$printLevel <- NULL
   }
   for(p in names(dddot)) {
      slot(control, p) <- dddot[[p]]
   }
   validObject(control)
   return(control)
}
