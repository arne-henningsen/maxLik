addDddotToControl <- function(control, ...) {
   ## add ... to the control
   ## only fetch openParam from ...
   dddot <- list(...)
   dddot <- dddot[names(dddot) %in% openParam(control)]
   if("iterlim" %in% names(dddot)) {
      slot(control, "iterlim") <- as.integer(dddot$iterlim)
      dddot$iterlim <- NULL
   }
   if(any(c("printLevel", "print.level") %in% names(dddot))) {
      slot(control, "printLevel") <- as.integer(dddot$iterlim)
      dddot$printLevel <- dddot$print.level <- NULL
   }
   for(p in names(dddot)) {
      slot(control, p) <- dddot[[p]]
   }
   validObject(control)
   return(control)
}
