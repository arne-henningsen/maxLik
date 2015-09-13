
## Function overwrite parameters of an existing MaxControl object using
## parameters supplied in a single list.
## We do not make it to a method: the signature would be indidistinguishable
## from add(maxControl, ...) where ... is a single list
addControlList <- function(x, y) {
   ## add list y to the control
   if(!inherits(x, "MaxControl")) {
      stop("'x' must be of class 'MaxControl'")
   }
   if(is.null(y)) {
      return(x)
   }
   if(!inherits(y, "list")) {
      stop("Control arguments to 'maxControl' must be supplied in the form of a list")
   }
   y <- y[names(y) %in% openParam(x)]
   if("iterlim" %in% names(y)) {
      slot(x, "iterlim") <- as.integer(y$iterlim)
      y$iterlim <- NULL
   }
   if("print.level" %in% names(y)) {
      slot(x, "printLevel") <- as.integer(y$print.level)
      y$print.level <- NULL
   }
   if("printLevel" %in% names(y)) {
      slot(x, "printLevel") <- as.integer(y$printLevel)
      y$printLevel <- NULL
   }
   for(p in names(y)) {
      slot(x, p) <- y[[p]]
   }
   validObject(x)
   return(x)
}
