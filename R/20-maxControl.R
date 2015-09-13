
### Default constructor of MaxControl object:
### take a list of parameters and overwrite the default values
maxControl.default <- function(...) {
   ## extract MaxControl arguments and add these to the control list
   ## 
   result <- new("MaxControl")
   dddot <- list(...)
   okNames <- c(slotNames("MaxControl"), "print.level")
   if(length(wrong <- names(dddot)[!(names(dddot) %in% okNames)]) > 0) {
      stop(wrong, " is not a supported argument by 'maxControl'")
   }
   if("print.level" %in% names(dddot)) {
      names(dddot)[names(dddot) == "print.level"] <- "printLevel"
   }
   for(s in names(dddot)) {
      val <- dddot[[s]]
      if(s == "printLevel") {
         val <- as.integer(val)
      }
      if(s == "iterlim") {
         val <- as.integer(val)
      }
      slot(result, s) <- val
   }
   validObject(result)
   return(result)
}

### Standard method for any arguments
setGeneric("maxControl",
           function(x, ...) standardGeneric("maxControl")
           )

### Method for missing arguments: just default values
setMethod("maxControl", "missing", maxControl.default)
