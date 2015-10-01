
### Default constructor of MaxControl object:
### take a list of parameters and overwrite the default values
maxControl.default <- function(...) {
   ## extract MaxControl arguments and add these to the control list
   ## 
   result <- new("MaxControl")
   dddot <- list(...)
   okNames <- c(slotNames("MaxControl"), "print.level",
                "temp", "tmax", "cand", "random.seed")
   if(length(wrong <- names(dddot)[!(names(dddot) %in% okNames)]) > 0) {
      stop(wrong, " is not a supported argument by 'maxControl'")
   }
   if("print.level" %in% names(dddot)) {
      names(dddot)[names(dddot) == "print.level"] <- "printLevel"
   }
   if("temp" %in% names(dddot)) {
      names(dddot)[names(dddot) == "temp"] <- "SANN_temp"
   }
   if("tmax" %in% names(dddot)) {
      names(dddot)[names(dddot) == "tmax"] <- "SANN_tmax"
   }
   if("cand" %in% names(dddot)) {
      names(dddot)[names(dddot) == "cand"] <- "SANN_cand"
   }
   if("random.seed" %in% names(dddot)) {
      names(dddot)[names(dddot) == "random.seed"] <- "SANN_randomSeed"
   }
   for(s in names(dddot)) {
      val <- dddot[[s]]
      if(s == "printLevel") {
         val <- as.integer(val)
      }
      if(s == "iterlim") {
         val <- as.integer(val)
      }
      if(s == "SANN_randomSeed") {
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

### Method for 'maxim' objects: fetch the stored MaxControl
setMethod("maxControl", "maxim", function(x, ...) x$control)

### Method for missing arguments: just default values
setMethod("maxControl", "missing", maxControl.default)
