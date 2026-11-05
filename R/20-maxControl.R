
### Default constructor of MaxControl object:
### take a list of parameters and overwrite the default values
maxControl.default <- function(...) {
   result <- new("MaxControl")
   result <- addControlDddot(result, ...)
   return(result)
}

### Standard method for any arguments
setGeneric("maxControl",
           function(x, ...) standardGeneric("maxControl")
           )

### Method for 'maxim' objects: fetch the stored MaxControl
### '.field()' is used rather than '$': a user-defined '$' method for
### 'maxLik' or 'maxim' objects would otherwise change how maxLik reads
### the components of its own objects.
setMethod("maxControl", "maxim", function(x, ...) .field(x, "control"))

### Method for missing arguments: just default values
setMethod("maxControl", "missing", maxControl.default)
