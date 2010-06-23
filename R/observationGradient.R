
### The function tests whether a given gradient is given observation-wise.  It tests essentially the # of row in the gradient
observationGradient <- function(g) {
   if(is.null(dim(g)))
       return(FALSE)
   if(nrow(g) == 1)
       return(FALSE)
   return(TRUE)
}
