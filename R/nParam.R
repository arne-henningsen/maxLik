## Return the #of parameters of model
nParam.maxim <- function(x, free=FALSE, ...) {
   if(free)
       sum(x$activePar)
   else
       length( x$estimate )
}
