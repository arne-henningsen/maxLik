## Return the #of parameters of model
nParam.maxim <- function(x, free=FALSE, ...) {
   if(free)
       sum( activePar( x ) )
   else
       length( x$estimate )
}
