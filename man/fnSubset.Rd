\name{fnSubset}
\alias{fnSubset}
\title{
  Call fnFull with variable and fixed parameters
}
\description{
  Combine variable parameters in \code{x} with \code{xFixed} and pass to
  \code{fnFull}.  Useful for optimizing over a subset of parameters
  without writing a separate function.  Values are combined by name if
  available.  Otherwise, \code{xFull} is constructed by position (the
  default).  
}
\usage{
fnSubset(x, fnFull, xFixed, xFull=c(x, xFixed), ...) 
}
\arguments{
  \item{x}{
    Variable parameters to be passed to \code{fnFull}.  
  }
  \item{fnFull}{
    Function whose first argument has length = length(xFull).  
  }
  \item{xFixed}{
    Parameters to be combined with \code{x} to construct the first
    argument for a call to \code{fnFull}.  
  }
  \item{xFull}{
    Prototype initial argument for \code{fnFull}.  
  }
  \item{\dots}{
    Optional arguments passed to \code{fnFull}.  
  }
}
\details{
  1.  Confirm that length(x) + length(xFixed) = length(xFull)

  2.  If \code{xFull} has names, match at least \code{xFixed} by name.
  Else xFull = c(x, xFixes), the default.  

  3.  fnFull(xFull, ...) 
}  
\value{
  value returned by \code{fnFull} 
}
%\references{ }
\author{
  Spencer Graves 
}
\seealso{
  \code{\link{optim}}
  \code{\link[dlm]{dlmMLE}}
  \code{\link{maxLik}}
  \code{\link{maxNR}} 
}
\examples{
##
## Test with 'optim'
##
fn <- function(x, x0)(x[2]-2*x[1]-x0)^2
fullEst <- optim(1:2, fn, x0=3)

# Fix the last component 
est4 <- optim(1, fnSubset, x0=3, fnFull=fn, xFixed=4)

# Fix the first component
fnSubset(1, fn, c(a=4), c(a=1, b=2), x0=3)
# After substitution:  xFull = c(a=4, b=1),
# so fn = (1-2*4-3)^2 = (-10)^2 = 100 

est4. <- optim(1, fnSubset, x0=3, fnFull=fn, xFixed=c(a=4), 
               xFull=c(a=1, b=2))
# At optim: xFull=c(a=4, b=10.9),
# so fn = (10.9-2*4-3)^2 = (-0.1)^2 = 0.01 

##
## Test with maxNR
##
#  fn2max = -fn
fn2max <- function(x, x0, ...)(-(x[2]-2*x[1]-x0)^2)
#  Need "..." here when called directly from maxNR,
#  because maxNR will also pass 'constantPar'  
# Fix the last component 
NR4 <- maxNR(fnSubset, start=1, x0=3, fnFull=fn2max, xFixed=4)
# Same thing using maxNR(..., activePar)
NR4. <- maxNR(fn2max, start=c(1, 4), x0=3, constantPar=2)

##
## Test with maxLik
##
# Same as maxNR
max4 <- maxLik(fnSubset, start=1, x0=3, fnFull=fn2max, xFixed=4)
# Same thing using constantPar in maxNR, called by maxLik 
max4 <- maxLik(fn2max, start=c(1, 4), x0=3, constantPar=2)

}
\keyword{optimize}% at least one, from doc/KEYWORDS
