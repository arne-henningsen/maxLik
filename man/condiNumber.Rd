\name{condiNumber}
\alias{condiNumber}
\alias{condiNumber.default}
\alias{condiNumber.maxLik}
\title{Print matrix condition numbers column-by-column}
\description{
  This function prints the condition number of a matrix while adding
  columns one-by-one.  This is useful for testing multicollinearity and
  other numerical problems.  This is a generic function with default
  method, and a method for \code{maxLik} objects.
}
\usage{
condiNumber(x, ...)
\method{condiNumber}{default}(x, exact = FALSE, norm = FALSE,
   print.level=1, digits = getOption( "digits" ), ... )
\method{condiNumber}{maxLik}(x, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{numeric matrix, condition numbers of which are to be printed}
  \item{exact}{logical, should condition numbers be exact or
    approximations (see \code{link{kappa}})}
  \item{norm}{logical, whether the columns should be normalised to have
    unit norm}
  \item{print.level}{numeric, positive value will output the numbers
    during the calculations.  Useful for interactive work.}
  \item{digits}{minimal number of significant digits to print
      (only relevant if argument \code{print.level} is larger than zero).}
  \item{\dots}{Further arguments to \code{condiNumber.default}
      are currently ignored;
      further arguments to \code{condiNumber.maxLik}
      are passed to \code{condiNumber.default}.} 
}
\details{
  Statistical model often fail because of strong correlation between
  the explanatory variables in the linear index (multicollinearity)
  or because the evaluated maximum of a non-linear model is
  virtually flat.  In both cases, the (near) singularity of the related
  matrices may give us a hint how to improve the results.

  \code{condiNumber} allows to inspect the matrices column-by-column and
  see which variables leads to a jump in the condition
  number (cause singularity).
  If the matrix column name does not immediately indicate the
  problem, one may run an OLS model by estimating this column
  using all the previous columns as explanatory variables.  Those
  columns that explain almost all the variation in the current one will
  have very high
  \eqn{t}{t}-values. 
}
\value{
  Invisible vector of condition numbers by column.  If the start values
  for \code{\link{maxLik}} are named, the condition numbers are named
  accordingly. 
}
\references{W. Greene, Advanced Econometrics, p ... }
\author{Ott Toomet \email{otoomet@ut.ee}}
\seealso{\code{\link{kappa}}}
\examples{
   set.seed(0)
   ## generate a simple multicollinear dataset
   x1 <- runif(100)
   x2 <- runif(100)
   x3 <- x1 + x2 + 0.000001*runif(100) # this is virtually equal to x1 + x2
   x4 <- runif(100)
   y <- x1 + x2 + x3 + x4 + rnorm(100)
   m <- lm(y ~ -1 + x1 + x2 + x3 + x4)
   print(summary(m)) # note the low t-values while R^2 is 0.88.
                     # This hints multicollinearity
   condiNumber(model.matrix(m)) # this _prints_ condition numbers.
                                # note the values 'explode' with x3
   ## we may test the results further:
   print(summary(lm(x3 ~ -1 + x1 + x2))) # Note the high t-values and R^2
}
\keyword{math}
\keyword{utilities}
\keyword{debugging} % is it debugging?
