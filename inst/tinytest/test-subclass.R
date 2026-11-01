### Test that maxLik methods do not dispatch to subclass methods internally

if(!requireNamespace("tinytest", quietly = TRUE)) {
   message("These tests require 'tinytest' package\n")
   q("no")
}
library(tinytest)  # for running it as a separate script
library(maxLik)

## Temporarily install S3 methods in the global environment, where they are
## found by generics called from the maxLik namespace.
## Any methods with the same name .GlobalEnv are restored afterwards.
with_s3_methods <- function(methods, code) {
   methodNames <- names(methods)
   ## Store the existing methods in .GlobalEnv into 'old'
   old <- lapply(methodNames, function(n) {
      if(exists(n, envir = .GlobalEnv, inherits = FALSE)) {
         get(n, envir = .GlobalEnv)
      } else {
         NULL
      }
   })
   names(old) <- methodNames
   ## Restore the methods in .GlobalEnv on exit
   on.exit({
      for(n in methodNames) {
         if(is.null(old[[n]])) {
            rm(list = n, envir = .GlobalEnv)
         } else {
            assign(n, old[[n]], envir = .GlobalEnv)
         }
      }
   }, add = TRUE)
   ## Assign the methods to .GlobalEnv
   for(n in methodNames) {
      assign(n, methods[[n]], envir = .GlobalEnv)
   }
   ## And compute with the updated .GlobalEnv
   force(code)
}


## Create coef.foo()-type functions that fail with the appropriate
## message
with_poison <- function(generics, code) {
   methods <- lapply(generics, function(g) {
      ## Store the generic name in a new environment in order to
      ## avoid being overwritten in a common environment later
      local({
         gg <- g
         function(...) stop(gg, ".foo was called", call. = FALSE)
      })
   })
   names(methods) <- paste0(generics, ".foo")
   with_s3_methods(methods, code)
}

## Run 'fun(object)' with each generic in 'generics' poisoned and compare
## the result with the same call on the plain object.
expect_no_dispatch <- function(generics, fun, object, reference, label) {
   for(g in generics) {
      result <- tryCatch(
         with_poison(g, fun(object)),
         error = function(e) e
      )
      expect_false(
         inherits(result, "error"),
         info = paste0(label, "() called ", g, ".foo()")
      )
      if(!inherits(result, "error")) {
         expect_equal(
            result, reference,
            info = paste0(label, "() unaffected by ", g, ".foo()")
         )
      }
   }
}

as_foo <- function(x) {
   class(x) <- c("foo", class(x))
   x
}

## ---------------------------------------------------------------- objects
set.seed(123)
y <- rnorm(20, mean = 1, sd = 2)
ll <- function(theta) {
   sum(dnorm(y, mean = theta[1], sd = theta[2], log = TRUE))
}
m <- maxLik(ll, start = c(1, 2))

## Observation-level likelihood and gradient, needed for estfun/bread/nobs
llInd <- function(theta) {
   dnorm(y, mean = theta[1], sd = theta[2], log = TRUE)
}
gradInd <- function(theta) {
   mu <- theta[1]
   sigma <- theta[2]
   cbind((y - mu) / sigma^2,
         -1 / sigma + (y - mu)^2 / sigma^3)
}
mObs <- maxLik(llInd, grad = gradInd, start = c(1, 2))

## A plain 'maxim' object
f <- function(theta) {
   -(theta[1] - 1)^2 - (theta[2] - 2)^2
}
mx <- maxNR(f, start = c(0, 0))

## Summary objects
s <- summary(m)
smx <- summary(mx)

## Subclassed versions
x <- as_foo(m)
xObs <- as_foo(mObs)
xmx <- as_foo(mx)
xs <- as_foo(s)
xsmx <- as_foo(smx)

## ------------------------------------------------------- summary.maxLik()
reference <- maxLik:::summary.maxLik(m)
expect_no_dispatch(
   c("coef", "activePar", "stdEr", "vcov", "hessian", "nParam",
     "maxValue", "maximType", "nIter", "returnCode", "returnMessage",
     "logLik", "AIC"),
   maxLik:::summary.maxLik, x, reference, "summary.maxLik"
)

## nParam() is only reached through .vcov() when the Hessian is absent
x_nohess <- x
x_nohess$hessian <- NULL
m_nohess <- m
m_nohess$hessian <- NULL
reference_nohess <- maxLik:::summary.maxLik(m_nohess)
expect_no_dispatch(
   "nParam",
   maxLik:::summary.maxLik, x_nohess, reference_nohess,
   "summary.maxLik without Hessian"
)

## The exported summary() path must also be unaffected
result <- tryCatch(
   with_poison(
      c("coef", "activePar", "stdEr", "vcov", "hessian", "nParam",
        "maxValue", "maximType", "nIter", "returnCode", "returnMessage"),
      summary(x)
   ),
   error = function(e) e
)
expect_false(
   inherits(result, "error"),
   info = "summary() must not use subclass accessor methods"
)
if(!inherits(result, "error")) {
   expect_equal(result, reference)
}

## ------------------------------------------------------------ print methods
expect_no_dispatch(
   c("coef", "activePar", "maxValue", "maximType", "nIter",
     "returnCode", "returnMessage", "vcov", "hessian", "nParam", "stdEr"),
   function(object) capture.output(maxLik:::print.maxLik(object)),
   x, capture.output(maxLik:::print.maxLik(m)), "print.maxLik"
)

expect_no_dispatch(
   c("coef", "activePar", "maxValue", "maximType", "nIter",
     "returnCode", "returnMessage", "vcov", "hessian", "nParam", "stdEr"),
   function(object) capture.output(maxLik:::print.maxim(object)),
   xmx, capture.output(maxLik:::print.maxim(mx)), "print.maxim"
)

expect_no_dispatch(
   c("maximType", "nIter", "returnCode", "returnMessage", "maxValue",
     "coef", "activePar", "vcov", "stdEr", "hessian", "nParam"),
   function(object) capture.output(maxLik:::print.summary.maxLik(object)),
   xs, capture.output(maxLik:::print.summary.maxLik(s)),
   "print.summary.maxLik"
)

expect_no_dispatch(
   c("maximType", "nIter", "returnCode", "returnMessage", "maxValue",
     "coef", "gradient", "hessian", "activePar"),
   function(object) capture.output(maxLik:::print.summary.maxim(object)),
   xsmx, capture.output(maxLik:::print.summary.maxim(smx)),
   "print.summary.maxim"
)

## ------------------------------------------------- logLik and AIC methods
expect_no_dispatch(
   c("maxValue", "activePar", "coef", "nParam", "vcov", "stdEr", "hessian"),
   maxLik:::logLik.maxLik, x, maxLik:::logLik.maxLik(m),
   "logLik.maxLik"
)

expect_no_dispatch(
   c("maxValue", "activePar", "coef", "logLik"),
   maxLik:::logLik.summary.maxLik, xs,
   maxLik:::logLik.summary.maxLik(s),
   "logLik.summary.maxLik"
)

expect_no_dispatch(
   c("logLik", "nParam", "activePar", "maxValue", "coef", "vcov",
     "hessian", "stdEr"),
   maxLik:::AIC.maxLik, x, maxLik:::AIC.maxLik(m), "AIC.maxLik"
)

## ------------------------------------------------------ inference methods
expect_no_dispatch(
   c("coef", "vcov", "stdEr", "activePar", "hessian", "nParam", "maxValue"),
   maxLik:::confint.maxLik, x, maxLik:::confint.maxLik(m),
   "confint.maxLik"
)

expect_no_dispatch(
   c("vcov", "nObs", "nobs", "stdEr", "activePar", "hessian", "nParam",
     "coef"),
   maxLik:::bread.maxLik, xObs, maxLik:::bread.maxLik(mObs),
   "bread.maxLik"
)

expect_no_dispatch(
   c("coef", "gradient", "vcov", "activePar"),
   maxLik:::estfun.maxLik, xObs, maxLik:::estfun.maxLik(mObs),
   "estfun.maxLik"
)

expect_no_dispatch(
   c("nObs", "nobs", "coef", "activePar"),
   maxLik:::nobs.maxLik, xObs, maxLik:::nobs.maxLik(mObs),
   "nobs.maxLik"
)

expect_no_dispatch(
   c("nObs", "nobs", "coef", "activePar"),
   maxLik:::nObs.maxLik, xObs, maxLik:::nObs.maxLik(mObs),
   "nObs.maxLik"
)

## activePar.default() must keep its fallback for arbitrary objects
lmFit <- lm(y ~ 1)
expect_equal(
   activePar(lmFit), rep(TRUE, length(coef(lmFit))),
   info = "activePar.default() fallback for arbitrary objects"
)

expect_no_dispatch(
   c("hessian", "activePar", "coef", "vcov", "nParam", "stdEr"),
   function(object) capture.output(maxLik:::condiNumber.maxLik(object)),
   x, capture.output(maxLik:::condiNumber.maxLik(m)),
   "condiNumber.maxLik"
)

## ---------------------------------------------------------- summary.maxim()
expect_no_dispatch(
   c("maxValue", "maximType", "nIter", "returnCode", "returnMessage",
     "gradient", "hessian", "coef", "activePar", "vcov", "stdEr", "nParam"),
   maxLik:::summary.maxim, xmx, maxLik:::summary.maxim(mx),
   "summary.maxim"
)

## ------------------------------------------------------- tidy and glance
if(requireNamespace("tibble", quietly = TRUE)) {
   expect_no_dispatch(
      c("summary", "coef", "vcov", "stdEr", "activePar", "hessian",
        "nParam", "logLik", "AIC"),
      maxLik:::tidy.maxLik, x, maxLik:::tidy.maxLik(m), "tidy.maxLik"
   )

   expect_no_dispatch(
      c("coef", "gradient", "maxValue", "nIter", "returnCode",
        "returnMessage", "activePar"),
      maxLik:::tidy.maxim, xmx, maxLik:::tidy.maxim(mx), "tidy.maxim"
   )

   expect_no_dispatch(
      c("logLik", "nObs", "nobs", "AIC", "nParam", "activePar",
        "maxValue", "coef", "vcov", "hessian", "stdEr"),
      maxLik:::glance.maxLik, xObs, maxLik:::glance.maxLik(mObs),
      "glance.maxLik"
   )

   expect_no_dispatch(
      c("maxValue", "nIter", "gradient", "returnCode", "returnMessage",
        "coef", "activePar"),
      maxLik:::glance.maxim, xmx, maxLik:::glance.maxim(mx),
      "glance.maxim"
   )
}

## ------------------------------------------------------ public S3 dispatch
## The public API must still dispatch to subclass methods.
expect_equal(
   with_s3_methods(
      list(coef.foo = function(object, ...) "foo"),
      coef(x)
   ),
   "foo",
   info = "coef(x) should still dispatch to coef.foo(x)"
)

expect_equal(
   with_s3_methods(
      list(maxValue.foo = function(x, ...) "foo"),
      maxValue(x)
   ),
   "foo",
   info = "maxValue(x) should still dispatch to maxValue.foo(x)"
)

expect_equal(
   with_s3_methods(
      list(returnCode.foo = function(x, ...) "foo"),
      returnCode(x)
   ),
   "foo",
   info = "returnCode(x) should still dispatch to returnCode.foo(x)"
)

expect_equal(
   with_s3_methods(
      list(nParam.foo = function(x, ...) "foo"),
      nParam(x)
   ),
   "foo",
   info = "nParam(x) should still dispatch to nParam.foo(x)"
)
