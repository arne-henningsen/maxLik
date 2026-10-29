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


## Create coef.foo()-type functions that fail with the appripriate
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

 
## Simple fitted model
set.seed(123)
y <- rnorm(20, mean = 1, sd = 2)
ll <- function(theta) {
   sum(dnorm(y, mean = theta[1], sd = theta[2], log = TRUE))
}
m <- maxLik(ll, start = c(1, 2))
## Reference results: correct summary with 'maxLik' object
x_plain <- m
reference <- maxLik:::summary.maxLik(x_plain)
## 'x': otherwise the same, except the class
x <- m
class(x) <- c("foo", class(x))

## Individual poisons: each of these generics is called somewhere in the
## summary -> stdEr -> vcov -> activePar/hessian/nParam chain.
for(g in c("coef", "activePar", "stdEr", "vcov", "hessian")) {
   result <- tryCatch(
      with_poison(g, maxLik:::summary.maxLik(x)),
      error = function(e) e
   )
   expect_false(
      inherits(result, "error"),
      info = paste0("summary.maxLik() called ", g, ".foo()")
   )
   if(!inherits(result, "error")) {
      expect_equal(
         result, reference,
         info = paste0("summary.maxLik() unaffected by ", g, ".foo()")
      )
   }
}

## nParam() is only reached through vcov.maxLik() when the Hessian is absent
## So test nParam() here separately w/o Hessian
x_nohess <- x
x_nohess$hessian <- NULL
x_nohess_plain <- x_plain
x_nohess_plain$hessian <- NULL
reference_nohess <- maxLik:::summary.maxLik(x_nohess_plain)

result <- tryCatch(
   with_poison("nParam", maxLik:::summary.maxLik(x_nohess)),
   error = function(e) e
)
expect_false(
   inherits(result, "error"),
   info = "summary.maxLik() called nParam.foo() when Hessian is missing"
)
if(!inherits(result, "error")) {
   expect_equal(
      result, reference_nohess,
      info = "summary.maxLik() unaffected by nParam.foo() without Hessian"
   )
}

## Poison all methods and test summary.maxLik() ...
result <- tryCatch(
   with_poison(
      c("coef", "activePar", "stdEr", "vcov", "hessian"),
      maxLik:::summary.maxLik(x)
   ),
   error = function(e) e
)
expect_false(
   inherits(result, "error"),
   info = "summary.maxLik() must not use subclass accessor methods"
)
if(!inherits(result, "error")) {
   expect_equal(result, reference)
}
## ... and test summary()
result <- tryCatch(
   with_poison(
      c("coef", "activePar", "stdEr", "vcov", "hessian"),
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

## The public S3 API must still dispatch to the subclass method.
expect_equal(
   with_s3_methods(
      list(coef.foo = function(object, ...) "foo"),
      coef(x)
   ),
   "foo",
   info = "coef(x) should still dispatch to coef.foo(x)"
)
