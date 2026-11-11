### Test that maxLik methods do not dispatch to subclass methods internally
##
## Unlike the optimizer tests in 'test-optimizers.R', these tests are
## deliberately also run on CRAN: they are fast, they compare a poisoned
## call with an unpoisoned one within the same run, and they test S3
## dispatch and accessor behaviour rather than numerical optimization.

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
   ## Test these generics
   c("coef", "activePar", "maxValue", "maximType", "nIter",
     "returnCode", "returnMessage", "vcov", "hessian", "nParam", "stdEr"),
   ## Run this function with poison
   function(object) capture.output(maxLik:::print.maxLik(object)),
   ## Use this for the function argument:
   x,
   ## Reference value: should get the same result
   capture.output(maxLik:::print.maxLik(m)),
   ## Informative label to be printed
   "print.maxLik"
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

## ---------------------------------------------------- hardening fixtures
## The generic list used by the branch-covering tests below.
hardeningGenerics <- c(
   "coef", "vcov", "stdEr", "hessian", "activePar", "nParam",
   "maxValue", "maximType", "nIter", "returnCode", "returnMessage",
   "gradient", "logLik", "AIC", "nobs", "nObs"
)

## ---------------------------------------------- fixed parameters
mFix <- maxLik(ll, start = c(1, 2), fixed = c(FALSE, TRUE))
xFix <- as_foo(mFix)

expect_no_dispatch(
   hardeningGenerics,
   maxLik:::summary.maxLik, xFix, maxLik:::summary.maxLik(mFix),
   "summary.maxLik with fixed parameters"
)

expect_no_dispatch(
   hardeningGenerics,
   function(object) capture.output(maxLik:::print.maxLik(object)),
   xFix, capture.output(maxLik:::print.maxLik(mFix)),
   "print.maxLik with fixed parameters"
)

expect_no_dispatch(
   hardeningGenerics,
   maxLik:::vcov.maxLik, xFix, maxLik:::vcov.maxLik(mFix),
   "vcov.maxLik with fixed parameters"
)

expect_no_dispatch(
   hardeningGenerics,
   maxLik:::stdEr.maxLik, xFix, maxLik:::stdEr.maxLik(mFix),
   "stdEr.maxLik with fixed parameters"
)

expect_no_dispatch(
   hardeningGenerics,
   maxLik:::confint.maxLik, xFix, maxLik:::confint.maxLik(mFix),
   "confint.maxLik with fixed parameters"
)

## ---------------------------------------------- missing Hessian
mNoHess <- m
mNoHess$hessian <- NULL
xNoHess <- as_foo(mNoHess)

expect_no_dispatch(
   hardeningGenerics,
   maxLik:::summary.maxLik, xNoHess, maxLik:::summary.maxLik(mNoHess),
   "summary.maxLik without Hessian"
)

expect_no_dispatch(
   hardeningGenerics,
   maxLik:::vcov.maxLik, xNoHess, maxLik:::vcov.maxLik(mNoHess),
   "vcov.maxLik without Hessian"
)

expect_no_dispatch(
   hardeningGenerics,
   maxLik:::stdEr.maxLik, xNoHess, maxLik:::stdEr.maxLik(mNoHess),
   "stdEr.maxLik without Hessian"
)

expect_no_dispatch(
   hardeningGenerics,
   maxLik:::confint.maxLik, xNoHess, maxLik:::confint.maxLik(mNoHess),
   "confint.maxLik without Hessian"
)

## ---------------------------------------------- singular Hessian
## Zero one symmetric row/column pair of an existing Hessian.  This keeps
## the matrix symmetric and produces exactly one zero eigenvalue.
mSing <- m
H <- mSing$hessian
H[1, ] <- 0
H[, 1] <- 0
mSing$hessian <- H
xSing <- as_foo(mSing)

expect_no_dispatch(
   hardeningGenerics,
   maxLik:::vcov.maxLik, xSing, maxLik:::vcov.maxLik(mSing),
   "vcov.maxLik with singular Hessian"
)

expect_no_dispatch(
   hardeningGenerics,
   maxLik:::stdEr.maxLik, xSing, maxLik:::stdEr.maxLik(mSing),
   "stdEr.maxLik with singular Hessian"
)

expect_no_dispatch(
   hardeningGenerics,
   maxLik:::summary.maxLik, xSing, maxLik:::summary.maxLik(mSing),
   "summary.maxLik with singular Hessian"
)

expect_no_dispatch(
   hardeningGenerics,
   maxLik:::confint.maxLik, xSing, maxLik:::confint.maxLik(mSing),
   "confint.maxLik with singular Hessian"
)

## ---------------------------------------------- precomputed varcovar
mVar <- m
mVar$varcovar <- diag(c(0.01, 4))
dimnames(mVar$varcovar) <- dimnames(mVar$hessian)
xVar <- as_foo(mVar)

expect_no_dispatch(
   hardeningGenerics,
   maxLik:::vcov.maxLik, xVar, maxLik:::vcov.maxLik(mVar),
   "vcov.maxLik with precomputed varcovar"
)

expect_no_dispatch(
   hardeningGenerics,
   maxLik:::stdEr.maxLik, xVar, maxLik:::stdEr.maxLik(mVar),
   "stdEr.maxLik with precomputed varcovar"
)

expect_no_dispatch(
   hardeningGenerics,
   maxLik:::summary.maxLik, xVar, maxLik:::summary.maxLik(mVar),
   "summary.maxLik with precomputed varcovar"
)

expect_no_dispatch(
   hardeningGenerics,
   maxLik:::confint.maxLik, xVar, maxLik:::confint.maxLik(mVar),
   "confint.maxLik with precomputed varcovar"
)

## ---------------------------------------------- constrained maxLik
## Synthetic constraints, with return code and message.
mCon <- m
mCon$constraints <- list(
   type = "linear equality constraints",
   code = 1,
   message = "convergence",
   outer.iterations = 2,
   barrier.value = 0.001
)
xCon <- as_foo(mCon)
sCon <- maxLik:::summary.maxLik(mCon)
xsCon <- as_foo(sCon)

expect_no_dispatch(
   hardeningGenerics,
   maxLik:::summary.maxLik, xCon, sCon,
   "summary.maxLik with constraints"
)

expect_no_dispatch(
   hardeningGenerics,
   function(object) capture.output(maxLik:::print.summary.maxLik(object)),
   xsCon, capture.output(maxLik:::print.summary.maxLik(sCon)),
   "print.summary.maxLik with constraints"
)

## Synthetic constraints without code/message, as for 'constrOptim'.
mConNoMsg <- m
mConNoMsg$constraints <- list(
   type = "linear inequality constraints",
   outer.iterations = 3,
   barrier.value = 0.01
)
xConNoMsg <- as_foo(mConNoMsg)
sConNoMsg <- maxLik:::summary.maxLik(mConNoMsg)
xsConNoMsg <- as_foo(sConNoMsg)

expect_no_dispatch(
   hardeningGenerics,
   maxLik:::summary.maxLik, xConNoMsg, sConNoMsg,
   "summary.maxLik with constraints without code/message"
)

expect_no_dispatch(
   hardeningGenerics,
   function(object) capture.output(maxLik:::print.summary.maxLik(object)),
   xsConNoMsg, capture.output(maxLik:::print.summary.maxLik(sConNoMsg)),
   "print.summary.maxLik with constraints without code/message"
)

## A summary object without 'NActivePar' must compute it from 'fixed'
sNoNAct <- s
sNoNAct$NActivePar <- NULL
xsNoNAct <- as_foo(sNoNAct)

expect_no_dispatch(
   hardeningGenerics,
   function(object) capture.output(maxLik:::print.summary.maxLik(object)),
   xsNoNAct, capture.output(maxLik:::print.summary.maxLik(sNoNAct)),
   "print.summary.maxLik without NActivePar"
)

## ---------------------------------------------- constrained maxim
mxCon <- mx
mxCon$constraints <- list(
   type = "linear equality constraints",
   code = 1,
   message = "convergence",
   outer.iterations = 2,
   barrier.value = 0.001
)
xmxCon <- as_foo(mxCon)
smxCon <- maxLik:::summary.maxim(mxCon)
xsmxCon <- as_foo(smxCon)

expect_no_dispatch(
   hardeningGenerics,
   maxLik:::summary.maxim, xmxCon, smxCon,
   "summary.maxim with constraints"
)

expect_no_dispatch(
   hardeningGenerics,
   function(object) capture.output(maxLik:::print.summary.maxim(object)),
   xsmxCon, capture.output(maxLik:::print.summary.maxim(smxCon)),
   "print.summary.maxim with constraints"
)

## ---------------------------------------------- maxim Hessian and last step
expect_no_dispatch(
   hardeningGenerics,
   function(object) maxLik:::summary.maxim(object, hessian = TRUE),
   xmx, maxLik:::summary.maxim(mx, hessian = TRUE),
   "summary.maxim with Hessian"
)

mxLast <- mx
mxLast$code <- 3
mxLast$names <- c("theta1", "theta2")
mxLast$last.step <- list(
   theta0 = c(0, 0),
   theta1 = c(1, 2),
   f0 = -5
)
xmxLast <- as_foo(mxLast)
smxLast <- maxLik:::summary.maxim(mxLast, unsucc.step = TRUE)
xsmxLast <- as_foo(smxLast)

expect_no_dispatch(
   hardeningGenerics,
   function(object) maxLik:::summary.maxim(object, unsucc.step = TRUE),
   xmxLast, smxLast,
   "summary.maxim with unsuccessful last step"
)

expect_no_dispatch(
   hardeningGenerics,
   function(object) capture.output(maxLik:::print.summary.maxim(object)),
   xsmxLast, capture.output(maxLik:::print.summary.maxim(smxLast)),
   "print.summary.maxim with unsuccessful last step"
)

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

## ------------------------------------------------ static no-dispatch guard
## The tests above check the call paths actually exercised by the test
## objects.  This block additionally inspects the method bodies, so a
## forbidden generic call is caught even if that branch is not executed.
## The functions to be inspected are derived from the package namespace,
## so newly added methods and accessors are covered automatically.  Only
## the list of forbidden generics and the list of exceptions below need to
## be maintained by hand.
if(requireNamespace("codetools", quietly = TRUE)) {
   ns <- asNamespace("maxLik")

   ## The public extractor generics that maxLik's own code must never call.
   ## This list is a policy decision and cannot be derived: whenever a new
   ## public extractor generic is added to the package, it has to be added
   ## here.  'hardeningGenerics' above is the shorter list used by the
   ## dynamic poison tests, where a generic is only useful if the test
   ## object actually reaches the corresponding code path.
   staticForbidden <- union(hardeningGenerics,
      c("objectiveFn", "storedParameters", "storedValues"))

   ## The internal accessors, derived from the namespace so that a new
   ## accessor is covered automatically.  R's own namespace bookkeeping is
   ## excluded: entries such as '.__C__maxLik', '.__T__show:methods',
   ## '.__S3MethodsTable_.' and '.__global__' all start with a dot
   ## followed by an underscore.  The remaining non-accessors are the
   ## namespace startup functions, '.packageName' and
   ## '.S3MethodsClasses'.
   ## No letter class like [a-z] is used on purpose: in POSIX regexes such
   ## a range is a range of collating elements, so its extent depends on
   ## LC_COLLATE.
   dotNames <- grep("^\\.[^_]", ls(ns, all.names = TRUE), value = TRUE)
   internalAccessors <- setdiff(dotNames,
      c(".onAttach", ".onLoad", ".onUnload", ".packageName",
        ".S3MethodsClasses"))
   ## Anything that is not a function (further bookkeeping) is skipped.
   internalAccessors <- internalAccessors[
      vapply(internalAccessors, function(name)
         is.function(get(name, envir = ns)), logical(1))
   ]

   ## All S3 methods registered by the package, derived from the namespace
   ## so that a new 'S3method()' entry in NAMESPACE is covered
   ## automatically.  'sumt' is a plain function rather than an S3 method
   ## for 'maxLik' or 'maxim' objects, and hence has to be listed.
   internalMethods <- unique(getNamespaceInfo(ns, "S3methods")[, 3])
   internalMethods <- union(internalMethods, "sumt")

   ## The S4 methods defined by the package, derived from the method
   ## tables of its namespace.  Entries that do not belong to the maxLik
   ## namespace are skipped, e.g. the 'ANY' method of the 'nobs' generic
   ## imported from 'stats'.
   s4Methods <- list()
   for(table in grep("^\\.__T__", ls(ns, all.names = TRUE), value = TRUE)) {
      generic <- sub(":.*$", "", sub("^\\.__T__", "", table))
      entries <- get(table, envir = ns)
      for(sig in ls(entries, all.names = TRUE)) {
         method <- get(sig, envir = entries)
         if(is.function(method) && identical(topenv(environment(method)), ns)) {
            s4Methods[[paste0(generic, "(", sig, ")")]] <- method
         }
      }
   }
   ## The S4 'show' methods for 'maxLik' and 'maxim' call the S3 generic
   ## 'print()' deliberately, see the comments in 'R/print.maxLik.R' and
   ## 'R/print.maxim.R'.  'print' is therefore not part of
   ## 'staticForbidden'.

   ## Exceptions: functions that may use a public generic on purpose.
   ## Any new exception needs a comment explaining why it is safe.
   allowed <- list(
      ## documented fallback for arbitrary objects
      "activePar.default" = "coef"
   )
   allowedFor <- function(name) {
      if(name %in% names(allowed)) allowed[[name]] else character()
   }

   ## Sanity checks: a filter that silently selects too little would let
   ## the loop below pass without checking anything.
   expect_true(length(internalAccessors) >= 20)
   expect_true(length(internalMethods) >= 40)
   expect_true(length(s4Methods) >= 5)
   expect_true(all(c(".field", ".vcov") %in% internalAccessors))
   expect_true(all(c("summary.maxLik", "print.summary.maxLik") %in%
                   internalMethods))
   expect_true(all(c("maxControl(maxim)", "show(maxLik)") %in%
                   names(s4Methods)))

   expectNoGenericCall <- function(name, fun, allowedCalls) {
      calls <- codetools::findGlobals(fun, merge = FALSE)$functions
      bad <- setdiff(intersect(calls, staticForbidden), allowedCalls)
      expect_equal(
         bad, character(),
         info = paste0(name, "() uses public generic(s): ",
                       paste(bad, collapse = ", "))
      )
   }

   for(name in c(internalAccessors, internalMethods)) {
      expectNoGenericCall(name, get(name, envir = ns), allowedFor(name))
   }
   for(signature in names(s4Methods)) {
      expectNoGenericCall(signature, s4Methods[[signature]], character())
   }

   ## Canary: 'activePar.default()' is the one known use of a public
   ## generic.  If this use is not found, the check itself is broken, e.g.
   ## because the behaviour of 'findGlobals()' has changed.
   canary <- codetools::findGlobals(get("activePar.default", envir = ns),
      merge = FALSE)$functions
   expect_true("coef" %in% canary)
} else {
   message("Skipping static no-dispatch checks: 'codetools' not available")
}

