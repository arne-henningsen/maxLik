### Does maxControl stuff behave?
### test for:
### 1. create maxControl object

library(maxLik)
library(testthat)

maxControl(tol=1e-4, lambdatol=1e-5, qrtol=1e-6, qac="marquardt",
           marquardt_lambda0=0.1, marquardt_lambdaStep=3, marquardt_maxLambda=1e10,
           nm_alpha=2, nm_beta=1, nm_gamma=4,
           sann_temp=5, sann_tmax=100, sann_randomSeed=1,
           iterlim=10, printLevel=3)
