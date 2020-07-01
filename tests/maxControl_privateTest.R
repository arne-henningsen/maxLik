### Does maxControl stuff behave?
### no need to test it on CRAN, hence private test
### 
### test for:
### 1. create maxControl object
### 2. SGA_batchSize NULL
### 3. negative batch size
### 4. more than 1 batch size
### SG_clip: NULL, negative, more than one

library(maxLik)
library(testthat)

### ---------- create maxControl object
maxControl(tol=1e-4, lambdatol=1e-5, qrtol=1e-6, qac="marquardt",
           marquardt_lambda0=0.1, marquardt_lambdaStep=3, marquardt_maxLambda=1e10,
           nm_alpha=2, nm_beta=1, nm_gamma=4,
           sann_temp=5, sann_tmax=100, sann_randomSeed=1,
           SGA_momentum=0.9, Adam_momentum1=0.5, Adam_momentum2=0.55,
           SG_learningRate=0.5, SG_batchSize=10, SG_clip=1000, 
           SG_patience=7, SG_patienceStep=10,
           iterlim=10, printLevel=3)

### ---------- SG_batchSize
maxControl(SG_batchSize=NULL)  # should work
try(maxControl(SG_batchSize=-1))  # should fail
try(maxControl(SG_batchSize=2:3))  # should fail

maxControl(SG_clip=NULL)  # works
try(maxControl(SG_clip=-1))  # fails
try(maxControl(SG_clip=2:3))  # fails

try(maxControl(Adam_momentum1=NA))  # should fail w/'NA in Adam_momentum'
