source("ycheck.R")
source("argcheck.R")
source("constrain.R")
source("fdParcheck.R")
source("glm.fda.R")
source("smooth.GLM.R")

printf <- function(...) cat(sprintf(...))


#  tests for function smooth.GLM

#  ----------------  normal link with no covariates  --------------------

library(pracma)
library(fda)

n       = 101
argvals = linspace(0,2*pi,n)
y0      = sin(2*argvals)
y02     = cos(2*argvals)
sig     = 0.2
y       = y0 + sig * randn(n,1)
y       = cbind(y, y02 + sig * randn(n,1))

basisobj = create.bspline.basis(c(0,2*pi),n+2)

basismat = eval.basis(argvals, basisobj)

Lfdobj = int2Lfd(2)
penmat = eval.penalty(basisobj,Lfdobj)

lambda  = 1e-1
lamRmat = lambda * penmat

glmresult = glm.fda(basismat, y, 'normal', lamRmat)

fdobj = fd(glmresult$coef,basisobj)

plotfit.fd(y, argvals, fdobj)

fdParobj = fdPar(basisobj, Lfdobj, lambda)

smoothResult = smooth.GLM(argvals, y, fdParobj ,family = 'normal')
plotfit.fd(y, argvals, fdobj)

smoothResult$df

smoothResult$gcv

smoothResult$SSE

smoothResult$dev

smoothResult$beta

#  ----------------  normal link with covariates  --------------------

n       = 101
argvals = linspace(0,2*pi,n)
y0      = sin(2*argvals)
sig     = 0.2
y       = y0 + sig * randn(n,1)
sigcov  = 0.1
covariates = randn(n,1)
beta    = 1
y       = y + covariates*beta

basisobj = create.bspline.basis(c(0,2*pi),11)

basismat = eval.basis(argvals, basisobj)

Lfdobj = int2Lfd(2)
penmat = eval.penalty(basisobj,Lfdobj)

lambda  = 1
lamRmat = lambda * penmat

fdParobj = fdPar(basisobj, Lfdobj, lambda)

smoothResult = smooth.GLM(argvals, y, fdParobj, family = 'normal', covariates=covariates)

plotfit.fd(y, argvals, fdobj)

smoothResult$beta

smoothResult$df

smoothResult$gcv

smoothResult$SSE

smoothResult$dev

#  ----------------  binomial link with no covariate  --------------------

n       = 501
argvals = linspace(0,1,n)
y0      = sin(4*pi*argvals)
sig     = 0.5
y = matrix(0,n,1)
y[y0 + sig*randn(n,1) >= 0.0] = 1

basisobj = create.bspline.basis(c(0,1),13)

basismat = eval.basis(argvals, basisobj)

Lfdobj = int2Lfd(2)
penmat = eval.penalty(basisobj,Lfdobj)

lambda  = 1e-4
lamRmat = lambda*penmat

glmResult = glm.fda(basismat, y, 'binomial', lamRmat)

fdobj = fd(glmResult$coef,basisobj)

plotfit.fd(y, argvals, fdobj)

fdParobj = fdPar(basisobj, Lfdobj, lambda)

smoothResult = smooth.GLM(argvals, y, fdParobj, family = 'binomial')

plot(fdobj)

smoothResult$beta

smoothResult$df

smoothResult$gcv

smoothResult$SSE

smoothResult$dev

# stats{}

#  ----------------  poisson link with no covariates  --------------------

n = 101
argvals = linspace(0,2*pi,n)'
y0 = sin(2*argvals)
sig = 0.2
y = y0 + sig.*randn(n,1)
y = exp(y)

basisobj = create.bspline.basis([0,2*pi],53)

lambda = 1e-1
Lfdobj = int2Lfd(2)
fdParobj = fdPar(basisobj, Lfdobj, lambda)

[fdobj, beta, df, gcv, SSE, dev, stats] = ...
            smooth.basis.GLM(argvals, y, fdParobj, ...
                             'family', 'poisson')

plotfit.fd(log(y), argvals, fdobj)

beta

df

gcv

SSE

dev

stats{}







