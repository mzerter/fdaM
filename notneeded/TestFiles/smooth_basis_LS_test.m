%  tests for function smooth_basis_LS

%  ----------------  A call with no covariates  ------------------------

n = 101;
argvals = linspace(0,2*pi,n)';
y0 = sin(2*argvals);
sig = 0.2;
y = y0 + sig.*randn(n,1);

basisobj = create_bspline_basis([0,2*pi],n+2);

lambda = 1e-2;
Lfdobj = int2Lfd(2);
fdParobj = fdPar(basisobj, Lfdobj, lambda);

[fdobj, df, gcv, SSE, penmat, y2cMap] = ...
                      smooth_basis_LS(argvals, y, fdParobj);

[fdobj1, df1, gcv1, SSE1, penmat1, y2cMap1] = ...
                      smooth_basis(argvals, y, fdParobj);
                  
plot(fdobj1-fdobj)

df1-df

gcv1-gcv

SSE1-SSE

max(max(abs(y2cMap1-y2cMap)))

%  ----------------  A call with covariates  ------------------------

n = 101;
argvals = linspace(0,2*pi,n)';
y0 = sin(2*argvals);
sig = 0.2;
y = y0 + sig.*randn(n,1);
sigcov = 0.1;
covariates = sigcov.*randn(n,1);
beta = 1;
y = y + covariates*beta;

basisobj = create_bspline_basis([0,2*pi],n+2);

lambda   = 1e-2;
Lfdobj   = int2Lfd(2);
fdParobj = fdPar(basisobj, Lfdobj, lambda);

[fdobj, beta, df, gcv, SSE, penmat, y2cMap] = ...
          smooth_basis_LS(argvals, y, fdParobj, 'covariates', covariates);

beta

plotfit_fd(y, argvals, fdobj)

df

gcv

SSE


