tvec = linspace(0,1,21)';

xvec = sin(2*pi*tvec);

sigerr = 0.001;

yvec = xvec + randn(21,1)*sigerr;

argvals    = tvec;
y          = xvec;
wtvec      = [];
covariates = [];
fdnames    = [];

nbasis = 11;
basisobj = create_bspline_basis([0,1],nbasis);
fdParobj = basisobj;

%  for lambda == 0

[fdobj, beta, df, gcv, SSE, penmat, y2cMap] = ...
                 smooth_basis_LS(argvals, y, fdParobj);

df

%  test for large value of lambda

lambda = 1e-4;

fdParobj = fdPar(basisobj, 2, lambda);

[fdobj, beta, df, gcv, SSE, penmat, y2cMap] = ...
                 smooth_basis_LS(argvals, y, fdParobj);

%  test with covariate

zvec = zeros(21,1);
zvec(11) = 1;
beta = 0.5;

yvec = yvec + zvec*beta;
y = yvec;

lambda = 1e-8;
fdParobj = fdPar(basisobj, 2, lambda);

[fdobj, beta, df, gcv, SSE, penmat, y2cMap] = ...
                 smooth_basis_LS(argvals, y, fdParobj, ...
                                 'covariates', zvec);

beta

%  test with multiple records

xmat = [sin(2*pi*tvec),cos(2*pi*tvec)];

sigerr = 0.1;

ymat = xmat + randn(21,2)*sigerr;

beta = [-1,1];
ymat(:,1) = ymat(:,1) + zvec*beta(1);
ymat(:,2) = ymat(:,2) + zvec*beta(2);

y = ymat;

lambda = 1e-4;
fdParobj = fdPar(basisobj, 2, lambda);

[fdobj, beta, df, gcv, SSE, penmat, y2cMap] = ...
                 smooth_basis_LS(argvals, y, fdParobj, ...
                                 'covariates', zvec);

beta

plotfit_fd(y, argvals, fdobj)



