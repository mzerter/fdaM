%  tests for function smooth_basis_LS

%  ----------------  normal link with no covariates  --------------------

n = 101;
argvals = linspace(0,2*pi,n)';
y0 = sin(2*argvals);
sig = 0.2;
y = y0 + sig.*randn(n,1);

basisobj = create_bspline_basis([0,2*pi],n+2);

lambda = 1e-2;
Lfdobj = int2Lfd(2);
fdParobj = fdPar(basisobj, Lfdobj, lambda);

[fdobj, beta, df, gcv, SSE, dev, stats] = ...
            smooth_basis_GLM(argvals, y, fdParobj, ...
                        'family', 'normal');
                  
plotfit_fd(y, argvals, fdobj)

beta

df

gcv

SSE

dev

stats{:}

%  ----------------  normal link with covariates  --------------------

n = 101;
argvals = linspace(0,2*pi,n)';
y0 = sin(2*argvals);
sig = 0.2;
y = y0 + sig.*randn(n,1);
sigcov = 0.1;
covariates = ones(n,1);
beta = 1;
y = y + covariates*beta;

basisobj = create_bspline_basis([0,2*pi],n+2);

lambda   = 1e-2;
Lfdobj   = int2Lfd(2);
fdParobj = fdPar(basisobj, Lfdobj, lambda);

[fdobj, beta, df, gcv, SSE, dev, stats] = ...
          smooth_basis_GLM(argvals, y, fdParobj, 'covariates', covariates);

beta

plotfit_fd(y, argvals, fdobj)

df

gcv

SSE

dev

stats{:}

%  ----------------  binomial link with no covariate  --------------------

n = 501;
argvals = linspace(0,1,n)';
y0 = sin(4*pi*argvals);
sig = 0.5;
y = y0 + sig.*randn(n,1);
y = y >= 0.0;

basisobj = create_bspline_basis([0,1],13);

lambda = 1e-4;
Lfdobj = int2Lfd(2);
fdParobj = fdPar(basisobj, Lfdobj, lambda);

varargin = cell(2,1);
varargin{1} = 'family';
varargin{2} = 'binomial';

[fdobj, beta, df, gcv, SSE, dev, stats] = ...
                      smooth_basis_GLM(argvals, y, fdParobj, ...
                                       'family', 'binomial');
                  
plot(fdobj)

beta

df

gcv

SSE

dev

stats{:}

%  ----------------  poisson link with no covariates  --------------------

n = 101;
argvals = linspace(0,2*pi,n)';
y0 = sin(2*argvals);
sig = 0.2;
y = y0 + sig.*randn(n,1);
y = exp(y);

basisobj = create_bspline_basis([0,2*pi],53);

lambda = 1e-1;
Lfdobj = int2Lfd(2);
fdParobj = fdPar(basisobj, Lfdobj, lambda);

[fdobj, beta, df, gcv, SSE, dev, stats] = ...
            smooth_basis_GLM(argvals, y, fdParobj, ...
                             'family', 'poisson');

plotfit_fd(log(y), argvals, fdobj)

beta

df

gcv

SSE

dev

stats{:}


