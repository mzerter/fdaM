%  Test problems for LMdensity

N = 101;

zmat = [ones(N,1), randn(N,1)];
ncov = size(zmat,2);

bvec = randn(size(zmat,2),1);

y0 = zmat*bvec;

%  re-define zmat so that residuals can be centered

zsum = sum(zmat)';

[Qz, Rz] = qr(zsum);
Qres = Qz(:,2:ncov);

zmatctr = zmat*Qres;

sigma = 1;

y = y0 + randn(N,1).*sigma;

Wrange  = [-2.0,2.5].*sigma;
Wnbasis = 5;
Wbasis  = create_bspline_basis(Wrange, Wnbasis);
Wcoef   = randn(Wnbasis,1);
Wfd     = fd(Wcoef, Wbasis);
Wlambda = 1e-2;
WfdPar  = fdPar(Wfd, int2Lfd(2), Wlambda);

zerobasmat = zerobasis(Wnbasis);

beta0 = Qres'*bvec;
sigma0 = sigma;
res0 = y - zmatctr*beta0;
res0 = res0 - mean(res0);

conv = 1e-2;
iterlim = 20;
active = 1:Wnbasis;
dbglev = 1;

[Wfdobj, C, hmat, Fstr, iter, iterhist] = ...
    density_fd(res0/sigma, WfdPar, conv, iterlim, active, dbglev);

Wcoef   = getcoef(Wfdobj);
Wfd     = fd(Wcoef, Wbasis);
Wlambda = 1e-2;
WfdPar  = fdPar(Wfd, int2Lfd(2), Wlambda);

[Wfdobj, beta, Cval, res, hmat, Fstr, iternum, iterhist] = ...
    LMdens_fd(y, WfdPar, zmat, beta0, sigma0, ...
                  conv, iterlim, active, dbglev);

cvec = getcoef(Wfdobj)
zerobasmat'*cvec

Cval

[beta, beta0]

plot(Wfdobj)

zfine = linspace(Wrange(1),Wrange(2),101)';

Wvec = eval_fd(zfine, Wfdobj);

Pvec = exp(Wvec)./Cval;

plot(zfine, Pvec, '-', res/sigma0, zeros(nobs,1), 'r.')




