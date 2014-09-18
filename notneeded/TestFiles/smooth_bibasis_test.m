snbasis = 15;
sbasis  = create_bspline_basis([0,1],snbasis);
tnbasis = 15;
tbasis  = create_bspline_basis([0,1],tnbasis);
lambdas = 0e-4;
fdPars  = fdPar(sbasis, 2, lambdas);
lambdat = 1e-3;
fdPart  = fdPar(tbasis, 2, lambdat);

ns = 15;
nt = 15;

sarg = linspace(0,1,ns)';
targ = linspace(0,1,nt)';

sigma = 0.2;
y0 = sin(2*pi*sarg)*cos(2*pi*targ)';
y  = y0 + randn(ns,nt).*sigma;

[bifdobj, df, gcv, coef, SSE, y2cMap] = ...
    smooth_bibasis(sarg, targ, y, fdPars, fdPart, fdnames);

bimat = eval_bifd(sarg, targ, bifdobj);

figure(1)
surf(bimat)
xlabel('t')
ylabel('s')
figure(2)
surf(y0)
xlabel('t')
ylabel('s')
figure(3)
surf(y)
xlabel('t')
ylabel('s')



