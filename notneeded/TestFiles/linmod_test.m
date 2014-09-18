%  ------------------------------------
%      set up problem  for linmod
%  ------------------------------------

rngs = [0,1];
rngt = [0,2*pi];

nfine = 501;
sfine = linspace(rngs(1),rngs(2),nfine)';
tfine = linspace(rngt(1),rngt(2),nfine)';

N = 100;

snbasis = 8;
sbasis  = create_bspline_basis(rngs,snbasis);
smat    = eval_basis(sfine, sbasis);

tnbasis = 10;
tbasis  = create_bspline_basis(rngt,tnbasis);
tmat    = eval_basis(tfine, tbasis);

alphanbasis = 20;
alphabasis  = create_bspline_basis(rngt,alphanbasis);
sigma       = 1;
alphacoef0  = sigma.*randn(alphanbasis,1);
alphafd0    = fd(alphacoef0, alphabasis);

% alphanbasis = 1;
% alphabasis  = create_constant_basis(rngt);
% alphacoef0  = 0;
% alphafd0    = fd(alphacoef0, alphabasis);

betacoef0 = randn(snbasis,tnbasis);
betafd0   = bifd(betacoef0,sbasis,tbasis);

xnbasis = 53;
xbasis  = create_bspline_basis(rngs, xnbasis);
sigmaX  = 10;
xfd     = fd(sigmaX.*randn(xnbasis,N),xbasis);

sigmaE  = 0.1;
enbasis = 53;
ebasis  = create_bspline_basis(rngt, enbasis);
efd     = fd(randn(enbasis,N).*sigmaE, ebasis);

Hmat     = inprod(xfd, sbasis);

xcoef    = (Hmat*betacoef0)';
xbetafd0 = fd(xcoef, tbasis);

yhatmat = eval_fd(tfine, alphafd0)*ones(1,N) + ...
          eval_fd(tfine, xbetafd0);
yhatfd0 = smooth_basis(tfine, yhatmat, ebasis);
ymat    = yhatmat + eval_fd(tfine, efd);
yfd     = smooth_basis(tfine, ymat, ebasis);

betacell    = cell(1,2);
betacell{1} = fdPar(alphafd0,  2, 1e4);
betacell{2} = bifdPar(betafd0, 2, 2, 0, 1e4);

linmodstr = linmod(xfd, yfd, betacell);

alphafd = linmodstr.alpha;
betafd  = linmodstr.beta;
yhatfd  = linmodstr.yhat;

display([getcoef(alphafd),alphacoef0])

subplot(1,1,1)
plot(alphafd)
lhdl = line(alphafd0);
set(lhdl, 'LineStyle', '--')

betacoef = getcoef(betafd);

display([betacoef', betacoef0'])

subplot(2,1,1)
plot(yhatfd)
subplot(2,1,2)
plot(yhatfd0)

bmat = smat*betacoef*tmat';

subplot(1,1,1)
contour(bmat)
axis('square')
xlabel('\fontsize{13} t')
ylabel('\fontsize{13} s')

surf(bmat)
xlabel('\fontsize{13} t')
ylabel('\fontsize{13} s')

%  ------------------------------------
%  set up problem  for fRegress
%  ------------------------------------

xfdcell = cell(1,1+snbasis);

xfdcell{1} = ones(N,1);
for j=2:(snbasis+1)
    xfdcell{j} = Hmat(:,j-1);
end

betacellfR = cell(1,1+snbasis);
betacellfR{1} = fdPar(alphafd0);
tfd0 = fd(ones(tnbasis,1),tbasis);
tfdPar = fdPar(tfd0);
for j=2:(snbasis+1)
    betacellfR{j} = tfdPar;
end

fRegressStr = fRegress(yfd, xfdcell, betacellfR);

betaestcell = fRegressStr.betahat;

yhatfd = fRegressStr.yhat;

subplot(2,1,1)
plot(yhatfd)
subplot(2,1,2)
plot(yhatfd0)

alphafdPar = betaestcell{1};
alphafd = getfd(alphafdPar);

display([getcoef(alphafd),alphacoef0])

subplot(1,1,1)
plot(alphafd)
lhdl = line(alphafd0);
set(lhdl, 'LineStyle', '--')

for j=1:snbasis
    display(getcoef(getfd(betaestcell{j+1}))')
end

display(betacoef0)




