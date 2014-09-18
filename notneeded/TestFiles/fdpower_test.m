%  code for testing pointwise product of two functional data objects

%  Spencer's test in fdarm-ch04.R

bspl2 = create_bspline_basis([0,1],2,2);
plot(bspl2)

getbasistype(bspl2)

tstFn0 = fd([-1; 2], bspl2);
tstFn1 = fd([ 1; 3], bspl2);

subplot(2,1,1)
plot(tstFn0)
subplot(2,1,2)
plot(tstFn1)

subplot(1,1,1)
fdsumobj = tstFn0 + tstFn1;
plot(fdsumobj)

subplot(1,1,1)
fddifobj = tstFn0 - tstFn1;
plot(fddifobj)

fdprdobj = tstFn0.*tstFn1;
plot(fdprdobj)

fdsqrobj = fdsumobj.^10;
subplot(2,1,1)
plot(fdsumobj)
subplot(2,1,2)
plot(fdsqrobj)

basisobj = create_bspline_basis([0,1],13);
% coefmat  = exp(randn(13,1));
coefmat  = randn(13,1);
fdobj    = fd(coefmat, basisobj);

a = 2;
fdpowerobj = fdobj.^a;
subplot(2,1,1)
plot(fdobj)
subplot(2,1,2)
plot(fdpowerobj)

