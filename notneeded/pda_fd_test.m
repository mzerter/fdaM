%  set up objects for examples

%  constant basis for estimating weight functions
cbasis = create_constant_basis([0,1]);
%  monomial basis: {1,t}  for estimating weight functions
mbasis = create_monomial_basis([0,1],2);
%  quartic spline basis with 54 basis functions for
%    defining functions to be analyzed
xbasis = create_bspline_basis([0,1],24,5);
%  set up functional parameter objects for weight bases
cfdPar = fdPar(cbasis);
mfdPar = fdPar(mbasis);
%  sampling points over [0,1]
tvec = linspace(0,1,101)';

%  Example 1:  a single first order constant coefficient unforced equation
%     Dx = -4.*x  for  x(t) = exp(-4t)

beta    = 4;
xvec    = exp(-beta.*tvec);
xfd     = smooth_basis(tvec, xvec, xbasis);
xfdcell = {xfd};
bwtcell = {cfdPar};

[bwtcellout, awtcellout, rescell] = pda_fd(xfdcell, bwtcell);

%  display weight coefficient for variable
bwtfd      = getfd(bwtcellout{1});
subplot(1,1,1)
plot(bwtfd)
title('Weight coefficient for variable')
disp(getcoef(bwtfd))
%  display residual functions
plot(rescell{1})
title('Residual function')

%  Example 2:  a single first order varying coefficient unforced equation
%     Dx(t) = -t.*x(t) or x(t) = exp(-t^2/2)

bvec    = tvec;
xvec    = exp(-tvec.^2./2);
xfd     = smooth_basis(tvec, xvec, xbasis);
xfdcell = {xfd};
bwtcell = {mfdPar};

[bwtcellout, awtcellout, rescell] = pda_fd(xfdcell, bwtcell);

%  display weight coefficient for variable
bwtfd      = getfd(bwtcellout{1});
subplot(1,1,1)
plot(bwtfd)
title('Weight coefficient for variable')
disp(getcoef(bwtfd))
%  display residual function
plot(rescell{1})
title('Residual function')

%  Example 3:  a single second order constant coefficient unforced equation
%     Dx(t) = -(2.*pi).^2.*x(t) or x(t) = sin(2.*pi.*t)
%  (2*pi)^2 = 39.4784
xvec    = sin(2.*pi.*tvec);
xfd     = smooth_basis(tvec, xvec, xbasis);
xfdcell = {xfd};
bwtcell = {cfdPar,cfdPar};

[bwtcellout, awtcellout, rescell] = pda_fd(xfdcell, bwtcell);

%  display weight coefficients
bwtfd1     = getfd(bwtcellout{1});
bwtfd2     = getfd(bwtcellout{2});
subplot(2,1,1)
plot(bwtfd1)
title('Weight coefficient for variable')
subplot(2,1,2)
plot(bwtfd2)
title('Weight coefficient for derivative of variable')
disp([getcoef(bwtfd1);getcoef(bwtfd2)]);
%  display residual function
subplot(1,1,1)
plot(rescell{1})
title('Residual function')

%  Example 4:  two first order constant coefficient unforced equations
%     Dx1(t) = x2(t) and Dx2(t) = -x1(t)
%   equivalent to  x1(t) = sin(2.*pi.*t)

xvec1     = sin(2.*pi.*tvec);
xvec2     = 2.*pi.*cos(2.*pi.*tvec);
xfd1      = smooth_basis(tvec, xvec1, xbasis);
xfd2      = smooth_basis(tvec, xvec2, xbasis);
xfdcell   = {xfd1;xfd2};
bwtcell   = cell(2,2,1);
bwtcell{1,1,1} = cfdPar;
bwtcell{1,2,1} = cfdPar;
bwtcell{2,1,1} = cfdPar;
bwtcell{2,2,1} = cfdPar;

[bwtcellout, awtcellout, rescell] = pda_fd(xfdcell, bwtcell);

%  display weight coefficients
bwtfd11    = getfd(bwtcellout{1,1});
bwtfd21    = getfd(bwtcellout{2,1});
bwtfd12    = getfd(bwtcellout{1,2});
bwtfd22    = getfd(bwtcellout{2,2});
subplot(2,2,1)
plot(bwtfd11)
title('Weight for variable 1 in equation 1')
subplot(2,2,2)
plot(bwtfd12)
title('Weight for variable 1 in equation 2')
subplot(2,2,3)
plot(bwtfd21)
title('Weight for variable 2 in equation 1')
subplot(2,2,4)
plot(bwtfd22)
title('Weight for variable 2 in equation 2')
disp(getcoef(bwtfd11))
disp(getcoef(bwtfd12))
disp(getcoef(bwtfd21))
disp(getcoef(bwtfd22))
%  display residual functions
figure(1)
plot(rescell{1})
title('Residual function for variable 1')
figure(2)
plot(rescell{2})
title('Residual function for variable 2')

%  Example 5:  a single first order constant coefficient equation
%     Dx = -4.*x  for  x(t) = exp(-4t) forced by u(t) = 2

beta    = 4;
alpha   = 2;
xvec0   = exp(-beta.*tvec);
intv    = (exp(beta.*tvec) - 1)./beta;
xvec    = xvec0.*(1 + alpha.*intv);
xfd     = smooth_basis(tvec, xvec, xbasis);
xfdcell = {xfd};
bwtcell = {cfdPar};
awtcell = {cfdPar};
ufdcell = {fd(1,cbasis)};

[bwtcellout, awtcellout, rescell] = ...
                  pda_fd(xfdcell, bwtcell, awtcell, ufdcell);

%  display weight coefficients
bwtfd      = getfd(bwtcellout{1});
awtfd      = getfd(awtcellout{1});
figure(1)
plot(bwtfd)
title('Weight for variable')
figure(2)
plot(awtfd)
title('Weight for forcing function')
%  display residual function
figure(3)
plot(rescell{1})
title('Residual function')

%  Example 6:  two first order equations forced by constant functions
%     Dx = -4.*x     for  x(t) = exp(-4t)       forced by u(t) =  2
%     Dx = -4.*t.*x  for  x(t) = exp(-4t.^2./2) forced by u(t) = -1

beta    = 4;
xvec10  = exp(-beta.*tvec);
alpha1  = 2;
alpha2  = -1;
xvec1   = xvec10 + alpha1.*(1-xvec10)./beta;
xvec20  = exp(-beta.*tvec.^2./2);
vvec    = exp(beta.*tvec.^2./2);
intv    = 0.01.*(cumsum(vvec) - 0.5.*vvec);
xvec2   = xvec20.*(1 + alpha2.*intv);
xfd1    = smooth_basis(tvec, xvec1, xbasis);
xfd2    = smooth_basis(tvec, xvec2, xbasis);
xfdcell = {xfd1; xfd2};
bwtcell = cell(2,2);
bwtcell{1,1} = cfdPar;
bwtcell{1,2} = cfdPar;
bwtcell{2,1} = cfdPar;
bwtcell{2,2} = mfdPar;
awtcell = {cfdPar; cfdPar};
ufdcell = {fd(1,cbasis); fd(1,cbasis)};

[bwtcellout, awtcellout, rescell] = ...
                     pda_fd(xfdcell, bwtcell, awtcell, ufdcell);

% display weight functions for variables
bwtfd11    = getfd(bwtcellout{1,1});
bwtfd21    = getfd(bwtcellout{2,1});
bwtfd12    = getfd(bwtcellout{1,2});
bwtfd22    = getfd(bwtcellout{2,2});
figure(1)
plot(bwtfd11)
title('weight on variable 1 in equation 1')
figure(2)
plot(bwtfd21)
title('weight on variable 2 in equation 1')
figure(3)
plot(bwtfd12)
title('weight on variable 1 in equation 2')
figure(4)
plot(bwtfd22)
title('weight on variable 2 in equation 2')
disp(getcoef(bwtfd11))
disp(getcoef(bwtfd21))
disp(getcoef(bwtfd12))
disp(getcoef(bwtfd22))
%  display weight functions for forcing functions
awtfd1     = awtcellout{1};
awtfd2     = awtcellout{2};
figure(1)
plot(awtfd1)
title('weight on forcing function in equation 1')
figure(2)
plot(awtfd2)
title('weight on forcing function in equation 2')
%  display residual functions
figure(1)
plot(rescell{1})
title('residual function for equation 1')
figure(2)
plot(rescell{2})
title('residual function for equation 2')



