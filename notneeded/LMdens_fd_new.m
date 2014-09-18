function [Wfdobj, beta, Cval, Fstr, iternum, iterhist] = ...
    LMdens_fd_new(y, WfdPar, zmat, beta0, sigma0, ...
                  conv, iterlim, active, dbglev)
%LMDENS_FD estimates a regression and the density of the residuals.
%  If betaO is [], or empty, then only the density is estimated. 

%  Arguments are:
%  Y       ... vector of dependent variable values
%  WFDPAR  ... functional parameter object specifying the initial log
%              density, the linear differential operator used to smooth
%              smooth it, and the smoothing parameter.
%  ZMAT    ... matrix of covariates
%  BETA0   ... initial model   coefficients
%  SIGMA0  ... initial standard error
%  CONV    ... convergence criterion
%  ITERLIM ... iteration limit for scoring iterations
%  ACTIVE  ... indices among 1:NBASIS of parameters to optimize 
%              Normally the first index is set inactive.
%  DBGLEV  ... level of output of computation history

%  Returns:
%  WFD      ... functional data basis object defining final density
%  BETA     ... final model coefficients
%  FSTR     ... Struct object containing
%               FSTR.F    ... final log likelihood
%               FSTR.NORM ... final norm of gradient
%  ITERNUM  ... Number of iterations
%  ITERHIST ... History of iterations

%  last modified 23 October 2005

if nargin < 2
    error('WFDPAR is not supplied.');
end

%  check WfdPar

if ~isa_fdPar(WfdPar) 
    if isa_fd(WfdPar) | isa_basis(WfdPar)
        WfdPar = fdPar(WfdPar);
    else
        error(['WFDPAR is not a functional parameter object, ', ...
                'not a functional data object, and ', ...
                'not a basis object.']);
    end
end

%  set up WFDOBJ

Wfdobj = getfd(WfdPar);
cvec0  = getcoef(Wfdobj);

%  set up LFDOBJ

Lfdobj = getLfd(WfdPar);
Lfdobj = int2Lfd(Lfdobj);
nderiv = getnderiv(Lfdobj);

%  set up BASIS

basisobj = getbasis(Wfdobj);
nbasis   = getnbasis(basisobj);
rangey   = getbasisrange(basisobj);

%  get basis for centered coefficients

zerobasmat = zerobasis(nbasis);
czerovec0  = zerobasmat'*cvec0;

nobs   = length(y);

%  set some default arguments

if nargin < 9, dbglev  = 1;         end
if nargin < 8, active  = 1:nbasis;  end
if nargin < 7, iterlim = 20;        end
if nargin < 6, conv    = 1e-2;      end

%  initialize some arrays

climit    = [-50,0;0,50]*ones(2,nbasis);
hmat      = zeros(nbasis,nbasis);
inactive  = ones(1,nbasis);
inactive(active) = 0;
inactive  = find(inactive);
ninactive = length(inactive);
dbgwrd    = dbglev > 1;

ind1      = 1:(nbasis-1);

%  Set up linear model

covwrd = ~isempty(beta0);
if covwrd
    %  Independent variables are present ... estimate linear model
    if size(zmat,1) ~= nobs 
        error('ZMAT must have as many rows as length(X)');
    end
    ncov = size(zmat,2);
    Zsum = sum(zmat)';
    res0 = (y - zmat * beta0);
    ind2 = (nbasis):(nbasis+ncov-1);
    hessmap = [[zerobasmat, zeros(nbasis,ncov)]; ...
               [zeros(ncov,nbasis-1), eye(ncov)]];
else
    %  No independent variables are present
    ncov = 0;
    zmat = [];
    Zsum = [];
    res0 = y;
    hessmap = zerobasmat;
end

npar = nbasis + ncov;

%  bring residuals out of range to range and set 
%    U0 = res./sigma0;

[res0, U0, indlo, indhi] = reschk(res0, rangey, sigma0);

%  initialize matrix Kmat defining penalty term

lambda = getlambda(WfdPar);
if lambda > 0 
    Kmat = lambda.*eval_penalty(basisobj, Lfdobj);
    Kzeromat = zerobasmat'*Kmat*zerobasmat;
end

%  evaluate log likelihood
%    and its derivatives with respect to these coefficients

[logl, Dlogl, Cval, Ephi] = ...
    loglfnLM(basisobj, cvec0, U0, zmat, sigma0, zerobasmat);
% disp(logl)
% disp(Dlogl)

%  compute initial badness of fit measures

Foldstr.f  =  -logl;
gvec       = -Dlogl;
if lambda > 0 
    gvec(ind1) = gvec(ind1) + 2.*(Kmat * czerovec0);
    Foldstr.f  = Foldstr.f  + cvec0' * Kmat * cvec0;
end
if ninactive > 0, gvec(inactive) = 0;  end
Foldstr.norm = sqrt(mean(gvec.^2));

%  compute the initial expected Hessian

hmat = -EHessfnLM(basisobj, cvec0, U0, zmat, Cval, Ephi, ...
    sigma0, zerobasmat);

% disp(hmat)

if lambda > 0 
    hmat(ind1,ind1) = hmat(ind1,ind1) + 2.*Kzeromat;
end

% if ninactive > 0 
%     hmat(inactive,:) = 0;
%     hmat(:,inactive) = 0;
%     hmat(inactive,inactive) = eye(ninactive);
% end

%  evaluate the initial update vector for correcting the initial bmat

deltac    = -hmat\gvec;
cosangle  = -gvec'*deltac/sqrt(sum(gvec.^2)*sum(deltac.^2));
disp(deltac);

%  initialize iteration status arrays

iternum = 0;
status = [iternum, Foldstr.f, -logl, Foldstr.norm];
if dbglev > 0
    fprintf('\nIteration  Criterion  Neg. Log L  Grad. Norm\n')  
    fprintf('\n%5.f     %10.4f %10.4f %10.4f\n', status);
end
iterhist = zeros(iterlim+1,length(status));
iterhist(1,:) = status;

%  quit if ITERLIM == 0

if iterlim == 0
    Fstr = Foldstr;
    iterhist = iterhist(1,:);
    beta = beta0;
    if length(indlo) > 0
        warning([num2str(length(indlo)),' lower residuals trimmed']);
    end
    indhi = find(res > rangey(2)*sigma0);
    if length(indhi) > 0
        warning([num2str(length(indhi)),' upper residuals trimmed']);
    end
    return;
end

%  -------  Begin iterations  -----------

STEPMAX  = 5;
MAXSTEP  = 100;
trial    = 1;
cvec     = cvec0;
czerovec = czerovec0;
beta     = beta0;
linemat  = zeros(3,5);

for iter = 1:iterlim
    iternum = iternum + 1;
    Fstr    = Foldstr;
    %  set initial switches
    dblwrd = [0,0]; limwrd = [0,0]; stpwrd = 0; ind = 0; ips = 0;
    %  normalize search direction vector
    sdg     = sqrt(sum(deltac.^2));
    deltac  = deltac./sdg;
    %  compute initial slope
    linemat(2,1) = sum(deltac.*gvec);
    %  return with error condition if initial slope is nonnegative
    if linemat(2,1) >= 0
        fprintf('Initial slope nonnegative.\n');
        ind = 3;
        iterhist = iterhist(1:(iternum+1),:);
        break;
    end
    %  return successfully if initial slope is very small
    if linemat(2,1) >= -1e-5;
        if dbglev > 1, fprintf('Initial slope too small\n'); end
        iterhist = iterhist(1:(iternum+1),:);
        break;
    end
    %  load up initial search matrix 
    linemat(1,1:4) = 0;
    linemat(2,1:4) = linemat(2,1);
    linemat(3,1:4) = Foldstr.f;
    %  output initial results for stepsize 0
    stepiter  = 0;
    if dbglev > 1
        fprintf('      %3.f %10.4f %10.4f %10.4f\n', ...
            [stepiter, linemat(:,1)']); 
    end
    %  first step set to trial
    linemat(1,5)  = trial;
    %  Main iteration loop for linesrch
    for stepiter = 1:STEPMAX
        %  ensure that step does not go beyond limits on parameters
        %  check the step size
%         disp(deltac)
%         disp(deltac(ind1));
%         disp(deltac(ind2));
%         disp(zerobasmat)
        deltac0 = [zerobasmat*deltac(ind1); deltac(ind2)];
        [linemat(1,5),ind,limwrd] = ...
            stepchk(linemat(1,5), cvec, deltac0, limwrd, ind, ...
            climit, active, dbglev);
        if linemat(1,5) <= 1e-9 
            %  Current step size too small ... terminate
            Fstr    = Foldstr;
            cvecnew = cvec;
            betanew = beta;
            gvecnew = gvec;
            if dbglev > 1
                fprintf('Stepsize too small: %10.4f', linemat(1,5));
            end
            break;
        end
        cvecnew = cvec + linemat(1,5).*zerobasmat*deltac(ind1);
        czerovecnew = zerobasmat'*cvecnew;
        if covwrd
            betanew = beta + linemat(1,5).*deltac(ind2);
            resnew  = y - zmat * betanew;
        else
            resnew  = y;
        end
        %  compute new function value and gradient
        [resnew, Unew, indlo, indhi] = ...
              reschk(resnew, rangey, sigma0);
        [logl, Dlogl, Cval, Ephi]  = ...
              loglfnLM(basisobj, cvecnew, Unew, zmat, ...
              sigma0, zerobasmat);
        Fstr.f  =  -logl;
        gvecnew = -Dlogl;
        if lambda > 0 
            gvecnew(ind1) = gvecnew(ind1) + ...
                2.*Kzeromat * czerovecnew;
            Fstr.f = Fstr.f + cvecnew' * Kmat * cvecnew;
        end
%         if ninactive > 0, gvecnew(inactive) = 0;  end
        Fstr.norm  = sqrt(mean(gvecnew.^2));
        %  update search matrix
        linemat(2,5) = sum(deltac.*gvecnew);
        linemat(3,5) = Fstr.f;
        %  output current results
        if dbglev > 1 
            fprintf('      %3.f %10.4f %10.4f %10.4f\n', ...
                [stepiter, linemat(:,5)']); 
        end
        %  compute next step
        [linemat,ips,ind,dblwrd] = ...
            stepit(linemat, ips, ind, dblwrd, MAXSTEP, dbglev);
        trial  = linemat(1,5);
        %  ind == 0 implies convergence
        if ind == 0 | ind == 5, break; end
        %  end iteration loop
    end
    
    %  update current parameter vectors
    
    cvec     = cvecnew;
    czerovec = zerobasmat'*cvec;
    gvec     = gvecnew;
    Wfdobj   = putcoef(Wfdobj, cvec);
    if covwrd
        beta = betanew;
        res  = y - zmat * beta;
    else
        res = y;
    end
    %  check residuals and truncate if needed
    [res, U, indlo, indhi] = reschk(res, rangey, sigma0);
    %  update and output iteration status
    status = [iternum, Fstr.f, -logl, Fstr.norm];
    iterhist(iter+1,:) = status;
    fprintf('%5.f     %10.4f %10.4f %10.4f\n', status);
    %  test for convergence
    if abs(Fstr.f-Foldstr.f) < conv
        iterhist = iterhist(1:(iternum+1),:);
        if length(indlo) > 0
            warning([num2str(length(indlo)),' lower residuals trimmed']);
        end
        indhi = find(res > rangey(2)*sigma0);
        if length(indhi) > 0
            warning([num2str(length(indhi)),' upper residuals trimmed']);
        end
        break;
    end
    %  exit loop if convergence
    if Fstr.f >= Foldstr.f,  
        break;  
    end
    %  compute the new Hessian
    hmat = EHessfnLM(basisobj, cvec, U, zmat, Cval, Ephi, ...
        sigma0, zerobasmat);
    if lambda > 0
        hmat(ind1,ind1) = hmat(ind1,ind1) + 2.*Kzeromat;
    end
%     if ninactive > 0 
%         hmat(inactive,:) = 0;
%         hmat(:,inactive) = 0;
%         hmat(inactive,inactive) = eye(ninactive);
%     end
    %  evaluate the update vector
    deltac    = -hmat\gvec;
    cosangle  = -gvec'*deltac/ ...
                     sqrt(sum(gvec.^2)*sum(deltac.^2));
%     if cosangle < 0
%         if dbglev > 1, disp('cos(angle) negative');  end
%         deltac = -gvec;
%     end
    Foldstr = Fstr;    
end

%  ------------------------------------------------------

function [logl, Dlogl, Cval, Ephi] = ...
               loglfnLM(basisobj, cvec, U, zmat, sigma, zerobasmat)
%  U is a vector containing standardized residuals.
%  The first NBASIS elements in the gradient are the 
%  derivatives with respect to CVEC defining the density.
%  The next NCOV elements, if required, are the derivatives
%  with respect to BETA defining the linear model
nbasis = getnbasis(basisobj);
nobs   = length(U);
%  CVEC derivatives
phimat = getbasismatrix(U, basisobj);
Cval   = normalize_phi(basisobj, cvec);
logl   = sum(phimat*cvec  - log(Cval));
Ephi   = expect_phi(basisobj, cvec, Cval);
Dlogl  = sum(phimat)' - nobs.*Ephi;
Dlogl  = zerobasmat'*Dlogl;
%  BETA derivatives
if ~isempty(zmat)
    Dphimat = getbasismatrix(U, basisobj, 1);
    DWvec   = Dphimat*cvec;
    Dlogl   = [Dlogl; -zmat'*DWvec./sigma];
else
    Dphimat = [];
end

%  ------------------------------------------------------

function  EHess = EHessfnLM(basisobj, cvec, U, zmat, ...
                            Cval, Ephi, sigma, zerobasmat) 
%  The upper left order NBASIS sub-matrix contains
%  the Hessian with respect to CVEC defining the density.
%  The lower right order NCOV sub-matrix, if required, 
%  is the expected Hessian with respect to BETA defining 
%  the linear model.  The off-diagonal submatrices 
%  contain the cross derivatives
nbasis  = getnbasis(basisobj);
nobs    = size(zmat,1);
Ephiphi = expect_phiphit(basisobj, cvec, Cval);
EHess   = -nobs*(Ephiphi - Ephi*Ephi');
EHess   = zerobasmat'*EHess*zerobasmat;
if ~isempty(zmat)
    ncov  = size(zmat,2);
    onesn = ones(nobs,1);
%  This code for exact hessian
%     D1phimat = getbasismatrix(U, basisobj, 1);
%     D2phimat = getbasismatrix(U, basisobj, 2);
%     D2Wvec   = D2phimat * cvec;
%     DcdlnL   = -zmat'*D1phimat./sigma;
%     DddlnL   = (zmat.*(D2Wvec*ones(1,ncov)))'*zmat./sigma^2;
%     EHess    = [EHess, DcdlnL'; DcdlnL, DddlnL];
%  This code is for expected hessian
    EDphi   = expect_phi(basisobj, cvec, Cval, 1);
    EDcdlnL = -zmat'*onesn*EDphi'./sigma;
    EDcdlnL = EDcdlnL*zerobasmat;
    ED2phi  = expect_phi(basisobj, cvec, Cval, 2);
    EDddlnL = zmat'*zmat.*(ED2phi'*cvec)./sigma^2;
    EHess   = [EHess, EDcdlnL'; EDcdlnL, EDddlnL];
end

%  ------------------------------------------------------

function  [res, U, indlo, indhi] = ...
                              reschk(res, rangey, sigma0) 
%RESCHK brings residuals outside of limits to limits
indlo = find(res < rangey(1)*sigma0);
if length(indlo) > 0
    res(indlo) = rangey(1)*sigma0;
end
indhi = find(res > rangey(2)*sigma0);
if length(indhi) > 0
    res(indhi) = rangey(2)*sigma0;
end
U   = res./sigma0;

%  ---------------------------------------------------------------

function Cval = normalize_phi(basisobj, cvec)

%  Computes integrals of
%      p(x) = exp phi'(x) * cvec
%  by numerical integration using Romberg integration

%  Arguments:
%  basisobj  ...  Basis function object with basis functions phi.
%  CVEC ... coefficient vector defining density, of length NBASIS
%  MU   ... mean values to be subtracted from variates
%  SIGMA .. standard deviation to define u = (x - mu)/sigma
%  RNG ...  vector of length 2 giving the interval over which the
%           integration is to take place.  Multiply a standard interval
%           like (-5,5) by sigma to make it scale free
%  JMAX ... maximum number of allowable iterations
%  EPS  ... convergence criterion for relative error

%  Return:
%  The integral of the function.

%  check arguments, and convert basis objects to functional data objects

if ~strcmp(class(basisobj), 'basis')
    error('First argument must be a basis function object.');
end

Cval = funcint(@sumpxfun, cvec, basisobj, 0, 0);

%  ---------------------------------------------------------------

function Ephi = expect_phi(basisobj, cvec, Cval, nderiv)
%  Computes expectations of basis functions with respect to density
%      p(x) = Cval^{-1} exp c'phi(x)
%  by numerical integration using Romberg integration

%  Arguments:
%  basisobj  ...  A basis function
%           object.  In the latter case, a functional data object
%           is created from a basis function object by using the
%           identity matrix as the coefficient matrix.
%           The functional data objects must be univariate.
%  CVEC ... coefficient vector defining density, of length NBASIS
%  CVAL ... normalizing constant defining density
%  MU   ... mean value to be subtracted from variates
%  SIGMA .. standard deviation to define u = (x - mu)/sigma
%  RNG ...  vector of length 2 giving the interval over which the
%           integration is to take place
%  UWRD ... If T, expectation is of (D PHI)*U
%  JMAX ... maximum number of allowable iterations
%  EPS  ... convergence criterion for relative error

%  Return:
%  A vector SS of length NBASIS of integrals of functions.

%  check arguments, and convert basis objects to functional data objects

if ~strcmp(class(basisobj),'basis')
    error('First argument must be a basis function object.');
end

EPS  = 1e-7; 
JMAX = 15;   

if nargin < 4, nderiv = 0;  end

Ephi = funcint(@fxpxfun, cvec, basisobj, nderiv, 0)./Cval;

%  ---------------------------------------------------------------

function Ephiphi = expect_phiphit(basisobj, cvec, Cval)

%  Computes expectations of cross product of basis functions with
%  respect to density
%      p(x) = Cval^{-1} exp c'phi(x)
%  by numerical integration using Romberg integration

%  Arguments:
%  basisobj  ...  A basis function
%           object.  In the latter case, a functional data object
%           is created from a basis function object by using the
%           identity matrix as the coefficient matrix.
%           The functional data objects must be univariate.
%  CVEC ... coefficient vector defining density
%  CVAL ... normalizing constant defining density
%  RNG ...  vector of length 2 giving the interval over which the
%           integration is to take place
%  JMAX ... maximum number of allowable iterations
%  EPS  ... convergence criterion for relative error

%  Return:
%  A matrix of order NBASIS of integrals of functions.

%  check arguments, and convert basis objects to functional data objects

if ~strcmp(class(basisobj),'basis')
    error('First argument must be a basis function object.');
end

EPS  = 1e-7; 
JMAX = 9;  

Ephiphi = funcint(@fxpxfxfun, cvec, basisobj, 0, 0)./Cval;

%  ---------------------------------------------------------------

function sumpx = sumpxfun(x, cvec, basisobj, nderiv1, nderiv2)
if nargin < 4, nderiv1 = 0; end
if nargin < 5, nderiv2 = 0; end
fx = full(eval_basis(x, basisobj));
wx = fx * cvec;
wx(wx < -50) = -50;
px = exp(wx);
sumpx = sum(px);

%  ---------------------------------------------------------------

function fxpx = fxpxfun(x, cvec, basisobj, nderiv1, nderiv2)
if nargin < 4, nderiv1 = 0; end
if nargin < 5, nderiv2 = 0; end
fx = full(eval_basis(x, basisobj));
wx = fx * cvec;
wx(wx < -50) = -50;
px   = exp(wx);
if nderiv1 == 0
    Dfx = fx;
else
    Dfx = getbasismatrix(x, basisobj, 1);
end
fxpx = Dfx' * px;

%  ---------------------------------------------------------------

function fxpxfx = fxpxfxfun(x, cvec, basisobj, nderiv1, nderiv2)
if nargin < 4, nderiv1 = 0; end
if nargin < 5, nderiv2 = 0; end
nbasis = getnbasis(basisobj);
phimat = full(eval_basis(x, basisobj));
wx = phimat * cvec;
wx(wx < -50) = -50;
px = exp(wx);
fxpxfx = (phimat.*(px*ones(1,nbasis)))' * phimat;



