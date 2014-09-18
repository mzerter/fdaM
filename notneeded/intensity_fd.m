function [Wfdobj, hmat, Fstr, iter, iterhist] = ...
    intensity_fd(x, Wfd0Par, conv, iterlim, dbglev)
%INTENSITY_FD estimates the intensity function \lambda(x) of a 
%  nonhomogeneous Poisson process from a sample of event times.

%  Arguments are:
%  X       ... data value array.
%  WFDPAR  ... functional parameter object specifying the initial log
%              density, the linear differential operator used to smooth
%              smooth it, and the smoothing parameter.
%  CONV    ... convergence criterion
%  ITERLIM ... iteration limit for scoring iterations
%  DBGLEV  ... level of output of computation history

%  Returns:
%  WFDOBJ   ... functional data basis object defining final 
%               log intensity.
%  FSTR     ... Struct object containing
%               FSTR.F    ... final log likelihood
%               FSTR.NORM ... final norm of gradient
%  ITERNUM  ... Number of iterations
%  ITERHIST ... History of iterations

%  last modified 20 July 2006

if nargin < 2
    error('WFDPAR is not supplied.');
end

%  check Wfd0Par

if ~isa_fdPar(Wfd0Par) 
    if isa_fd(Wfd0Par) || isa_basis(Wfd0Par)
        Wfd0Par = fdPar(Wfd0Par);
    else
        error(['WFDPAR is not a functional parameter object, ', ...
                'not a functional data object, and ', ...
                'not a basis object.']);
    end
end

%  set up WFDOBJ

Wfdobj = getfd(Wfd0Par);

%  set up LFDOBJ

Lfdobj = getLfd(Wfd0Par);
Lfdobj = int2Lfd(Lfdobj);

%  set up LAMBDA

lambda = getlambda(Wfd0Par);

%  set up BASISOBJ

basisobj = getbasis(Wfdobj);
nbasis   = getnbasis(basisobj);
rangex   = getbasisrange(basisobj);
active   = 1:nbasis;

%  check X for being a vector

[N,m] = size(x);

if N == 1 && m > 1
    x = x';
    N = m;
elseif m > 1 && N > 1
    error('Argument X is not a vector.');
end

%  check for values outside of the range of WFD0

inrng = find(x >= rangex(1) & x <= rangex(2));
if (length(inrng) ~= N)
    disp([length(inrng), N])
    disp([rangex(1),rangex(2),min(x),max(x)])
    warning('Some values in X out of range and not used.');
end

x    = x(inrng);

%  set some default arguments and constants

if nargin < 5, dbglev  = 1;        end
if nargin < 4, iterlim = 20;       end
if nargin < 3, conv    = 1e-4;     end

%  set up some arrays

climit = [-50,0;0,400]*ones(2,nbasis);
cvec0  = getcoef(Wfdobj);
dbgwrd = dbglev > 1;

%  initialize matrix Kmat defining penalty term

if lambda > 0
    Kmat = lambda.*eval_penalty(basisobj, Lfdobj);
end

%  evaluate log likelihood
%    and its derivatives with respect to these coefficients

[logl, Dlogl] = loglfninten(x, basisobj, cvec0);

%  compute initial badness of fit measures

Foldstr.f  = -logl;
gvec       = -Dlogl;
if lambda > 0
    gvec      = gvec           + 2.*(Kmat * cvec0);
    Foldstr.f = Foldstr.f + cvec0' * Kmat * cvec0;
end
Foldstr.norm = sqrt(mean(gvec.^2));

%  compute the initial expected Hessian

hmat = Varfninten(basisobj, cvec0);
if lambda > 0
    hmat = hmat + 2.*Kmat;
end

%  evaluate the initial update vector for correcting the initial bmat

deltac   = -hmat\gvec;

%  initialize iteration status arrays

iternum = 0;
status = [iternum, Foldstr.f, -logl, Foldstr.norm];
if dbglev > 0
    fprintf('\nIteration  Criterion  Neg. Log L  Grad. Norm\n')
    fprintf('\n%5.f     %10.4f %10.4f %10.4f\n', status);
end
iterhist = zeros(iterlim+1,length(status));
iterhist(1,:)  = status;

%  quit if ITERLIM == 0

if iterlim == 0
    Fstr = Foldstr;
    iterhist = iterhist(1,:);
    return;
end

%  -------  Begin iterations  -----------

STEPMAX = 5;
MAXSTEP = 400;
trial   = 1;
cvec    = cvec0;
linemat = zeros(3,5);

for iter = 1:iterlim
    iternum = iternum + 1;
    %  take optimal stepsize
    dblwrd = [0,0]; limwrd = [0,0]; ind = 0;
    %  compute slope
    Fstr = Foldstr;
    linemat(2,1) = sum(deltac.*gvec);
    %  normalize search direction vector
    sdg     = sqrt(sum(deltac.^2));
    deltac  = deltac./sdg;
    linemat(2,1) = linemat(2,1)/sdg;
    %  return with error condition if initial slope is nonnegative
    if linemat(2,1) >= 0
        disp('Initial slope nonnegative.')
        iterhist = iterhist(1:(iternum+1),:);
        break;
    end
    %  return successfully if initial slope is very small
    if linemat(2,1) >= -1e-5;
        if dbglev>1, disp('Initial slope too small'); end
        iterhist = iterhist(1:(iternum+1),:);
        break;
    end
    %  load up initial search matrix 
    linemat(1,1:4) = 0;
    linemat(2,1:4) = linemat(2,1);
    linemat(3,1:4) = Foldstr.f;
    %  output initial results for stepsize 0
    stepiter  = 0;
    if dbglev>1
        fprintf('      %3.f %10.4f %10.4f %10.4f\n', ...
            [stepiter, linemat(:,1)']);
    end
    ips = 0;
    %  first step set to trial
    linemat(1,5)  = trial;
    %  Main iteration loop for linesrch
    for stepiter = 1:STEPMAX
        %  check the step size
        [linemat(1,5),ind,limwrd] = ...
            stepchk(linemat(1,5), cvec, deltac, limwrd, ind, ...
            climit, active, dbgwrd);
        if linemat(1,5) <= 1e-9
            %  Current step size too small ... terminate
            Fstr    = Foldstr;
            cvecnew = cvec;
            gvecnew = gvec;
            if dbglev > 1
                fprintf('Stepsize too small:  %10.4f\n', linemat(1,5));
            end
            break;
        end
        cvecnew = cvec + linemat(1,5).*deltac;
        %  compute new function value and gradient
        [logl, Dlogl] = loglfninten(x, basisobj, cvecnew);
        Fstr.f  = -logl;
        gvecnew = -Dlogl;
        if lambda > 0
            gvecnew = gvecnew + 2.*Kmat * cvecnew;
            Fstr.f = Fstr.f + cvecnew' * Kmat * cvecnew;
        end
        Fstr.norm = sqrt(mean(gvecnew.^2));
        linemat(3,5) = Fstr.f;
        %  compute new directional derivative
        linemat(2,5) = sum(deltac.*gvecnew);
        linemat(3,5) = Fstr.f;
        %  output current results
        if dbglev > 1
            fprintf('      %3.f %10.4f %10.4f %10.4f\n', ...
                [stepiter, linemat(:,5)']);
        end
        %  compute next step
        [linemat,ips,ind,dblwrd] = ...
            stepit(linemat, ips, ind, dblwrd, MAXSTEP, dbgwrd);
        trial  = linemat(1,5);
        %  ind == 0 implies convergence
        if ind == 0 || ind == 5, break; end
        %  end iteration loop
    end
    
    %  update current parameter vectors
    
    cvec   = cvecnew;
    gvec   = gvecnew;
    Wfdobj = putcoef(Wfdobj, cvec);
    status = [iternum, Fstr.f, -logl, Fstr.norm];
    iterhist(iter+1,:) = status;
    if dbglev > 0
        fprintf('%5.f     %10.4f %10.4f %10.4f\n', status);
    end
    %  test for convergence
    if abs(Fstr.f-Foldstr.f) < conv
        iterhist = iterhist(1:(iternum+1),:);
        break;
    end
    %  exit loop if convergence
    if Fstr.f >= Foldstr.f, break; end
    %  compute the new Hessian
    hmat = Varfninten(basisobj, cvec);
    if lambda > 0
        hmat = hmat + 2.*Kmat;
    end
    %  evaluate the update vector
    deltac    = -hmat\gvec;
    cosangle  = -gvec'*deltac/sqrt(sum(gvec.^2)*sum(deltac.^2));
    if cosangle < 0
        if dbglev > 1, disp('cos(angle) negative'); end
        deltac = -gvec;
    end
    Foldstr = Fstr;
end

%  ---------------------------------------------------------------

function [logl, Dlogl] = loglfninten(x, basisobj, cvec)
Cval    = normalize_phi(basisobj, cvec);
phimat  = full(eval_basis(x, basisobj));
logl    = sum(phimat * cvec) - Cval;
EDW     = expect_phi(basisobj, cvec);
Dlogl   = (sum(phimat)' - EDW);

%  ---------------------------------------------------------------

function  Varphi = Varfninten(basisobj, cvec)
Varphi  = expect_phiphit(basisobj, cvec);

%  ---------------------------------------------------------------

function Cval = normalize_phi(basisobj, cvec)

%  Computes integral from RNG(1) to x(nobs) of
%      \exp [phi'(x) * cvec]
%  by numerical integration using Romberg integration

%  Arguments:
%  BASISOBJ   ...  Basis function object with basis functions phi.
%  CVEC       ... coefficient vector defining density, of length NBASIS

%  Return:
%  The integral of the function.

%  check arguments, and convert basis objects to functional data objects

if ~strcmp(class(basisobj), 'basis')
    error('First argument must be a basis function object.');
end

Cval = funcint(@sumpxfun, cvec, basisobj);

%  ---------------------------------------------------------------

function Ephi = expect_phi(basisobj, cvec)
%  Computes inner products of basis functions with 
%      p(x) = exp [c'phi(x)],
%  by numerical integration using Romberg integration

%  Arguments:
%  BASISOBJ  ...  A basis object
%           object.  In the latter case, a functional data object
%           is created from a basis function object by using the
%           identity matrix as the coefficient matrix.
%           The functional data objects must be univariate.
%  CVEC ... coefficient vector defining density, of length NBASIS
%  JMAX ... maximum number of allowable iterations
%  EPS  ... convergence criterion for relative error

%  Return:
%  A vector SS of length NBASIS of inner products.

%  check arguments, and convert basis objects to functional data objects

if ~strcmp(class(basisobj),'basis')
    error('First argument must be a basis function object.');
end

Ephi = funcint(@fxpxfun, cvec, basisobj);

%  ---------------------------------------------------------------

function Ephiphi = expect_phiphit(basisobj, cvec)

%  Computes expectations of cross product of basis functions with
%  respect to density
%      p(x) = exp c'phi(x)
%  by numerical integration using Romberg integration

%  Arguments:
%  BASISOBJ  ...  A basis object
%           object.  In the latter case, a functional data object
%           is created from a basis function object by using the
%           identity matrix as the coefficient matrix.
%           The functional data objects must be univariate.
%  CVEC ... coefficient vector defining density

%  Return:
%  A matrix of order NBASIS of integrals of functions.

%  check arguments, and convert basis objects to functional data objects

if ~strcmp(class(basisobj),'basis')
    error('First argument must be a basis function object.');
end

Ephiphi = funcint(@fxpxfxfun, cvec, basisobj);

%  ---------------------------------------------------------------

function sumpx = sumpxfun(x, cvec, basisobj, nderiv1, nderiv2)
fx = full(eval_basis(x, basisobj));
wx = fx * cvec;
wx(wx < -50) = -50;
px = exp(wx);
sumpx = sum(px);

%  ---------------------------------------------------------------

function fxpx = fxpxfun(x, cvec, basisobj, nderiv1, nderiv2)
fx = full(eval_basis(x, basisobj));
wx = fx * cvec;
wx(wx < -50) = -50;
px   = exp(wx);
fxpx = fx' * px;

%  ---------------------------------------------------------------

function fxpxfx = fxpxfxfun(x, cvec, basisobj, nderiv1, nderiv2)
nbasis = getnbasis(basisobj);
fx = full(eval_basis(x, basisobj));
wx = fx * cvec;
wx(wx < -50) = -50;
px   = exp(wx);
fxpxfx = fx'*((px*ones(1,nbasis)).*fx);

