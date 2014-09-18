function [res, Dres, SSE, DSSE, PENSSE, fdobj, df, gcv] = ...
    profPDA_SISO2(bvec, fitcell, gradwrd)
% profPDA estimates a homogeneous linear differential equation
%  directly from discrete noisy observations of a process.  
%  In this version forcing functions can be accommodated, but not
%  multiple variables. 
%
%profPDA works with the basis function expansions of the
%  estimates of the coefficient functions a_k(t) and b_j(t) 
%  in the possibly nonhomogeneous linear differential operator
%
%    Lx(t) = a_1(t)u_1(t) + ... + a_k(t)u_K(t) + 
%       b_0(t)x(t) + ... + b_{m-1}(t)D^{m-1}x(t) + 
%       \exp(b_m(t))D^m x(t)
%
%  of order m = NORDER that minimizes in a least squares sense the residual
%  functions f(t) = Lx(t).  
%
%  Arguments:
%  BVEC       ... A super-vector containing vectors of coefficients defining 
%                 the weight functions to be estimated.  
%                 Coefficients are stacked on top of each other as follows:
%                 coefficients for b_0 first, b_1 next, and so on, and
%                 then continues on to the forcing function coefficients
%                 to be estimated.  See BVEC2LFD for more details.
%  FITCELL    ... A cell array of length 1.    
%                 FITCELL{1} contains a struct object that
%                 contains the following information:
%       Y          ... Matrix of curve values to be fit.
%       BASISOBJ   ... A basis object for representing the curve(s)
%                 that are the solutions of the differential equations
%                 that fit the data
%       BASISMAT   ... Matrix of basis matrix values.
%       BMAT       ... Weighted cross-product matrix for basis matrix values.
%       DMAT       ... Weighted product of basis matrix and Y.
%       BWTCELL    ... Cell object for the weight functions for the 
%                 homogeneous part of the equation.
%       AWTCELL    ... Cell object for the weight functions for the 
%                 forcing functions
%       UFDCELL    ... Cell object containing functional data objects for
%                 the forcing functions
%       LAMBDA     ... The smoothing parameter controlling the penalty on
%                 the roughness penalty.  In order to estimate the 
%                 DIFE, this should be large but not too large.  
%                 If it is too small, the roughness is ignored, and
%                 the differential equation is badly estimated.
%                 If it is too large, though, the data are ignored and
%                 the DIFE is also badly estimated.  
%       DORDER     ... The order of the differential equation
%  GRADWRD    ... Returns the gradient in DSSE if nonzero, and []
%                 otherwise.
%  Returns:
%  SSE     ...  The error sum of squares.  This is what is
%               required for numerical optimization software since it
%               is the objective function being minimized.
%  DSSE    ...  The gradient of SSE.  Computed only if DERIVS is nonzero,
%               otherwise is returned as an empty matrix.
%  PENSSE  ...  The penalized error sum of squares.  
%               It is PENSSE = SSE + lambda.*C'KC, where C is the 
%               coefficient vector or matrix defining the curve(s) 
%               fitting the data and K is the penalty matrix corresponding
%               to the estimated DIFE.                 
%  FDOBJ   ...  Functional data object fitting the data
%  DF      ...  A measure of the equivalent degrees of freedom in the
%               smoothing matrix.
%  GCV     ...  The generalized cross-validation criterion

%  Last modified 30 March 2006

if nargin < 3
    gradwrd = 1;
end

%  check number of variables

fitstruct = fitcell{1};
bwtcell  = fitstruct.bwtcell;
awtcell  = fitstruct.awtcell;
ufdcell  = fitstruct.ufdcell;
y        = fitstruct.y;
basisobj = fitstruct.basisobj;
basismat = fitstruct.basismat;
Bmat     = fitstruct.Bmat;
Dmat     = fitstruct.Dmat;
lambda   = fitstruct.lambda;

J = size(bwtcell,1);

if J > 1 & length(bwtcell) > 1
    error('This version cannot handle multiple variables.');
end

%  get number of data points and number of curves

[n, ncurves] = size(y);
if n == 1 & ncurves > 1
    %  transpose Y if n == 1
    y = y';
    n = ncurves;
    ncurves = 1;
end

%  check basis

if ~isa_basis(basisobj)
    error('BASISOBJ is not a basis object.');
end

nbasis   = getnbasis(basisobj);
onebasis = ones(1,nbasis);
ncurves  = size(Dmat,2);

%  check LAMBDA

if lambda < 0
    warning ('Value of LAMBDA was negative, and 0 used instead.');
    lambda = 0;
end

%  transfer parameters from vector BVEC to cell objects
%  BWTCELL and AWTCELL

nderiv = size(bwtcell,length(size(bwtcell)));

%  calculate number of forcing functions

if isempty(awtcell) | isempty(ufdcell)
    nforce = 0;
    awtcell = {};
else
    nforce = size(ufdcell,2);
end

%  loop through derivatives to transfer coefficients from
%  BVEC to the appropriate derivative weight functions

[fitcell, npar] = bvec2fitcell(bvec, fitcell);
fitstruct = fitcell{1};
bwtcell  = fitstruct.bwtcell;
awtcell  = fitstruct.awtcell;
    
%  compute the penalty matrix and penalty vector

[penmat, penvec, DR, Ds] = ...
    eval_Rs(bwtcell, awtcell, ufdcell, basisobj, gradwrd);

%  check for ill conditioning in the linear equations

Bnorm   = sqrt(sum(sum(Bmat.^2)));
pennorm = sqrt(sum(sum(penmat.^2)));
condno  = pennorm/Bnorm;
if lambda*condno > 1e12
    lambda = 1e12/condno;
    warning(['LAMBDA reduced to ',num2str(lambda),...
            ' to prevent overflow']);
end

%  update the coefficient matrix with the penalty matrix

Mmat  = Bmat + lambda .* penmat;

%  update the right side vector with the penalty vector

if ~isempty(awtcell) & ~isempty(ufdcell)
    %  if the linear differential operator is nonhomogeneous
    %  use PENVEC to alter the right side of the equation.
    Dmat = Dmat + lambda.*(penvec*ones(1,ncurves));
end

%  compute inverse of Mmat

if is_diag(Mmat)
    Mmatinv = diag(1./diag(Mmat));
else
    Mmatinv = inv(Mmat);
end

%  compute degrees of freedom of smooth

df = sum(diag(Mmatinv*Bmat));

%  solve normal equations for each observation

coef = Mmatinv * Dmat;

%  compute error sum of squares

yhat = basismat * coef;
res  = y - yhat;
SSE  = sum(sum(res.^2));

%  compute gradient

if gradwrd
    nparam = length(bvec);
    Dres = zeros(n, nparam);
    DSSE = zeros(nparam,1);
    PhiMmatinv = basismat * Mmatinv;
    m2 = 0;
    for j=1:nderiv
        bfdParj = bwtcell{j};
        if getestimate(bfdParj) 
            basisj  = getbasis(getfd(bfdParj));
            nbasisj = getnbasis(basisj);
            m1 = m2 + 1;
            m2 = m2 + nbasisj;
            for m = m1:m2
                DRm = squeeze(DR(:,:,m));
                if nforce > 0
                    DAm = lambda.*PhiMmatinv* ...
                        (DRm*coef - Ds(:,m)*ones(1,ncurves));
                else
                    DAm = lambda.*PhiMmatinv*DRm*coef;
                end
                Dres(:,m) = 2.*DAm;
                DSSE(m)   = 2.*(res'*DAm);
            end
        end
    end
    
    %  update penalized least squares for terms for 
    %  smoothing forcing function coefficients
    
    for k=1:nforce
        afdPark = awtcell{k};
        if getestimate(afdPark) 
            basisk  = getbasis(getfd(afdPark));
            nbasisk = getnbasis(basisk);
            m1 = m2 + 1;
            m2 = m2 + nbasisk;
            for m = m1:m2
                DAm = lambda.*PhiMmatinv*Ds(:,m);
                Dres(:,m) = -2.*DAm;
                DSSE = DSSE -2.*(res'*DAm);
            end
        end
    end
else
    Dres = [];
end

%  compute  GCV index

if df < n
    gcv = (SSE/n)/((n - df)/n)^2;
else
    gcv = NaN;
end

fdobj = fd(coef, basisobj);

%  penalized least squares

PENSSE = SSE + lambda.*sum(diag(coef'*penmat*coef));

%  update least squares SSE for terms for 
%  smoothing derivative weight coefficients

% m2 = 0;
% for j=1:nderiv
%     bfdParj = bwtcell{j};
%     if getestimate(bfdParj) 
%         basisj  = getbasis(getfd(bfdParj));
%         nbasisj = getnbasis(basisj);
%         m1 = m2 + 1;
%         m2 = m2 + nbasisj;
%         lambdaj = getlambda(bfdParj);
%         if lambdaj > 0
%             Lfdobjj = getLfd(bfdParj);
%             penmatj = eval_penalty(basisj, Lfdobjj);
%             bcoefj  = bvec(m1:m2);
%             termj   = lambdaj.*bcoefj'*penmatj*bcoefj;
%             SSE     = SSE + termj;
%         end
%     end
% end

%  update penalized least squares for terms for 
%  smoothing forcing function coefficients

% nforce = length(ufdcell);
% for k=1:nforce
%     afdPark = awtcell{k};
%     if getestimate(afdPark) 
%         basisk  = getbasis(getfd(afdPark));
%         nbasisk = getnbasis(basisk);
%         m1 = m2 + 1;
%         m2 = m2 + nbasisk;
%         lambdak = getlambda(afdPark);
%         if lambdak > 0
%             Lfdobjk = getLfd(afdPark);
%             penmatk = eval_penalty(basisk, Lfdobjk);
%             bcoefk  = bvec(m1:m2);
%             SSE     = SSE + lambdak.*bcoefk'*penmatk*bcoefk;
%         end
%     end
% end

%  --------------------------------------------------------------
%  update derivatives of least squares SSE for terms for 
%  smoothing derivative weight coefficients
%  --------------------------------------------------------------

% if gradwrd
%     
%     m2 = 0;
%     for j=1:nderiv
%         bfdParj = bwtcell{j};
%         if getestimate(bfdParj) 
%             basisj  = getbasis(getfd(bfdParj));
%             nbasisj = getnbasis(basisj);
%             m1 = m2 + 1;
%             m2 = m2 + nbasisj;
%             lambdaj = getlambda(bfdParj);
%             if lambdaj > 0
%                 Lfdobjj = getLfd(bfdParj);
%                 penmatj = eval_penalty(basisj, Lfdobjj);
%                 bcoefj  = bvec(m1:m2);
%                 termj   = lambdaj.*penmatj*bcoefj;
%                 DSSE(m1:m2) = DSSE(m1:m2) + 2*termj;
%             end
%         end
%     end
%     
%     %  update penalized least squares for terms for 
%     %  smoothing forcing function coefficients
%     
%     nforce = length(ufdcell);
%     for k=1:nforce
%         afdPark = awtcell{k};
%         if getestimate(afdPark) 
%             basisk  = getbasis(getfd(afdPark));
%             nbasisk = getnbasis(basisk);
%             m1 = m2 + 1;
%             m2 = m2 + nbasisk;
%             lambdak = getlambda(afdPark);
%             if lambdak > 0
%                 Lfdobjk = getLfd(afdPark);
%                 penmatk = eval_penalty(basisk, Lfdobjk);
%                 bcoefk  = bvec(m1:m2);
%                 DSSE(m1:m2) = DSSE(m1:m2) + 2*lambdak.*penmatk*bcoefk;
%             end
%         end
%     end    
% end
% 
%  ----------------------------------------------------------------

function temp = eval_Lphi(jvar, bwtcell, basisobj)
%  Evaluates the linear combination of basis function derivatives 
%  for variable JVAR in differential operator L
%  defined by weight cell array BWTCELL.
%  The evaluation is at the quadrature points.

%  Last modified 30 March 2006

Dorder    = size(bwtcell, 2);
nbasis    = getnbasis(basisobj);
onebas    = ones(1,nbasis);

%  the function term in the operator

bfd       = getfd(bwtcell{jvar,1});
bbasis    = getbasis(bfd);
bcoef     = getcoef(bfd);
bbasismat = getvalues(bbasis);
bmat      = (bbasismat*bcoef)*onebas;
basismat  = getvalues(basisobj);
temp      = basismat.*bmat;

%  the terms for the derivatives

for jderiv = 2:Dorder
    Dbfd       = getfd(bwtcell{jvar, jderiv});
    Dbbasis    = getbasis(Dbfd);
    Dbcoef     = getcoef(Dbfd);
    Dbbasismat = getvalues(Dbbasis);
    Dbmat      = (Dbbasismat*Dbcoef)*onebas;
    Dbasismat  = getvalues(basisobj,jderiv-1);
    temp       = temp + Dbasismat.*Dbmat;
end

%  --------------------------------------------------------------

function auvec = eval_au(quadvals, awtcell, ufd)
%  evaluate weight function times forcing function
afd       = getfd(awtcell);
abasis    = getbasis(afd);
acoef     = getcoef(afd);
abasismat = getvalues(abasis);
avec      = abasismat*acoef;
uvec      = eval_fd(quadvals(:,1), ufd);
auvec     = avec.*uvec.*sqrt(quadvals(:,2));

