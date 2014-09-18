function [bsplinemat,P] = naturalbsplineM(x, breaks, norder, nderiv, sparsewrd)
%  NATURALBSPLINEM  Computes values or derivative values of natural
%  B-spline basis functions as well as projection matrix to convert
%  B-splines to natural B-splines
%  Arguments:
%  X         ... Argument values for which function values are computed
%  BREAKS    ... Increasing knot sequence spanning argument range
%  NORDER    ... Order of B-spline (one greater than degree) max = 19
%                Default 4.
%  NDERIV    ... Order of derivative required, default 0.
%  SPARSEWRD ... if 1, return in sparse form
%  Return:
%  BSPLINEMAT ... length(X) times number of basis functions matrix
%                 of Bspline values

%  last modified 13 August 2011 by Kris Villez based on bsplinepen in the
%  FDA toolbox

%  check dimensions of X and set up as a row vector

sizex = size(x);
if sizex(1) > 1 && sizex(2) > 1
    error('Argument X is not a vector.');
end
x = x(:);

n = length(x);

%  set default argument values

if nargin < 5, sparsewrd = 1;  end
if nargin < 4, nderiv    = 0;  end
if nargin < 3, norder    = 4;  end

ndcon       =   2 ; %norder-2    ; % constraint at derivative 2

% first compute B-splines second derivative coefficients at first and last knot.
bsplinemat  =   bsplineM(breaks, breaks, norder, ndcon, sparsewrd)   ;

IndexL  =   find(bsplinemat(1,:)~=0)        ;   % left side constraint
IndexR  =   find(bsplinemat(end,:)~=0)      ;   % right side constraint
IndexR  =   flipud(IndexR(:))               ;
CoeffL  =   bsplinemat(1,IndexL)            ;   % left side constraint coeff
CoeffR  =   bsplinemat(end,IndexR)          ;   % right side constraint coeff
nL      =   length(CoeffL)                  ;
nR      =   length(CoeffR)                  ;

% projection matrix
nb      =   length(breaks) + norder - 2     ;

PL      =   sparse(eye(nL))         ;
PL(1,:) =   -CoeffL(:)/CoeffL(1)    ;
PL      =   PL(:,2:end)             ;
PR      =   sparse(eye(nR))         ;
PR(1,:) =  -CoeffR(:)/CoeffR(1)     ;
PR      =   PR(:,2:end)             ;

P           =   zeros(nb,nb-2)  ;
vr          =   1:nL            ;
vc          =   1:nL-1          ;
P(vr,vc)    =   PL              ;
vr          =   nL+1:nb-nR      ;
vc          =   nL:nb-nR-1      ;
P(vr,vc)    =   eye(length(vc)) ;
vr          =   nb+1-(1:nR)     ;
vc          =   nb-1-(1:nR-1)   ;
P(vr,vc)    =   PR              ;

% now compute B-splines desired values/ derivative.
bsplinemat  =   bsplineM(x, breaks, norder, nderiv, sparsewrd)   ;

bsplinemat  =   bsplinemat*P    ;

if sparsewrd, bsplinemat = sparse(bsplinemat); end
