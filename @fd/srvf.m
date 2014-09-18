function srvfd = srvf(fdParobj)
%  The signed square root of the first derivative of the functions in
%  FDPAROBJ.

%  Caution!  The square root velocity transformation has infinite slope
%     at points where the derivative of a function crosses the zero axis.
%  This can only be approximated by a B-spline basis system, and only if
%  the basis is sufficiently rich or powerful to compute a large slope
%  value.  Be sure you use lots of knots!  This function will use ONLY the
%  basis that you use in the fdPar argument.

%  Last modified 14 November 2012

if strcmp(class(fdParobj), 'fd')
    fdParobj = fdPar(fdParobj, 2, 1e-10);
end

fdParobj = fdParcheck(fdParobj);
fdobj    = getfd(fdParobj);
basisobj = getbasis(fdobj);

%  test the basis for being of B-spline type

if ~strcmp('bspline',getbasistype(basisobj))
    error('FDOBJ does not have a spline basis.');
end

nbasis   = getnbasis(basisobj);
rangeval = getbasisrange(basisobj);

%  Number of points at which to evaluate the signed square root.  

nmesh = max([10*nbasis+1,501]);

%  determine number curves and variables

coefmat = getcoef(fdobj);
coefd   = size(coefmat);
ncurve  = coefd(2);
if length(coefd) == 2
    nvar = 1;
else
    nvar = coefd(3);
end

if nvar > 1
    error('multivariate svrf not implemented.');
end

%  evaluate derivative of the function over this mesh

tmesh = linspace(rangeval(1),rangeval(2),nmesh);
Dfmat = eval_fd(tmesh, fdobj, 1);

srvfmat = zeros(nmesh,ncurve);
for icurve=1:ncurve
    Dfveci = Dfmat(:,icurve);
    srvfveci = zeros(nmesh,1);
    posind = find(Dfveci > 0);
    srvfveci(posind) = sqrt(Dfveci(posind));
    negind = find(Dfveci < 0);
    srvfveci(negind) = -sqrt(-Dfveci(negind));
    srvfmat(:,icurve) = srvfveci;               
end
    
srvfd = smooth_basis(tmesh, srvfmat, fdParobj);
