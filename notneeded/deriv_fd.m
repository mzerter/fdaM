function derivfdobj = deriv_fd(fdobj, Lfdobj)
%  DERIV_FD applies linear differential operator LFD 
%  to functional data object FDOBJ
%  and returns the result as functional data object 
%  DERIVFDOBJ.

%  Last modified 24 September 2005

if ~strcmp(class(fdobj), 'fd') 
    error('Argument  FD not a functional data object.')
end

if nargin < 2, Lfdobj = int2Lfd(1);  end

Lfdobj   = int2Lfd(Lfdobj);

basisobj = getbasis(fdobj);
nbasis   = getnbasis(basisobj);
rangeval = getbasisrange(basisobj);

evalarg  = linspace(rangeval(1), rangeval(2), 10*nbasis+1)';
Lfdmat   = eval_fd(evalarg, fdobj, Lfdobj);

Lfdcoef  = project_basis(Lfdmat, evalarg, basisobj);

Dfdnames = getnames(fdobj);
Dfdnames{3} = ['D',Dfdnames{3}];

derivfdobj = fd(Lfdcoef, basisobj, Dfdnames);
