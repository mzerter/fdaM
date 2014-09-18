function [felsplobj,laplacefd] = ....
                  smooth_FEM_fd_Covar(data, fdobj, lambda, desmat)
% SMOOTH_FEM_FD_COVAR Compute a new solution for a FELspline
%  problem including covariates
%
%     Arguments:
% FELSPLOBJ a FELspline object, constructed by SOLVE_FELSPLINE.
% LAMBDA    a scalar smoothing parameter, MUST be positive
% DATA      (optional) a n-by-2 new set of observations.  DATA(:,1)
%           indexes the points (FELSPLOBJ.POINTS) at which the 
%           values in DATA(:,2) were observed.
% DESMAT    a design matrix
%
%     Output:
% FELSPLOBJ  ...  A FD object of the FEM type defined by the coefficient
%                 vector resulting from smoothing
% LAPLACEFD  ...  A FD object of the FEM type for the value of the 
%                 Laplace operator
%

% Last modified on 9 February 2011 by Laura Sangalli and Jim Ramsay

%  assign defaults to missing arguments

if nargin < 4,  desmat = [];     end
if nargin < 3,  lambda = 1e-12;  end

%  check arguments

if ~isa_fd(fdobj)
   error('FDOBJ is not a FD object');
end

if  ~isa(lambda,'double')
   error('LAMBDA is not numeric')
elseif size(lambda) ~= [1 2]
   error('LAMBDA is not a scalar')
end

%  check data argument

if nargin<3
   data = getdata(fdobj);
elseif size(data,2)~=2
   if size(data,1)~=2
      error('DATA is not a n-by-2 array')
   else
      data = data';
   end
end

%  Construct penalty matrix and 'b' vector for Ax=b.

basisobj  = getbasis(fdobj);
params    = getbasispar(basisobj);
numnodes  = size(params.nodes,1);

nodeStruct.order     = params.order;
nodeStruct.nodes     = params.nodes;
nodeStruct.nodeindex = params.nodeindex;
nodeStruct.J         = params.J;
nodeStruct.metric    = params.metric;

%  ---------------------------------------------------------------
% construct mass matrix K0 
%  ---------------------------------------------------------------

K0 = mass(nodeStruct);

%  ---------------------------------------------------------------
% construct stiffness matrix K1
%  ---------------------------------------------------------------

K1 = stiff1(nodeStruct);

%  ---------------------------------------------------------------
% construct projection matrix on the space spanned by the columns of the 
% design matrix desmat
% ATTENZIONE: is it possible to get the projection matrix without 
% having to compute it as below here?
%  ---------------------------------------------------------------

if ~isempty(desmat)
    Q = qr(desmat, 0);
    H = Q * Q';
else
    H = [];
end

%  ---------------------------------------------------------------
% construct the block diagonal matrix L, having upper left block given by I-H 
% and zero otherwise
%  ---------------------------------------------------------------

inddata = data(:,1);
penalty = zeros(numnodes,1);
penalty(inddata) = 1;
indnodes = 1:numnodes;

%  ---------------------------------------------------------------
% construct vector b for system Ax=b
%  ---------------------------------------------------------------
    
b = sparse(zeros(numnodes*2,1));
b(inddata) = data(:,2);

L = sparse(indnodes,indnodes,penalty);
if ~isempty(desmat)
    L(inddata,inddata) = L(inddata,inddata) - H;
    b(1:numnodes,:)  = L * b(1:numnodes,:);
end

%  ---------------------------------------------------------------
% construct matrix A for system Ax=b.
%  ---------------------------------------------------------------
    
A  = [ L    -lambda*K1; ...
      K1            K0];

% solve system

bigsol = A\b;  

solution = bigsol(indnodes);
s = bigsol(indnodes+numnodes);
    
% Make FELspline object
    
felsplobj = fd(solution, basisobj);
laplacefd = fd(s, basisobj);