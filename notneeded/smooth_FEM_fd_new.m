function [felsplobj,laplacefd] = smooth_FEM_fd_new(data,fdobj,lambda)
% SMOOTH_FEM_FD Compute a new solution for a FELspline
%  problem 
%
%     Arguments:
% FELSPLOBJ a FELspline object, constructed by SOLVE_FELSPLINE.
% LAMBDA    a scalar smoothing parameter
% DATA      (optional) a n-by-2 new set of observations.  DATA(:,1)
%           indexes the points (FELSPLOBJ.POINTS) at which the 
%           values in DATA(:,2) were observed.
%
%     Output:
% FELSPLOBJ  ...  A FD object of the FEM type defined by the coefficient
%                 vector resulting from smoothing
% LAPLACEFD  ...  A FD object of the FEM type for the value of the 
%                 Laplace operator 
%
% Last modified on 26 August 2010 by Jim Ramsay
% Last modified on 8 February 2011 by Laura Sangalli

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
   data=getdata(fdobj);
elseif size(data,2)~=2
   if size(data,1)~=2
      error('DATA is not a n-by-2 array')
   else
      data=data';
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
% construct the penalty matrix P with ones on diagonal at data points
%  ---------------------------------------------------------------

penalty = zeros(numnodes,1);
penalty(data(:,1)) = 1;
indnodes = 1:numnodes;
P = sparse(indnodes,indnodes,penalty);

%  ---------------------------------------------------------------
% construct the solution
%  ---------------------------------------------------------------

%  ---------------------------------------------------------------
% construct vector b for system Ax=b
%  ---------------------------------------------------------------
    
b = sparse(zeros(numnodes*2,1));
b(data(:,1)) = data(:,2);
    
%  ---------------------------------------------------------------
% construct matrix A for system Ax=b.
%  ---------------------------------------------------------------
    
A = [ P  -lambda*K1; ...
     K1         K0];
    
% solve system
    
bigsol = A\b;
    
solution = bigsol(indnodes);
s = bigsol(indnodes+numnodes);
    
% Make FELspline object
    
felsplobj = fd(solution, basisobj);
laplacefd = fd(s, basisobj);

