glm.fda <- function(Xmat, Ymat, distr, lamRmat, Wtvec=NULL, Bvec0=NULL, addterm=NULL) 
{
  
  #GLM.FDA Fits a generalized linear model with regularization.
  #  This function is called by function smooth.GLM
  #  Arguments
  #
  #  XMAT    ... a N by NBASIS matrix of values  
  # of basis functions
  #  YMAT    ... May be
  #                a N by NCURVE matrix of data to be fitted
  #                or, in the binomial case with local sample sizes M.i,
  #                a cell array of length 2, the first of which cantains
  #                the matrix above containing observed frequencies,
  #                and the second of which contains the corresponding
  #                sample sizes.  Note that in the binary or Bernoulli case,
  #                Y is a matrix of 1's and 0's and the M's are
  #                taken to be 1's.
  #  DISTR   ... a string indicating which of the five GLM family members
  #              is assumed
  #              'normal' or 'gaussian' or 'Gaussian'
  #              'binomial' or 'binary' or 'Bernoulli'
  #              'poisson'
  #              'gamma'
  #              'inverse gaussian' or 'inverse Gaussian'
  #              or a cell array of length(N) with each cell containing
  #              a specification of the GLM family of a single observation.
  #  LAMRMAT ... a \lambda R, that is, a roughness penalty matrix R of
  #              order equal to the number of basis functions used or number
  #              of columns of Xmat multiplied by a scalar roughness
  #              penalty parameter \lambda
  #  WTVEC   ... a vector of prior weights, such as the inverses of the
  #              relative variance of each observation.
  #  BVEC0   ... starting values for regresscion coefficients
  #  ADDTERM ... a addterm with a coefficient fixed at 1.0.
  #
  #  Returns
  #  BVEC     ... Final estimate of coefficients
  #  DEVIANCE ... Deviance values
  #
  #   Last modified 15 September 2014 by Jim Ramsay
  
  #--------------------------------------------------------------------------
  #                    Check arguments
  #--------------------------------------------------------------------------
  
  # if nargin < 3
  #     error('Number of arguments is less than 3.')
  # end
  # 
  # #  set default argument values
  # 
  # if nargin < 5, Wtvec    = [] end
  # if nargin < 6, Bvec0    = [] end
  # if nargin < 7, addterm  = [] end
  # 
  # if nargin < 4
  #     error('Less than four arguments for function glm.fda.')
  # end
  
  #  dimensions of Xmat
  eps      = 0.1
  
  Size_Xmat = dim(Xmat)
  N = Size_Xmat[1]
  nbasis = Size_Xmat[2]
  
  Size_Ymat = dim(Ymat)
  Ntmp = Size_Ymat[1]
  ncurve = Size_Ymat[2]
  
  # if iscell(Ymat)
  # 
  #     [Ntmp, ncurve] = size(Ymat{1})
  # else
  #     [Ntmp, ncurve] = size(Ymat)
  # end
  
  if (N != Ntmp)
  {
    stop('XMAT and YMAT do not have the same number of rows.')
  }
  
  
  #  define default weight vector WTVEC and check for positivity
  
  if (is.null(Wtvec))
  {
    Wtvec = matrix(1,N,1)
  }
  
  if (any(Wtvec <= 0))
  {
    stop('Non-positive values of WTVEC found.')
  }
  
  #--------------------------------------------------------------------------
  #  Process YMAT and define anonymous functions according to the
  #  distribution of Ymat
  #     devFn   ... the deviance or loss function,
  #                 called after convergence is achieved
  #     stdFn   ... the scale factor multiplying D eta
  #                 called second inside loop
  #     linkFn  ... link function, eta = linkFn(mu),
  #                 called prior to loop, maps data space into real line
  #     DlinkFn ... derivative of the link function wrt to mu
  #                 called first inside loop
  #     IlinkFn ... the inverse of the link function IlinkFn[eta] = mu,
  #                 called last inside loop, maps eta into data space
  # Then set a starting value for the mean mu, avoiding boundary values.
  #--------------------------------------------------------------------------
  
  M = NULL
  if (is.character(distr))
  {
    #  --------------------------------------------------------------------
    #    All observations are in the same family, distr is a string
    #  --------------------------------------------------------------------
    switch(distr,
           'normal' =
{           
  #  Note  Ymat can be any real number, no restrictions
  devFn   = function(mu,Ymat) (Ymat - mu)^2
  stdFn   = function(mu)  matrix(1,dim(mu)[1],dim(mu)[2])
  linkFn  = function(mu)  mu
  DlinkFn = function(mu)  matrix(1,dim(mu)[1],dim(mu)[2])
  IlinkFn = function(eta) eta
  mu      = Ymat
},

'binomial' =
{
  if (is.numeric(Ymat))
  {
    #  If YMAT a matrix, M is taken to be 1 (set below)
    #  and it must be a binary matrix containing only 0's and 1's
    if (any(any(Ymat < 0 | Ymat > 1)))
    {
      stop('For binomial case, YMAT a single column but contains values other than 0 or 1.')
    }
    M = matrix(1,N,ncurve)
  }            
  #             else if (iscell(Ymat) && length(Ymat) == 2)
  #                 #  If YMAT is a cell array of length 2, then first cell
  #                 #  contains a matrix containing the number of successes and
  #                 #  the second cell either contains a matrix of the same
  #                 #  size as the matrix in Ymat{1} or a single positive
  #                 #  integer.
  #                 #  These values or this value is the number of trials M
  #                 #  for a binomial or bernoulli distribution.
  #                 #  M must be a positive integer.
  #                 Freq = Ymat{1}
  #                 M    = Ymat{2}
  #                 if length(M) == 1
  #                     M = M*ones(N,ncurve)
  #                 end
  #                 if ~all(size(M) == size(Freq))
  #                     error(['DISTR is binomial and matrix M is not the same ', ...
  #                         'size as matrix FREQ'])
  #                 end
  #                 if any(any(M < 0))
  #                     error(['DISTR is binomial and one or more values in M ', ...
  #                         'have nonpositive values'])
  #                 end
  #                 if any(any(floor(M) ~= M))
  #                     error(['DISTR is binomial and one or more values in M ', ...
  #                         'have noninteger values.'])
  #                 end
  #                 #  Redefine YMAT is the proportion of sucesses
  #                 Ymat = Freq./M
  else
  {                
    stop('DISTR is binomial and YMAT has incorrect dimensions or is of wrong type.')
  }
  
  devFn   = function(mu,Ymat) 2*M*(Ymat*log((Ymat+(Ymat==0))/mu) + (1-Ymat) * log((1-Ymat+(Ymat==1))/(1-mu)))
  stdFn   = function(mu)  sqrt(mu * (1-mu)/M)
  linkFn  = function(mu)  log(mu/(1-mu))
  DlinkFn = function(mu)  1/(mu * (1-mu))
  loBnd   = log(eps)
  upBnd   = -loBnd
  IlinkFn = function(eta) 1/(1 + exp(-constrain(eta,loBnd,upBnd)))
  mu      = (M*Ymat + 0.5)/(M + 1)
},

'poisson' = 
{
  #  Note Ymat must not contain negative numbers
  if (any(any(Ymat < 0)))
  {
    stop('DISTR is poisson and YMAT contains negative values')
  }
  
  devFn   = function(mu,Ymat) 2*(Ymat*(log((Ymat+(Ymat==0))/mu)) - (Ymat - mu))
  stdFn   = function(mu)  sqrt(mu)
  linkFn  = function(mu)  log(mu)
  DlinkFn = function(mu)  1./mu
  loBnd   = log(eps)
  upBnd   = -loBnd
  IlinkFn = function(eta) exp(constrain(eta,loBnd,upBnd))
  mu      = Ymat + 0.25
},
'gamma' = 
{
  #  Note  Ymat must contain only positive numbers
  if (any(Ymat <= 0))
  {
    stop('DISTR is gamma and Y contains nonpositive values')
  }
  
  devFn   = function(mu,Ymat) 2*(-log(Ymat/mu) + (Ymat - mu)/mu)
  stdFn   = function(mu) mu
  linkFn  = function(mu)  1./mu
  DlinkFn = function(mu) -1./mu.^2
  loBnd   = eps
  upBnd   = 1/loBnd
  IlinkFn = function(eta) 1./constrain(eta,loBnd,upBnd)
  mu      = max(Ymat, eps)
},
'inverse gaussian'=
{
  #  Note  Ymat must contain only positive numbers
  if (any(any(Ymat <= 0)))
  {
    error('DISTR is inverse gaussian and Y contains nonpositive values')
  }
  devFn   = function(mu,Ymat) ((Ymat - mu)/mu)^2/ Ymat
  stdFn   = function(mu)  mu.^(3/2)
  loBnd   = eps^(1/2)
  upBnd   = 1/loBnd
  linkFn  = function(mu) apply(mu,1:2, function(x_input) constrain(x_input,loBnd,upBnd)^(-2))
  DlinkFn = function(mu)  -2*mu^(-3)
  IlinkFn = function(eta) apply(eta,1:2, function(x_input) constrain(x_input,loBnd,upBnd)^(-1/2))
  mu      = Ymat
},        
stop('Distribution name is invalid.')
    )
  }
# else if iscell(distr) && length(distr) == N
#     #  --------------------------------------------------------------------
#     #    Observations can be in different families, distr is a cell array.
#     #  --------------------------------------------------------------------
#     mu      = zeros(N,1)
#     loBnd   = zeros(N,1)
#     upBnd   = zeros(N,1)
#     devFn   = cell(N,1)
#     stdFn   = cell(N,1)
#     linkFn  = cell(N,1)
#     DlinkFn = cell(N,1)
#     IlinkFn = cell(N,1)
#     #  Dealing with the presence of some binomial observations Ymat has
#     #  to be a cell with N rows and 2 columns for all data.  Ugh!
#     binomwrd = iscell(Ymat) && all(size(Ymat) == [N,2])
#     for i=1N
#         distri = distr{i}
#         if ~ischar(distri)
#             error('A distribution specification is not a string.')
#         end
#         switch distri
#             case 'normal'
#                 #  Note  Ymat can be any real number, no restrictions
#                 devFn{i}   = @(mu,Ymat) (Ymat - mu).^2
#                 stdFn{i}   = @(mu)  ones(size(mu))
#                 linkFn{i}  = @(mu)  mu
#                 DlinkFn{i} = @(mu)  ones(size(mu))
#                 IlinkFn{i} = @(eta) eta
#                 mu(i,) = Ymat(i,)
#             case 'binomial'
#                 if all(isnumeric(Ymat(i,)))
#                     #  If YMAT a matrix, M is taken to be 1 (set below)
#                     #  and it must be a binary matrix containing only
#                     #0's and 1's
#                     if any(Ymat(i,) < 0 | Ymat(i,) > 1)
#                         error(['For binomial case, YMAT a single column but ', ...
#                             ' contains values other than 0 or 1.'])
#                     end
#                 elseif binomwrd
#                     Freqi = Ymat{i,1}
#                     Mi    = Ymat{i,2}
#                     if length(Mi) == 1
#                         Mi = Mi*ones(1,ncurve)
#                     end
#                     if ~all(size(Mi) == size(Freqi))
#                         error(['DISTR is binomial and matrix M is not the same ', ...
#                             'size as matrix FREQ'])
#                     end
#                     if any(any(Mi < 0))
#                         error(['DISTR is binomial and one or more values in M ', ...
#                             'have nonpositive values'])
#                     end
#                     if any(any(floor(Mi) ~= Mi))
#                         error(['DISTR is binomial and one or more values in M ', ...
#                             'have noninteger values.'])
#                     end
#                     #  Redefine YMAT is the proportion of sucesses
#                     Ymat(i,) = (Freqi./Mi)
#                 else
#                     error(['DISTR is binomial and YMAT has incorrect dimensions ', ...
#                         ' or is of wrong type.'])
#                 end
#                 devFn{i}   = @(mu,Ymat) 2*M.*(Ymat.*log((Ymat+(Ymat==0))./mu) + ...
#                     (1-Ymat).*log((1-Ymat+(Ymat==1))./(1-mu)))
#                 stdFn{i}   = @(mu)  sqrt(mu.*(1-mu)./M)
#                 linkFn{i}  = @(mu)  log(mu./(1-mu))
#                 DlinkFn{i} = @(mu)  1./(mu.*(1-mu))
#                 loBnd(i)   = log(eps)
#                 upBnd(i)   = -loBnd(i)
#                 IlinkFn{i} = @(eta) 1./(1 + exp(-constrain(eta,loBnd,upBnd)))
#                 mu(i)      = (M(i).*Ymat(i) + 0.5)./(M(i) + 1)
#             case 'poisson'
#                 #  Note Ymat must not contain negative numbers
#                 if any(Ymat(i,) < 0)
#                     error('DISTR is poisson and YMAT contains negative values')
#                 end
#                 devFn{i}   = @(mu,Ymat) 2*(Ymat.*(log((Ymat+(Ymat==0))./mu)) - ...
#                     (Ymat - mu))
#                 stdFn{i}   = @(mu)  sqrt(mu)
#                 linkFn{i}  = @(mu)  log(mu)
#                 DlinkFn{i} = @(mu)  1./mu
#                 loBnd(i)   = log(eps)
#                 upBnd(i)   = -loBnd(i)
#                 IlinkFn{i} = @(eta) exp(constrain(eta,loBnd,upBnd))
#                 mu(i,)    = Ymat(i,) + 0.25
#             case 'gamma'
#                 #  Note  Ymat must contain only positive numbers
#                 if any(Ymat(i) <= 0)
#                     error('DISTR is gamma and Y contains nonpositive values')
#                 end
#                 devFn{i}   = @(mu,Ymat) 2*(-log(Ymat./mu) + (Ymat - mu)./mu)
#                 stdFn{i}   = @(mu) mu
#                 linkFn{i}  = @(mu)  1./mu
#                 DlinkFn{i} = @(mu) -1./mu.^2
#                 loBnd(i)   = eps
#                 upBnd(i)   = 1/loBnd(i)
#                 IlinkFn{i} = @(eta) 1./constrain(eta,loBnd,upBnd)
#                 mu(i,)    = max(Ymat(i,), eps)
#             case 'inverse gaussian'
#                 #  Note  Ymat must contain only positive numbers
#                 if any(Ymat(i,) <= 0)
#                     error(['DISTR is inverse gaussian and Y contains ', ...
#                         'nonpositive values'])
#                 end
#                 devFn{i}   = @(mu,Ymat) ((Ymat - mu)./mu).^2./ Ymat
#                 stdFn{i}   = @(mu)  mu.^(3/2)
#                 loBnd(i)   = eps.^(1/2)
#                 upBnd(i)   = 1/loBnd(i)
#                 linkFn{i}  = @(mu)  constrain(mu,loBnd,upBnd).^(-2)
#                 DlinkFn{i} = @(mu)  -2*mu.^(-3)
#                 IlinkFn{i} = @(eta) constrain(eta,loBnd,upBnd).^(-1/2)
#                 mu(i,)    = Ymat(i,)
#             otherwise
#                 error('Distribution name is invalid.')
#         end
#     end
else
{
  stop('DISTR is neither a string or a cell array of length N.')
}


#--------------------------------------------------------------------------
#                   Initialize mu and eta from Ymat.
#--------------------------------------------------------------------------

# compute eta = E(y) from mu

if (is.character(distr))
{
  eta = linkFn(mu)
}
# else
#     {    
#       eta  = matrix(0,N,nurve)
#       Deta = matrix(0,N,nurve)
#       stdm = matrix(0,N,nurve)
#       for (i in 1:N)
#           {
#             linkFni  = linkFn{i}
#             eta(i,) = linkFni(mu(i,))
#           }
#       
#     }

#--------------------------------------------------------------------------
#                        Set up for iterations
#--------------------------------------------------------------------------

iter     = 0
iterLim  = 100

seps     = sqrt(eps)
convcrit = 1e-6
sqrtwt   = sqrt(Wtvec)

#  set up starting value Bvec0 if required

if (is.null(Bvec0))
{
  Bvec0 = matrix(0,nbasis,ncurve)
}

Bvec = Bvec0

# Enforce limits on mu to guard against an inverse linkFn that doesn't map
# into the support of the distribution.

switch(distr,
       'binomial' = 
{
  # mu is a probability, so order one is the natural scale, and eps is a
  # reasonable lower limit on that scale (plus it's symmetric).
  muLims = c(eps,1-eps)
},
# Here we don't know the natural scale for mu, so make the lower limit
# small.  This choice keeps mu^4 from underflowing.  No upper limit.
'poisson'           = {muLims = realmin.^.25},
'gamma'             = {muLims = realmin.^.25},
'inverse gaussian'  = {muLims = realmin.^.25}
)
#--------------------------------------------------------------------------
#                       Start of GLM iteration loop
#--------------------------------------------------------------------------

while (iter <= iterLim)
{
  iter = iter+1
  printf("Iteration %d \n",iter)
  
  
  # Compute adjusted dependent variable for least squares fit
  
  if (is.character(distr))
  {
    Deta = DlinkFn(mu)
    stdm = stdFn(mu)
  }
  #     else
  #         {
  #           for (i in 1:N)
  #             {
  #               DlinkFni  = DlinkFn{i}
  #               stdFni    = stdFn{i}
  #               mui       = mu(i,)
  #               Deta(i,) = DlinkFni(mui)
  #               stdm(i,) = stdFni(mui)
  #             }
  #         }
  
  Zvec = eta + (Ymat - mu) * Deta
  
  # Compute IRLS weights the inverse of the variance function
  
  sqrtw = (sqrtwt %*% matrix(1,1,ncurve))/(abs(Deta)*stdm)
  
  # Compute coefficient estimates for this iteration - the IRLS step
  
  Bvec.old   = Bvec
  if (!is.null(addterm))
  {
    Ymattmp = Zvec - addterm
  }
  else
  {
    Ymattmp = Zvec
  }
  
  Ymatw   = Ymattmp*sqrtw
  Xmatw   = Xmat*(sqrtwt %*% matrix(1,1,nbasis))
  
  if (is.null(lamRmat))
  {
    Mmat = t(Xmatw)%*%Xmatw
  }
  else
  {
    Mmat = t(Xmatw) %*% Xmatw + lamRmat
  }
  
  Bvec    = solve(Mmat,(t(Xmatw) %*% Ymatw))
  
  if (!is.null(addterm))
  {
    eta = Xmat %*% Bvec + addterm
  }
  else
  {
    eta = Xmat %*% Bvec
  }
  
  if (is.character(distr))
  {
    mu = IlinkFn(eta)
  }
  #     else
  #         {
  #           for (i in 1:N)
  #             IlinkFni  = IlinkFn{i}
  #             mu(i,) = IlinkFni(eta(i,))
  #         }
  
  
  # Force mean in bounds, in case the linkFn function is faulty
  
  if (is.character(distr))
  {
    # Hidden function for reusing code
    f_ = function() 
    {
      if (any(any(mu < muLims(1))))
      {
        for (j in 1:m)
        {
          mu[,j] = max(mu[,j],muLims[1])
        }
      }             
    }
    
    switch(distr,
           'binomial'=
{ 
  if (any(any(mu < muLims[1] | muLims[2] < mu)))
  {
    for (j in 1:m)
    {
      mu[,j] = max(min(mu[,j],muLims[2]),muLims[1])
    }
  }
},
'poisson' =
{
  f_()
},
'gamma' = 
{
  f_()
},
'inverse gaussian' = 
{
  f_()
}
    )
  }
#     else
#       {
#         for (i in 1:N)
#             distri = distr{i}
#             switch distri
#                 case 'binomial'
#                     if any(any(mu(i,) < muLims(1) | muLims(2) < mu(i,)))
#                         for j=1m
#                             mu(i,j) = max(min(mu(i,j),muLims(2)),muLims(1))
#                         end
#                     end
#                 case {'poisson' 'gamma' 'inverse gaussian'}
#                     if any(any(mu(i,) < muLims(1)))
#                         for j=1m
#                             mu(i,j) = max(mu(i,j),muLims(1))
#                         end
#                     end
#             end
#         end
#       }


# Check stopping conditions

if (max(max(abs(Bvec-Bvec.old))) < convcrit*max(max(abs(Bvec.old)))) break

#     pause

}
#--------------------------------------------------------------------------
#                    end of GLM iteration loop
#--------------------------------------------------------------------------


if (iter > iterLim) warning('smooth.GLM:','Iteration ','Iteration limit reached.')


# if (nargout > 1)
#     # Sum components of deviance to get the total deviance.
#     if ischar(distr)
#         di       = devFn(mu,Ymat)
#         Deviance = sum((Wtvec*ones(1,ncurve)).*di)
#     else
#         Deviance = zeros(N,ncurve)
#         for i=1N
#             devFni = devFn{i}
#             di = devFni(mu(i,),Ymat(,i))
#             Deviance(i,) = sum((Wtvec(i)*ones(1,ncurve)).*di)
#         end
#     end
# end
# if nargout > 2
#     stats = []
# end

if (is.character(distr))
{
  di = devFn(mu,Ymat)
  Deviance = sum((Wtvec %*% matrix(1,1,ncurve)) * di)
}
else
{
  for (i in 1:N)
  {
    devFni  = devFn[i]
    di      = devFni(mu[i,],Ymat[,i])
    Deviance[i,] = sum((Wtvec[i] %*% matrix(1,1,ncurve)) * di)
  }
}

return(list(Bvec = Bvec, Deviance = Deviance))
}