function [Bvec,Deviance,stats] = ...
              glm_fda(Xmat, Ymat, distr, lamRmat, Wtvec, Bvec0, addterm)
%GLM_FDA Fits a generalized linear model with regularization. 
%  Arguments:
%
%  YMAT    ... a N by NCURVE matrix of data to be fitted
%  XMAT    ... a N by NBASIS matrix of values of basis functions 
%  DISTR   ... a string indicating which of the five GLM family members
%              is assumed:
%              'normal' or 'gaussian' or 'Gaussian'
%              'binomial' or 'binary' or 'Bernoulli'
%              'poisson'
%              'gamma'
%              'inverse gaussian' or 'inverse Gaussian'
%  LAMRMAT ... a \lambda R, that is, a roughness penalty matrix R of 
%              order equal to the number of basis functions used or number
%              of columns of Xmat multiplied by a scalar roughness 
%              penalty parameter \lambda
%  WTVEC   ... a vector of prior weights, such as the inverses of the
%              relative variance of each observation.
%  BVEC0   ... starting values for regresscion coefficients
%  ADDTERM ... a addterm with a coefficient fixed at 1.0.
%
%  Returns:
%  DFE      ... degrees of freedom for error
%  S        ... theoretical or estimated dispersion parameter
%  SFIT     ... estimated dispersion parameter
%  SE       ... standard errors of coefficient estimates B
%  COEFCORR ... correlation matrix for B
%  COVB     ... estimated covariance matrix for B
%  T        ... t statistics for B
%  P        ... p-values for B
%  RESID    ... residuals
%  RESIDP   ... Pearson residuals
%  RESIDD   ... deviance residuals
%  RESIDA   ... Anscombe residuals
%
%   Last modified 3 February 2013 by Jim Ramsay

%--------------------------------------------------------------------------
%                    Check arguments
%--------------------------------------------------------------------------

if nargin < 3
    error('Number of arguments is less than 3.');
end

%  set default argument values

if nargin < 4, Wtvec    = []; end
if nargin < 5, addterm  = []; end
if nargin < 6, Bvec0    = []; end
if nargin < 7, addterm  = []; end
   
%  dimensions of Xmat

[N, nbasis] = size(Xmat);
if iscell(Ymat)
    [Ntmp, ncurve] = size(Ymat{1});
else
    [Ntmp, ncurve] = size(Ymat);
end
if N ~= Ntmp
    error('XMAT and YMAT do not have the same number of rows.');
end

%  define default weight vector WTVEC and check for positivity

if isempty(Wtvec)
    Wtvec = ones(N,1);
end

if any(Wtvec <= 0)
    error('Non-positive values of WTVEC found.');
end

%--------------------------------------------------------------------------
%                         Process YMAT
%--------------------------------------------------------------------------

M = [];
switch distr
case 'normal'
    %  Note:  Ymat can be any real number, no restrictions
case 'binomial'
    if isnumeric(Ymat) 
        %  If YMAT a matrix, M is taken to be 1 (set below)
        %  and it must be a binary matrix containing only 0's and 1's
        if any(any(Ymat < 0 | Ymat > 1))
            error(['For binomial case, YMAT a single column but ', ...
                   ' contains values other than 0 or 1.']);
        end
    elseif iscell(Ymat) && length(Ymat) == 2
        %  if YMAT is a cell array of length 2, then first cell contains
        %  a matrix containing the number of successes and the second 
        %  cell either contains a matrix of the same size as the matrix
        %  in cell {1} or a single positive integer.
        %  these values or this value is the number of trials M 
        %  for a binomial or bernoulli distribution.
        %  M must be a positive integer.
        Freq = Ymat{1};
        M    = Ymat{2};
        if length(M) == 1
            M = M*ones(N,ncurve);
        end
        if ~all(size(M) == size(Freq))
            error(['DISTR is binomial and matrix M is not the same ', ...
                   'size as matrix FREQ']);
        end
        if any(any(M < 0))
            error(['DISTR is binomial and one or more values in M ', ...
                   'have nonpositive values']);
        end
        if any(any(floor(M) ~= M))
            error(['DISTR is binomial and one or more values in M ', ...
                   'have noninteger values.']);
        end
        %  Redefine YMAT is the proportion of sucesses
        Ymat = Freq./M;
    else
        error(['DISTR is binomial and YMAT has incorrect dimensions ', ...
               ' or is of wrong type.']);
    end
case 'poisson'
    %  Note: Ymat must not contain negative numbers
    if any(any(Ymat < 0))
        error('DISTR is poisson and YMAT contains negative values');
    end
case 'gamma'
    %  Note:  Ymat must contain only positive numbers
    if any(Ymat <= 0)
        error('DISTR is gamma and Y contains nonpositive values');
    end
case 'inverse gaussian'
    %  Note:  Ymat must contain only positive numbers    
    if any(any(Ymat <= 0))
        error(['DISTR is inverse gaussian and Y contains ', ...
               'nonpositive values']);
    end
otherwise
    error('Distribution name is invalid.');
end
%  Set M to vector of 1's in Bernoulli case
if isempty(M),  M = ones(N,ncurve);  end

%--------------------------------------------------------------------------
%  Define anonymous functions according to the distribution of Ymat
%     devFn   ... the deviance or loss function, 
%                 called after convergence is achieved
%     stdFn   ... the scale factor multiplying D eta
%                 called second inside loop        
%     linkFn  ... link function, eta = linkFn(mu),
%                 called prior to loop, maps data space into real line
%     DlinkFn ... derivative of the link function wrt to mu
%                 called first inside loop
%     IlinkFn ... the inverse of the link function IlinkFn[eta] = mu,
%                 called last inside loop, maps eta into data space
%--------------------------------------------------------------------------

switch distr
case 'normal'
    devFn   = @(mu,Ymat) (Ymat - mu).^2;
    stdFn   = @(mu)  ones(size(mu));
    linkFn  = @(mu)  mu;
    DlinkFn = @(mu)  ones(size(mu));
    IlinkFn = @(eta) eta;
case 'binomial'
    devFn   = @(mu,Ymat) 2*M.*(Ymat    .*log((Ymat+(Ymat==0))./mu) + ...
                               (1-Ymat).*log((1-Ymat+(Ymat==1))./(1-mu)));
    stdFn   = @(mu)  sqrt(mu.*(1-mu)./M);
    linkFn  = @(mu)  log(mu./(1-mu));
    DlinkFn = @(mu)  1./(mu.*(1-mu));
    loBnd   = log(eps); 
    upBnd   = -loBnd;
    IlinkFn = @(eta) 1./(1 + exp(-constrain(eta,loBnd,upBnd)));
case 'poisson'
    devFn   = @(mu,Ymat) 2*(Ymat      .*(log((Ymat+(Ymat==0))./mu)) - ...
                           (Ymat - mu));
    stdFn   = @(mu)  sqrt(mu);
    linkFn  = @(mu)  log(mu);
    DlinkFn = @(mu)  1./mu;
    loBnd   = log(eps); 
    upBnd   = -loBnd;
    IlinkFn = @(eta) exp(constrain(eta,loBnd,upBnd));
case 'gamma'
    devFn   = @(mu,Ymat) 2*(-log(Ymat./mu) + (Ymat - mu)./mu);
    stdFn   = @(mu) mu;
    linkFn  = @(mu)  1./mu;
    DlinkFn = @(mu) -1./mu.^2;
    loBnd   = eps; 
    upBnd   = 1/loBnd;
    IlinkFn = @(eta) 1./constrain(eta,loBnd,upBnd);
case 'inverse gaussian'
    devFn   = @(mu,Ymat) ((Ymat - mu)./mu).^2./ Ymat;
    stdFn   = @(mu)  mu.^(3/2);
    loBnd   = eps.^(1/2); 
    upBnd   = 1/loBnd;
    linkFn  = @(mu)  constrain(mu,loBnd,upBnd).^(-2);
    DlinkFn = @(mu)  -2*mu.^(-3);
    IlinkFn = @(eta) constrain(eta,loBnd,upBnd).^(-1/2);
otherwise
    error('Distribution name is invalid.');
end

%%%%%%%%%%%%%%%%%% stopped here on 3 Feb 2013  %%%%%%%%%%%%%%%%%%%%%%%%%%
%  continue to allow for YMAT to be a N by ncurve matrix
%  first task ... fix function constrain ... make it matrixBnd

%--------------------------------------------------------------------------
%                   Initialize mu and eta from Ymat.
%--------------------------------------------------------------------------

% Find a starting value for the mean mu, avoiding boundary values

switch distr
    case 'normal'
        mu = Ymat;
    case 'binomial'
        mu = (M.*Ymat + 0.5)./(M + 1);
    case 'poisson'
        mu = Ymat + 0.25;
    case {'gamma' 'inverse gaussian'}
        mu = max(Ymat, eps); 
    otherwise
        mu = Ymat;
end

% compute eta = E(y) from mu

eta = linkFn(mu);

% figure(2)
% subplot(2,1,1)
% plot(eta)
% ylabel('eta')
% subplot(2,1,2)
% plot(mu)
% ylabel('mu')
% 
% pause


%--------------------------------------------------------------------------
%                    Set up for iterations
%--------------------------------------------------------------------------

iter     = 0;
iterLim  = 100;
warned   = false;
seps     = sqrt(eps);
convcrit = 1e-6;

%  set up starting value Bvec0 if required

if isempty(Bvec0)
    Bvec0 = zeros(p,1);
end
Bvec = Bvec0;

% Enforce limits on mu to guard against an inverse linkFn that doesn't map 
% into the support of the distribution.

switch distr
case 'binomial'
    % mu is a probability, so order one is the natural scale, and eps is a
    % reasonable lower limit on that scale (plus it's symmetric).
    muLims = [eps 1-eps];
case {'poisson' 'gamma' 'inverse gaussian'}
    % Here we don't know the natural scale for mu, so make the lower limit
    % small.  This choice keeps mu^4 from underflowing.  No upper limit.
    muLims = realmin.^.25;
end

%--------------------------------------------------------------------------
%                    start of GLM iteration loop
%--------------------------------------------------------------------------

sqrtwt = sqrt(Wtvec);

while iter <= iterLim
    iter = iter+1;

    % Compute adjusted dependent variable for least squares fit
    
    Deta = DlinkFn(mu);
    stdm = stdFn(mu);
    Zvec = eta + (Ymat - mu).*Deta;

%     figure(3)
%     subplot(2,1,1)
%     plot(Ymat - mu)
%     ylabel('Ymat - mu')
%     subplot(2,1,2)
%     plot(Deta)
%     ylabel('Deta')
    
    % Compute IRLS weights the inverse of the variance function
    
    sqrtw = sqrtwt./(abs(Deta).*stdm);
 
%     figure(1)
%     subplot(2,1,1)
%     plot(Zvec)
%     ylabel('Zvec')
%     subplot(2,1,2)
%     plot(sqrtw)
%     ylabel('sqrtw')

%     % Check sqrtw
%     
%     wtol = max(sqrtw)*eps^(2/3);
%     wtst = (sqrtw < wtol);
%     if any(wtst)
%         wtst = wtst & (sqrtw ~= 0);
%         if any(wtst)
%             sqrtw(wtst) = wtol;
%             if ~warned
%                 warning(...
%                   ['Weights are ill-conditioned.   Data may be badly ' ...
%                    'scaled, or\nthe linkFn function may be inappropriate.']);
%             end
%             warned = true;
%         end
%     end

    % Compute coefficient estimates for this iteration - the IRLS step
    
    Bvec_old   = Bvec;
    if ~isempty(addterm)
        Ymattmp = Zvec - addterm;
    else
        Ymattmp = Zvec;
    end
    Ymatw   = Ymattmp.*sqrtw;
    Xmatw   = full(Xmat).*(sqrtw*ones(1,p));
    Mmat    = Xmatw'*Xmatw + lamRmat;
    Bvec    = Mmat\Xmatw'*Ymatw;
    if ~isempty(addterm)
        eta = Xmat*Bvec + addterm;
    else
        eta = Xmat*Bvec;
    end
    mu  = IlinkFn(eta);

    % Force mean in bounds, in case the linkFn function is faulty
    
    switch distr
    case 'binomial'
        if any(mu < muLims(1) | muLims(2) < mu)
            mu = max(min(mu,muLims(2)),muLims(1));
        end
    case {'poisson' 'gamma' 'inverse gaussian'}
        if any(mu < muLims(1))
            mu = max(mu,muLims(1));
        end
    end

%     figure(2)
%     subplot(2,1,1)
%     plot(eta)
%     ylabel('eta')
%     subplot(2,1,2)
%     plot(mu)
%     ylabel('mu')
    
    % Check stopping conditions
    
    disp(['Iteration ',num2str(iter)])
    disp(['Conv. = ',num2str(max(abs(Bvec-Bvec_old)))])
    
    if (~any(abs(Bvec-Bvec_old) > convcrit * max(seps, abs(Bvec_old)))) 
        break; 
    end
    
%     pause
    
end

%--------------------------------------------------------------------------
%                    end of GLM iteration loop
%--------------------------------------------------------------------------

if iter > iterLim
    warning('Iteration limit reached.');
end

if nargout > 1
    % Sum components of deviance to get the total deviance.
    di       = devFn(mu,Ymat);
    Deviance = sum(Wtvec.*di);
end

%--------------------------------------------------------------------------
%                Compute additional statistics 
%--------------------------------------------------------------------------

if nargout > 2
    % Compute the sum of squares used to estimate dispersion, and the
    % Anscombe residuals.
    switch(distr)
    case 'normal'
        ssr = sum(Wtvec.*(Ymat - mu).^2);
        anscresid = Ymat - mu;
    case 'binomial'
        ssr = sum(Wtvec.*(Ymat - mu).^2./...
            (mu.*(1 - mu)./M));
        t = 2/3;
        anscresid = beta(t,t) * ...
            (betainc(Ymat,t,t)-betainc(mu,t,t))./...
            ((mu.*(1-mu)).^(1/6)./sqrt(M));
    case 'poisson'
        ssr = sum(Wtvec.*(Ymat - mu).^2./...
            mu);
        anscresid = 1.5 * ((Ymat.^(2/3) - mu.^(2/3))./...
            mu.^(1/6));
    case 'gamma'
        ssr = sum(Wtvec.*((Ymat - mu)./mu).^2);
        anscresid = 3 * (Ymat.^(1/3) - mu.^(1/3))./mu.^(1/3);
    case 'inverse gaussian'
        ssr = sum(Wtvec.*((Ymat - mu)./...
            mu.^(3/2)).^2);
        anscresid = (log(Ymat) - log(mu))./mu;
    end

    % Compute residuals, using original count scale for binomial
    
    if (isequal(distr, 'binomial'))
        resid = (Ymat - mu).*N;
    else
        resid  = Ymat - mu;
    end

    dfe = max(N - p, 0);
    stats.beta = Bvec;
    stats.dfe  = dfe;
    if dfe > 0
        stats.sfit = sqrt(ssr / dfe);
    else
        stats.sfit = NaN;
    end
    stats.s = stats.sfit;

    % Find coefficient standard errors and correlations
    if ~isnan(stats.s) 
        [Q,R] = qr(Xmat,0);
        RI = R\eye(p);
        C = RI * RI';  %  correct ... doesn't allow for lamRmat > 0
        C = C * stats.s^2; 
        se = sqrt(diag(C)); 
        se = se(:);   
        stats.covb = zeros(p,p);
        stats.covb = C;
        C = C./(se * se');
        stats.se = zeros(p,1); 
        stats.se = se;
        stats.coeffcorr = C;
        stats.t = NaN(p,1); 
        stats.t = Bvec./se;
        stats.p = 2 * tcdf(-abs(stats.t), dfe);
    else
        stats.se = NaN(size(Bvec),class(Bvec));
        stats.coeffcorr = NaN(length(Bvec),class(Bvec));
        stats.t = NaN(size(Bvec),class(Bvec));
        stats.p = NaN(size(Bvec),class(Bvec));
    end
    stats.resid  = resid;
    stats.residp = (Ymat - mu)./(stdFn(mu) + (Ymat==mu));
    stats.residd = sign(Ymat - mu).*sqrt(max(0,di));
    stats.resida = anscresid;
end

%  ------------------------------------------------------------------------

function x = constrain(x, loBnd, upBnd)
if x < loBnd, x = loBnd; end
if x > upBnd, x = upBnd; end


