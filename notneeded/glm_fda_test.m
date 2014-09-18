
%  test function glm_fda.m for function smooth_GLM.m

nbasis   = 13;
basisobj = create_bspline_basis([0,1],nbasis);
Rmat     = eval_penalty(basisobj,2);
lambda   = 1e-4;
lamRmat  = lambda.*Rmat;
fdParobj = fdPar(basisobj, 2, lambda);

N       = 101;
Tvec    = linspace(0,1,N)';
Wtvec   = ones(N,1);
addterm = [];
Xmat    = eval_basis(Tvec,basisobj);

Bvec0 = [];

%  specify distribution

distr = 'normal';
distr = 'binomial';
distr = 'poisson';
distr = 'gamma';
distr = 'inverse gaussian';

%  set up data for this distribution

switch distr
case 'normal'
    Ymat = zeros(N,2);
    eta1  = sin(2*pi*Tvec);
    Ymat(:,1) = eta1 + randn(N,1)*0.2;
    eta2  = cos(2*pi*Tvec);
    Ymat(:,2) = eta2 + randn(N,1)*0.2;
case 'binomial'
    M = 5;
    eta1 = 3*sin(2*pi*Tvec)*ones(1,M);
    eta2 = 3*cos(2*pi*Tvec)*ones(1,M);
    Umat1 = rand(N,M);
    Umat2 = rand(N,M);
    Ymattmp = zeros(N,2);
    for i=1:N
        Uveci = Umat1(i,:);
        Ymati = 1./(1+exp(-eta1(i,:)));
        Ycnti = zeros(1,M);
        Ycnti(Ymati >= Uveci) = 1;
        Ymattmp(i,1) = sum(Ycnti);
        Uveci = Umat2(i,:);
        Ymati = 1./(1+exp(-eta2(i,:)));
        Ycnti = zeros(1,M);
        Ycnti(Ymati >= Uveci) = 1;
        Ymattmp(i,2) = sum(Ycnti);
    end
    Ymat  = Ymattmp;
    Wtvec = M*ones(N,1);
case 'poisson'
    Ymat      = zeros(N,2);
    eta1      = sin(2*pi*Tvec);
    Ymat(:,1) = exp(eta1 + randn(N,1)*0.2);
    eta2      = cos(2*pi*Tvec);
    Ymat(:,2) = exp(eta2 + randn(N,1)*0.2);
case 'gamma'
    Ymat      = zeros(N,2);
    eta1      = sin(2*pi*Tvec).^2 + 1;
    Ymat(:,1) = 1./(eta1 + randn(N,1)*0.2);
    eta2      = cos(2*pi*Tvec).^2 + 1;
    Ymat(:,2) = 1./(eta2 + randn(N,1)*0.2);
case 'inverse gaussian'
    Ymat      = zeros(N,2);
    eta1   = sin(2*pi*Tvec).^2 + 1;
    Ymat(:,1) = 1./sqrt(eta1 + randn(N,1)*0.2);
    eta2   = cos(2*pi*Tvec).^2 + 1;
    Ymat(:,2) = 1./sqrt(eta2 + randn(N,1)*0.2);
otherwise
    error('Distribution name is invalid.');
end

%  do the smoothing and display deviances

[Bvec,Deviance,stats] = ...
              glm_fda(Xmat, Ymat, distr, lamRmat, Wtvec, Bvec0, addterm);
[fdobj, Bvec, df, gcv, Deviance, dev, stats] = ...
                       smooth_GLM(Tvec, Ymat, fdParobj, ...
                      'family', distr, 'weight', Wtvec);
disp(['degrees of freedom for error    = ',num2str(stats.dfe)])
disp(['Deviance                        = ',num2str(Deviance)])
disp(['dispersion parameter            = ',num2str(stats.s)])

%  plot the smooths

figure(1)
subplot(1,1,1)
Yhat = Xmat*Bvec;

switch distr
    case 'normal'
        plot(Tvec, Ymat, 'o', Tvec, Yhat, '-')
        xlabel('\fontsize{13} t')
        ylabel('\fontsize{13} \mu(t)')
        title(['\fontsize{16} ',distr])
    case 'binomial'
        Phat = 1./(1+exp(-Yhat));
        plot(Tvec, Ymat./M, 'o', Tvec, Phat, 'b-')
        xlabel('\fontsize{13} t')
        ylabel('\fontsize{13} \mu(t)')
        title(['\fontsize{16} ',distr])
   case 'poisson'
        plot(Tvec, Ymat, 'o', Tvec, exp(Yhat), 'b-')
        xlabel('\fontsize{13} t')
        ylabel('\fontsize{13} \mu(t)')
        title(['\fontsize{16} ',distr])
    case 'gamma'
        plot(Tvec, Ymat, 'o', Tvec, 1./Yhat, 'b-')
        xlabel('\fontsize{13} t')
        ylabel('\fontsize{13} \mu(t)')
        title(['\fontsize{16} ',distr])
    case 'inverse gaussian'
        plot(Tvec, Ymat, 'o', Tvec, 1./sqrt(Yhat), 'b-')
        xlabel('\fontsize{13} t')
        ylabel('\fontsize{13} \mu(t)')
        title(['\fontsize{16} ',distr])
end

%       'dfe'       degrees of freedom for error
%       's'         theoretical or estimated dispersion parameter
%       'sfit'      estimated dispersion parameter
%       'se'        standard errors of coefficient estimates B
%       'coeffcorr' correlation matrix for B
%       'covb'      estimated covariance matrix for B
%       't'         t statistics for B
%       'p'         p-values for B
%       'resid'     residuals
%       'residp'    Pearson residuals
%       'residd'    deviance residuals
%       'resida'    Anscombe residuals


