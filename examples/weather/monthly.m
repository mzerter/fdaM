addpath ('..\..')

%  Last modified  14 December 2012

%  -----------------------------------------------------------------------
%                      Monthly Weather data
%  -----------------------------------------------------------------------

%  ----------------   input the data and set up labels -------------------

fid = fopen('monthtemp.dat','rt');
tempvec = fscanf(fid,'%f');
tempmat = reshape(tempvec, [12, 35]);

fid = fopen('monthprec.dat','rt');
precvec = fscanf(fid,'%f');
precmat = reshape(precvec, [12, 35]);

meteonames =   [  'St. Johns    '; 'Charlottetown'; 'Halifax      '; ...
                  'Sydney       '; 'Yarmouth     '; 'Fredericton  '; ...
                  'Arvida       '; 'Montreal     '; 'Quebec City  '; ...
                  'Schefferville'; 'Sherbrooke   '; 'Kapuskasing  '; ...
                  'London       '; 'Ottawa       '; 'Thunder Bay  '; ...
                  'Toronto      '; 'Churchill    '; 'The Pas      '; ...
                  'Winnipeg     '; 'Prince Albert'; 'Regina       '; ...
                  'Beaverlodge  '; 'Calgary      '; 'Edmonton     '; ...
                  'Kamloops     '; 'Prince George'; 'Prince Rupert'; ...
                  'Vancouver    '; 'Victoria     '; 'Dawson       '; ...
                  'Whitehorse   '; 'Frobisher Bay'; 'Inuvik       '; ...
                  'Resolute     '; 'Yellowknife  '];

monthletter = ['J'; 'F'; 'M'; 'A'; 'M'; 'J'; 'J'; 'A'; 'S'; 'O'; 'N'; 'D'];

months = ['Jan'; 'Feb'; 'Mar'; 'Apr'; 'May'; 'Jun';
          'Jul'; 'Aug'; 'Sep'; 'Oct'; 'Nov'; 'Dec'];

%  indices for weather stations in each of four climate zones

atlindex = [1,2,3,4,5,6,7,8,9,10,11,13,14,16];
pacindex = [25,26,27,28,29];
conindex = [12,15,17,18,19,20,21,22,23,24,30,31,35];
artindex = [32,33,34];

%  -----------------  set up argument values and weights  ---------------

monthtime = linspace(0.5, 11.5, 12)';  %  mid points of months
weeks     = linspace(0,12,53)';        %  month values roughly in weeks

%  -------------------  set up the fourier basis object -----------------

nbasis = 12;
monthbasis = create_fourier_basis([0,12], nbasis);

%  ------------  set up the harmonic acceleration operator  -------------

%  The operator is Lx(t) = (pi/6)^2 Dx(t) + D^3x(t)
%  This operator has the shifted sinusoid as its null space, and it will
%  be used in situations where we want to smooth towards a function
%  that is a combination of a constant, a sine, and a cosine, with period
%  12 months.

Lbasis  = create_constant_basis([0,12]);  %  create a constant basis
Lcoef   = [0,(pi/6)^2,0];                 %  set up three coefficients
wfd     = fd(Lcoef,Lbasis);   % define an FD object for weight functions
wfdcell = fd2cell(wfd);       % convert the FD object to a cell object
harmaccelLfd = Lfd(3, wfdcell);  %  define the operator object

%  --------------  make temperature fd object, don't smooth  ------------

monthfdPar = fdPar(monthbasis, harmaccelLfd, 1e-8);
tempfd = smooth_basis(monthtime, tempmat, monthfdPar);

tempfd_fdnames{1} = 'Months';
tempfd_fdnames{2} = 'Station';
tempfd_fdnames{3} = 'Deg C';
tempfd = putnames(tempfd, tempfd_fdnames);

%  plot temperature functions

plot(tempfd);
title('Temperature Functions');

%  plot each curve along with data

plotfit_fd(tempmat, monthtime, tempfd, meteonames)

%  --------------  make precipitation fd object, don't smooth  ------------

precfd = smooth_basis(monthtime, precmat, monthfdPar);

precfd_fdnames{1} = 'Months';
precfd_fdnames{2} = 'Station';
precfd_fdnames{3} = 'mm';
precfd = putnames(precfd, precfd_fdnames);

%  plot precipitation functions

plot(precfd);
title('Precipitation Functions')

%  ---------------  plot temperature for 4 stations --------------------

stnindex = [8, 24, 27, 34];

plot(tempfd(stnindex))
axis([0,12,-35,25])
title('\fontsize{16} Selected Temperature Functions')
legend(meteonames(stnindex,:))

%  this plot is for the FDA lecture series

xfine = linspace(0,12,101)';
tempfine = eval_fd(xfine, tempfd(stnindex));

plot(xfine, tempfine, '-', monthtime, tempmat(:,stnindex), 'o')
xlabel('\fontsize{12} Month')
ylabel('\fontsize{12} Deg C')
title('\fontsize{16} Temperatures')
legend(meteonames(stnindex,:))

print -dpsc2 'c:/MyFiles/talks/fdacourse/figs/fourtemps.ps'

%  plot harmonic acceleration values

xfine = linspace(0,12,101)';
tempLfine = eval_fd(xfine, tempfd(stnindex), harmaccelLfd);

plot(xfine, tempLfine, '-', [0,12], [0,0], ':')
xlabel('\fontsize{12} Month')
ylabel('\fontsize{12} L(Deg C)')
title('\fontsize{16} Harmonic Accelerations')
legend(meteonames(stnindex,:))

print -dpsc2 'c:/MyFiles/talks/fdacourse/figs/fourLtemps.ps'

% plot temperature for Montreal (Figure 1.13 in book}

tempD2fine = eval_fd(tempfd(8), xfine,     2);
tempD2mat  = eval_fd(tempfd(8), monthtime, 2);

subplot(1,1,1)
plot(xfine, tempfine(:,1), '-')
xlabel('\fontsize{12} Month')
ylabel('\fontsize{12} deg C')
title('\fontsize{16} Montreal Temperature')
text(monthtime, tempmat(:,8), monthletter)
axis([0,12,-15,25])
axis('square')

print -dpsc2 'c:/MyFiles/talks/fdacourse/figs/Mtlplot1.ps'

subplot(1,1,1)
plot(tempfine(:,1), tempD2fine(:,1), '-', [-15,25],[0,0],':')
xlabel('\fontsize{12} deg C')
ylabel('\fontsize{12} D^2 deg C')
text(tempmat(:,8), tempD2mat, monthletter)
axis([-15,25,-6,9])
axis('square')

print -dpsc2 'c:/MyFiles/talks/fdacourse/figs/Mtlplot2.ps'

%  ----------------------------------------------------------------------
%                Descriptive Statistics Functions
%  ----------------------------------------------------------------------

%  --  compute and plot mean and standard deviation of temperature ------

tempmeanfd   = mean(tempfd);

tempstddevfd = std(tempfd);

plot(tempmeanfd);
title('Mean Temperature');

plot(tempstddevfd);
title('Std. dev. of Temperature');
axis([0 12 0 10])

%  -----  plot the temperature variance bivariate function  --------

tempvarbifd = var(tempfd);

tempvarmat = eval_bifd(tempvarbifd, weeks, weeks);

subplot(1,1,1);
surf(tempvarmat);
xlabel('Weeks')
ylabel('Weeks')
zlabel('Covariance')
title('Temperature Variance-Covariance Function')

%  -------------  plot the correlation function  -----------------------

tempstddev = sqrt(diag(tempvarmat));
tempcormat = tempvarmat./(tempstddev*tempstddev');

subplot(1,1,1);
surf(tempcormat);
xlabel('Weeks')
ylabel('Weeks')
zlabel('Covariance')
title('Temperature Correlation Function')
axis([0,53,0,53,0,1])

%  -----  plot the precipitation variance bivariate function  --------

precvarbifd = var(precfd);

precvarmat = eval_bifd(precvarbifd, weeks, weeks);

subplot(1,1,1);
surf(precvarmat);
xlabel('Weeks')
ylabel('Weeks')
zlabel('Covariance')
title('Precipitation Variance-Covariance Function')

%  -------------  plot the correlation function  -----------------------

precstddev = sqrt(diag(precvarmat));
preccormat = precvarmat./(precstddev*precstddev');

subplot(1,1,1);
surf(preccormat);
xlabel('Weeks')
ylabel('Weeks')
zlabel('Covariance')
title('Precipitation Correlation Function')
axis([0,53,0,53,0,1])

%  -----  compute and plot the covariance between temp. and prec.  --

covbifd = var(tempfd, precfd);
covmat  = eval_bifd(covbifd,weeks,weeks);

contour(covmat)
xlabel('Months')
ylabel('Months')
title('Covariance')

surf(covmat)
xlabel('Months')
ylabel('Months')
zlabel('Covariance')

%  -----  compute and plot the correlation between temp. and prec.  --

cormat  = covmat./(tempstddev*precstddev');

contour(cormat)
xlabel('Months')
ylabel('Months')
title('Correlation')

surf(cormat)
xlabel('Months')
ylabel('Months')
zlabel('Correlation')

%  ----------------------------------------------------------------------
%               Principal Components Analysis of temperature  
%  ----------------------------------------------------------------------

%  ------------------  center and plot the temperature data  ------------

tempcenterfd = center(tempfd);

plot(tempcenterfd);
title('Centered Temperature');

%  Penalize harmonic acceleration

nharm  = 4;
lambda = 1e-3;

temppcastr = pca_fd(tempfd, nharm,monthfdPar);

%  plot harmonics

subplot(1,1,1)
plot_pca_fd(temppcastr);

tempharm = temppcastr.harmfd;
tempvarprop = temppcastr.varprop;
monthfine = linspace(0,12,101)';
tempharmfine = eval_fd(tempharm,monthfine);
tempharmmat  = eval_fd(tempharm,monthtime);
tempharmfine(:,3) = -tempharmfine(:,3);
tempharmmat(:,3)  = -tempharmmat(:,3);

for j=1:4
    subplot(2,2,j)
    plot(monthfine,tempharmfine(:,j),'g-', ...
         monthtime,tempharmmat(:,j), 'bo', ...
         [0,12],[0,0],'r:')
    xlabel('Month')
    ylabel('Component')
    title(['PC ',num2str(j),' (',num2str(round(tempvarprop(j)*1000)/10),'%)'])
    axis([0,12,-0.6,0.6])
end
 
print -dpsc2 'c:/Myfiles/talks/fdacourse/figs/temppca1.ps'

con = 10.*ones(1,4);
tempmeanfine = eval_fd(tempmeanfd,monthfine);
tempmeanvec  = eval_fd(tempmeanfd,monthtime);

for j=1:4
    subplot(2,2,j)
    yplus = tempmeanvec + con(j).*tempharmmat(:,j);
    ymins = tempmeanvec - con(j).*tempharmmat(:,j);
    plot(monthfine,tempmeanfine,'b-',...
         monthtime,yplus,'g--',...
         monthtime,ymins,'r--')
    xlabel('Month')
    ylabel('Component')
    title(['PC ',num2str(j),' (',num2str(round(tempvarprop(j)*1000)/10),'%)'])
    axis([0,12,-20,20])
end

print -dpsc2 'c:/Myfiles/talks/fdacourse/figs/temppca2.ps'

%  plot log eigenvalues,
%     passing a line through those from 5 to 12

tempharmeigval = temppcastr.values;
x = ones(8,2);
x(:,2) = reshape((5:12),[8,1]);
y = log10(tempharmeigval(5:12));
c = x\y;
subplot(1,1,1)
plot(1:12,log10(tempharmeigval(1:12)),'-o', 1:12, c(1)+ c(2).*(1:12), ':')
xlabel('Eigenvalue Number')
ylabel('Log10 Eigenvalue')

%  plot pca scores

subplot(1,1,1)
tempharmscr = temppcastr.harmscr;
plot(tempharmscr(:,1),tempharmscr(:,2), 'o')
xlabel('Scores on Harmonic 1')
ylabel('Scores on Harmonic 2')
text(tempharmscr(:,1),tempharmscr(:,2),meteonames)

print -dpsc2 'c:/Myfiles/talks/fdacourse/figs/temppcascores.ps'

%  plot VARIMAX rotated PCA solution

temppcastr = varmx_pca(temppcastr);

tempharm = temppcastr.harmfd;
tempvarprop = temppcastr.varprop;
tempharmfine = eval_fd(tempharm,monthfine);
tempharmmat  = eval_fd(tempharm,monthtime);
tempharmfine(:,3) = -tempharmfine(:,3);
tempharmmat(:,3)  = -tempharmmat(:,3);

con = 5.*ones(1,4);
tempmeanfine = eval_fd(tempmeanfd,monthfine);
tempmeanvec  = eval_fd(tempmeanfd,monthtime);

for j=1:4
    subplot(2,2,j)
    yplus = tempmeanvec + con(j).*tempharmmat(:,j);
    ymins = tempmeanvec - con(j).*tempharmmat(:,j);
    plot(monthfine,tempmeanfine,'b-',...
         monthtime,yplus,'g--',...
         monthtime,ymins,'r--')
    xlabel('Month')
    ylabel('Component')
    title(['Rot. PC ',num2str(j),' (',num2str(round(tempvarprop(j)*1000)/10),'%)'])
    axis([0,12,-20,20])
end

%  ----------------------------------------------------------------------
%  Carry out some linear models.  These analyses use the older function
%  LINMOD, now replaced by the more powerful function FREGRESS for all
%  but functional-to-functional regressions.  To see FREGRESS in action,
%  do the analyses in file daily.m.
%  ----------------------------------------------------------------------

%  ----------------------------------------------------------------------
%          Linear model for temperature using climate zones  
%          This is a functional one-way analysis of variance
%  ----------------------------------------------------------------------

%  setup design matrix.  Make it full rank and use atlantic zone
%     as baseline group

region_names = [ ...
'Canada  '; 'Atlantic'; 'Pacific '; 'Contintl'; 'Arctic  '];

zmat = zeros(35,4);
zmat(:       ,1) = 1;
zmat(pacindex,2) = 1;
zmat(conindex,3) = 1;
zmat(artindex,4) = 1;

p = 4;
xfdcell = cell(1,p);
for j=1:p
    xfdcell{j} = zmat(:,j);
end

%  set up the basis for the regression functions

nbetabasis = 13;
betabasis  = create_fourier_basis([0,12], nbetabasis);

%  set up the functional parameter object for the regression fns.

betafd    = fd(zeros(nbetabasis,p), betabasis);
estimate  = 1;
lambda    = 0;
betafdPar = fdPar(betafd, harmaccelLfd, lambda, estimate);
betacell = cell(1,p);
for j=1:p
    betacell{j} = betafdPar;
end

%  estimate linear model

fRegressStruct = fRegress(tempfd, xfdcell, betacell);

betaestcell = fRegressStruct.betahat;
tempyhatfd  = fRegressStruct.yhat;

%  plot four regression functions

for j=1:p
    plot(getfd(betaestcell{j}))
    title(['\fontsize{16} ',region_names(j+1,:)])
    pause
end

%  plot approximation functions (there are only 4 distinct ones)

plot(tempyhatfd)
title('Model Functions')

%  compute residual functions

tempresfd = tempfd - tempyhatfd;

%  compute mean of temperatures

tempmeanfd = mean(tempfd);

%  Compute the error sum of squares function

tempresmat = eval_fd(tempresfd, monthtime);
SSE = sum((tempresmat').^2)';

subplot(1,1,1)
plot(monthtime, SSE)
xlabel('Month')
ylabel('SSE')

%  Compute error sum of squares about mean

tempresmat0 = eval_fd(tempfd, monthtime) - ...
              eval_fd(tempmeanfd, monthtime)*ones(1,35);
SSY = sum((tempresmat0').^2)';

%  Compute squared multiple correlaton and F-ratio functions

RSQ = (SSY - SSE)./SSY;
Fratio = ((SSY - SSE)./3)./(SSE./31);

%  Plot these functions, 
%    along with 0.05 critical value for F for 3 and 31 df

subplot(1,2,1)
plot(monthtime, RSQ)
xlabel('Month')
title('R^2')
axis([0,12,0,1])
axis('square')
subplot(1,2,2)
plot(monthtime, Fratio, '-', [0,12], [2.9,2.9], '--')
xlabel('Month')
title('F')
axis([0,12,0,40])
axis('square')

%  Set up a design matrix having a column for the grand mean, and
%    a column for each climate zone effect. Add a dummy contraint
%    observation

zmat = zeros(35,5);
zmat(:       ,1) = 1;
zmat(atlindex,2) = 1;
zmat(pacindex,3) = 1;
zmat(conindex,4) = 1;
zmat(artindex,5) = 1;
%  attach a row of 0, 1, 1, 1, 1
z36    = ones(1,5);
z36(1) = 0;
zmat = [zmat; z36];

%  set up a new fd by adding a zero function

coef = getcoef(tempfd);
coef36 = [coef,zeros(13,1)];
tempfd36 = putcoef(tempfd, coef36);

p = 5;
xfdcell = cell(1,p);
for j=1:p
    xfdcell{j} = zmat(:,j);
end

betafd    = fd(zeros(nbetabasis,p), betabasis);
estimate  = 1;
lambda    = 0;
betafdPar = fdPar(betafd, harmaccelLfd, lambda, estimate);
betacell = cell(1,p);
for j=1:p
    betacell{j} = betafdPar;
end

%  fit linear model

fRegressStruct = fRegress(tempfd36, xfdcell, betacell);

betaestcell = fRegressStruct.betahat;
tempyhatfd  = fRegressStruct.yhat;

%  plot five regression functions

for j=1:p
    plot(getfd(betaestcell{j}))
    title(['\fontsize{16} ',region_names(j,:)])
    pause
end

%  ----------------------------------------------------------------------
%       now model log mean precipitation as a function of temperature  
%  ----------------------------------------------------------------------

%  Compute log precipitation

logannprec = log10(sum(precmat)');
logannprec = logannprec - mean(logannprec);

%  set up the basis for the regression functions

nbetabasis = 5;
betabasis  = create_fourier_basis([0,12], nbetabasis);
estimate  = 1;
lambda    = 0;
betafdPar = fdPar(betabasis, harmaccelLfd, lambda, estimate);

p = 2;
betacell = cell(1,p);
constbasis = create_constant_basis([0,12]);
betacell{1} = fdPar(constbasis);
betacell{2} = betafdPar;

tempcell{1} = ones(35,1);
tempcell{2} = tempfd;

fRegressStruct = fRegress(logannprec, tempcell, betacell);

betaestcell = fRegressStruct.betahat;
logannprechat  = fRegressStruct.yhat;

%  plot five regression functions

for j=1:p
    plot(getfd(betaestcell{j}))
    pause
end

%  compute fitted values and plot Y against Yhat

plot(logannprechat, logannprec, 'o', ...
     logannprechat, logannprechat, '--')
xlabel('Predicted log annual precipitation')
ylabel('Log annual precipitation')

%  compute squared correlation

covmat = cov(logannprechat,logannprec);
RSQ    = covmat(1,2)^2/(covmat(1,1)*covmat(2,2))

%  -------------------------------------------------------------------
%              Predict log precipitation from temperature   
%   The regression coefficient function is now bivariate:  \beta(s,t)
%  -------------------------------------------------------------------

%  set up functional data object for log precipitation

%  Model is y(t) = \alpha(t) + \int_0^T x((s)( \beta(s,t) ds + e(t)

%  set up smoothed log precipitation: the dependent variable

lnprecfd = smooth_basis(monthtime, log10(precmat), monthfdPar);

lnprecfd_fdnames{1} = 'Months';
lnprecfd_fdnames{2} = 'Station';
lnprecfd_fdnames{3} = 'log_{10} mm';
lnprecfd = putnames(lnprecfd, lnprecfd_fdnames);

%  plot log precipitation functions

plot(lnprecfd);
title('Log Precipitation Functions')

%  set up bases for variation over s and t

nsbasis = 5;
ntbasis = 5;
sbasis  = create_fourier_basis([0,12], nsbasis);
tbasis  = create_fourier_basis([0,12], ntbasis);

%  set up the cell array containing :
%    ... the fdPar   object for the estimate of \alpha(t)
%    ... the bifdPar object for the estimate of \beta(s,t)

betacell = cell(1,2);
%  fdpar object for \alpha
lambda0   = 1e-8;
alphafdPar = fdPar(betabasis, harmaccelLfd, lambda0);
%  bifdPar object for \beta
betabifd = bifd(zeros(13,13), betabasis, betabasis);
lambdas   = 1e-4;
lambdat   = 1e-4;
betabifdPar = bifdPar(betabifd, ...
                      harmaccelLfd, harmaccelLfd, ...
                      lambdas, lambdat);
%  place these two objects in beta cell array
betacell{1} = alphafdPar;
betacell{2} = betabifdPar;

%  compute the linear model:  tempfd is x, lnprecfd is y

linmodStruct = linmod(tempfd, lnprecfd, betacell);

%  retrieve the two estimated coefficient functions

alphafd  = linmodStruct.alpha;   %  The intercept function
betabifd = linmodStruct.beta;    %  The bivariate regression function

%  plot the intercept function estimate

plot(alphafd);
xlabel('\fontsize{13} Month (t)')
ylabel('\fontsize{13} \alpha(s)')

%  plot the regression function estimate as a surface

betamat = eval_bifd(weeks, weeks, betabifd);

subplot(1,1,1)
surf(weeks, weeks, betamat)
xlabel('\fontsize{13} Month (t)')
ylabel('\fontsize{13} Month (s)')
title('Regression Function \beta(s,t)')
axis([0,12,0,12,-0.4,0.4])

%  plot the regression function as a function of s 
%    for selected values of t.  Press any key to advance plot

for t=1:2:53
    plot(weeks, betamat(:,t), '-', [0,12], [0,0], '--')
    xlabel('\fontsize{13} Month (s)')
    ylabel('\fontsize{13} \beta(s,t)')
    title(['\fontsize{16} t = ',num2str(weeks(t))])
    axis([0,12,-0.4,0.4])
    pause
end

%  Get fitted functions

lnprechatfd = linmodStruct.yhat;

% Compute mean function as a benchmark for comparison

lnprecmeanfd = mean(lnprecfd);

%  Plot actual observed, fitted, and mean log precipitation for
%      each weather station, 

lnprechat0 = eval_fd(lnprecmeanfd,weeks);
for i=1:35
    lnpreci    = eval_fd(lnprecfd(i),    weeks);
    lnprechati = eval_fd(lnprechatfd(i), weeks);
    SSE = sum((lnpreci-lnprechati).^2);
    SSY = sum((lnpreci-lnprechat0).^2);
    RSQ = (SSY-SSE)/SSY;
    plot(weeks, lnpreci, 'o', weeks, lnprechati, 'b-', ...
                              weeks, lnprechat0, 'b--')
    xlabel('\fontsize{13} Month')
    ylabel('\fontsize{13} Log Precipitation')
    title(['\fontsize{16} ',meteonames(i,:),'  R^2 = ',num2str(RSQ)])
    legend('\fontsize{13} Data',      ...
           '\fontsize{13} Model fit', ...
           '\fontsize{13} Mean',      ...
           'Location', 'NorthWest')
    pause
end

%  ---------------------------------------------------------------
%                   Register temperature data   
%  ---------------------------------------------------------------

%  register the first derivative of temperature

%  first smooth the temperature functions 

harmaccelfdPar = fdPar(tempfd, harmaccelLfd, 1e-3);
stempfd = smooth_fd(tempfd, harmaccelfdPar);

%  set up the basis for the warping function

nbasis = 5;
wbasis = create_fourier_basis([0,12],nbasis);

%  set up parameters for the registration function

Lfd      = int2Lfd(3);
lambda   = 1e-3;

index = 1:35;
Dtempfd = deriv(stempfd(index),1);
y0fd = mean(Dtempfd);
index = [1, 8, 24, 27, 32];  %  register a subset of the stations
yfd  = Dtempfd(index);
xfine = linspace(0,12,101)';
ofine = ones(101,1);
y0vec = eval_fd(y0fd, xfine);
yvec  = eval_fd(yfd, xfine);

cvec0    = zeros(nbasis,length(index));
Wfd0     = fd(cvec0, wbasis);
WfdPar   = fdPar(Wfd0, Lfd, lambda);
periodic = 1;

[yregfd, Wfd, shift] = registerfd(y0fd, yfd, WfdPar, periodic);

yregmat = eval_fd(yregfd,xfine);
warpmat = monfn(xfine, Wfd);
warpmat = ofine*shift' + 12.*warpmat./(ofine*warpmat(101,:));

for i = 1:length(index)
   subplot(1,2,1)
   plot(xfine, yvec(:,i),    '-',  ...
        xfine, y0vec,        '--', ...
        xfine, yregmat(:,i), '-');
   axis('square')
   title(meteonames(index(i),:))
   subplot(1,2,2)
   plot(xfine, warpmat(:,i),   '-', ...
        xfine, xfine+shift(i), '--')
   axis('square')
   title(['Shift = ',num2str(shift(i))])
   pause
end




