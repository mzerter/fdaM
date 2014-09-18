%  Simulating data for Zephirin Yepdo and illustrating the analysis of them
%  Jim Ramsay, 26 July 2010

%  First, load some sample precipitation and temperature anomalies that
%  he has provided

load YepdoPrec.txt
load YepdoTemp.txt

%  Compute the means of these data for each of four regions

YepdoPrecMean = mean(YepdoPrec(:,2:5));
YepdoTempMean = mean(YepdoTemp(:,2:5));

disp(['Precipitation means:  ',num2str(YepdoPrecMean)])
disp(['Temperature   means:  ',num2str(YepdoTempMean)])

% Precipitation means:  -0.0202    129.2193    115.7105    103.2208
% Temperature   means:  -2.0217     -4.8508      1.2801     -1.9313

%  Compute the standard deviation of these data for each of four regions

YepdoPrecStdDev = sqrt(var(YepdoPrec(:,2:5)));
YepdoTempStdDev = sqrt(var(YepdoTemp(:,2:5)));

disp(['Precipitation Std. Dev.:  ',num2str(YepdoPrecStdDev)])
disp(['Temperature   Std. Dev.:  ',num2str(YepdoTempStdDev)])

% Precipitation Std. Dev.:  10.3108    217.3568     89.6749     89.1387
% Temperature   Std. Dev.:   3.3201      5.6122      2.0821      1.9676

%  Now we use these means and standard deviations to generate some 
%  random Gaussian data over 40 years for each of 12 months and 
%  four regions.  The data are stored in a matrix with
%  dimensions 'Month' and, 'Year/Region'  Year varies inside regions.

%  Precipitation data

PrecData = zeros(12,40*4);
m = 0;
for iregion = 1:4
    for jyear = 1:40
        m = m + 1;
        PrecData(:,m) = YepdoPrecMean(iregion) + ...
                       randn(12,1).*YepdoPrecStdDev(iregion);
    end
end

%  Temperature data

TempData = zeros(12,40*4);
m = 0;
for iregion = 1:4
    for jyear = 1:40
        m = m + 1;
        TempData(:,m) = YepdoTempMean(iregion) + ...
                       randn(12,1).*YepdoTempStdDev(iregion);
    end
end

%  Set up some time points within the year interval [0.5, 11.5]

Month = linspace(0.5,11.5,12)';

%  Now these data are converted to functional data objects using 
%  12 piecewise linear basis functions.  This interpolates the data,
%  and does not therefore lose any information in the original data

%  Set up the polygonal spline basis object

polygBasis = create_polygon_basis(Month);

% plot the basis

plot(polygBasis)  

%  Smooth the data using this basis

%  Precipitation

Precfd = smooth_basis(Month, PrecData, polygBasis);

%  Temperature

Tempfd = smooth_basis(Month, TempData, polygBasis);

%  Plot the data.  There will be 40 times 4 plots to look at, with
%  regions varying inside years

%  Plot precipitation

plotfit_fd(TempData, Month, Precfd)

%  Plot temperature

plotfit_fd(TempData, Month, Tempfd)

%  Now prepare to carry out the regression analysis.

%  First set up a design matrix for the data over four regions.  
%  The simplest method is to choose one of the regions as a baseline
%  region, and make an indicator variable for the remaining three
%  stations.  Here we choose the first region as the baseline.
%  The design matrix will have 40 times 4 rows and 4 columns.
%  The first column indicates the mean effect, and has 1 for all values.
%  The second codes the effect of being in the second region, and has
%  0's for the first 12 entries, 1's for the second 12, and 0's for the
%  remainder.  And so on for the third and fourth columns.

Z = zeros(40*4,4);
Z(   :   ,1) = 1;
Z( 41: 80,2) = 1;
Z( 81:120,3) = 1;
Z(121:160,4) = 1;

%  The model in this case is
%  Prec(t) = \mu(t) + \alpha_j(t) + \beta(t) Temp(t) + e(t)
%  where j = 2,3 and 4.

%  Now set up basis objects for each of the five regression functions
%  I've used five B-spline basis functions for each function.

betacell = cell(1,5);
betacell{1} = create_bspline_basis([0.5,11.5], 5);
betacell{2} = create_bspline_basis([0.5,11.5], 5);
betacell{3} = create_bspline_basis([0.5,11.5], 5);
betacell{4} = create_bspline_basis([0.5,11.5], 5);
betacell{5} = create_bspline_basis([0.5,11.5], 5);

%  Also set up a cell array for the independent variables: there are
%  two:  The design matrix (multivariate) and temperature (functional)

xfdcell = cell(1,5);
xfdcell{1} = Z(:,1);
xfdcell{2} = Z(:,2);
xfdcell{3} = Z(:,3);
xfdcell{4} = Z(:,4);
xfdcell{5} = Tempfd;

%  carry out the regression analysis

fRegressStruct = fRegress(Precfd, xfdcell, betacell);

%  Select the cell array containing the estimated regression coefficients

betaestcell = fRegressStruct.betahat;

%  plot each coefficient

figure(1)
plot(getfd(betaestcell{1}))
xlabel('Month')
ylabel('\mu(t)')
title('Mean  effect \mu(t)')

figure(2)
plot(getfd(betaestcell{2}))
xlabel('Month')
ylabel('\alpha_2(t)')
title('Region 2 effect \alpha_2(t)')

figure(3)
plot(getfd(betaestcell{3}))
xlabel('Month')
ylabel('\alpha_3(t)')
title('Region 3 effect \alpha_3(t)')

figure(4)
plot(getfd(betaestcell{4}))
xlabel('Month')
ylabel('\alpha_4(t)')
title('Region 4 effect \alpha_4(t)')

figure(5)
plot(getfd(betaestcell{5}))
xlabel('Month')
ylabel('\beta(t)')
title('Temperature effect \beta(t)')

%  Select the functional data object for the fit

yhatfd = fRegressStruct.yhat;

plot the Precipitation and its fit for each curve

figure(6)
for i=1:160
    phdl = plot(Precfd(i));
    set(phdl, 'LineWidth', 2, 'LineStyle', '-', 'color', 'b')
    lhdl = line(yhatfd(i));
    set(lhdl, 'LineWidth', 2, 'LineStyle', '--', 'color', 'r')
    title(['Curve ',num2str(i)])
    pause
end
