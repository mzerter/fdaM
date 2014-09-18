rng = [0,4*pi];
argvals = linspace(0,4*pi,51)';
y0 = argvals + sin(argvals);
sig = 0.5;

ncurve = 5;
y = y0*ones(1,ncurve) + randn(51,ncurve).*sig;

y1 = y(:,1:2);
y2 = y(:,2:3);
y = zeros(51,2,2);
y(:,:,1) = y1;
y(:,:,2) = y2;
ncurve = 2;
nvar   = 2;

wbasis = create_bspline_basis(rng, 53);

lambda = 10^(-1);

Wfd0 = smooth_basis(argvals, log(y0*ones(1,ncurve)+0.1), ...
                    fdPar(wbasis,2,lambda));
% Wfd0 = fd(zeros(53,ncurve,nvar), wbasis);
fdParobj = fdPar(Wfd0, 2, lambda);

[Wfd, beta] = smooth_monotone(argvals, y, fdParobj);

[Wfd] = smooth_pos(argvals, y, fdParobj);

% disp(beta)

disp(Fstr)

yfit = zeros(51,ncurve);
for icurve=1:ncurve
%     yfit(:,icurve) = ...
%         beta(1,icurve) + beta(2,icurve).*monfn(argvals, Wfd(icurve));
    yfit(:,icurve) = ...
        exp(eval_fd(argvals, Wfd(icurve)));
end

rmse = sqrt(mean((y - yfit).^2));

disp(rmse)

subplot(1,1,1)
for icurve=1:ncurve
    plot(argvals, y(:,icurve), 'o', ...
        argvals, yfit(:,icurve), 'b-', ...
        argvals, y0, 'g--')
    title(['Curve ',num2str(icurve)])
    pause
end
