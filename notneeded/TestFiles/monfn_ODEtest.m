n  = 201;
x  = linspace(0,1,n)';
y0 = exp(5*x) - 1;
% sig = 0.05;
% y = y + randn(n, 1).*sig;

% breaks = linspace(0,1,11);
% norder = 6;
% nbasis = length(breaks) + norder - 2;
% wbasis = create_bspline_basis([0,1],nbasis,norder,breaks);
% wfdPar = fdPar(wbasis, 1, 1e-5);
% Wfd = smooth_monotone(x, y, wfdPar);

wbasis = create_monomial_basis([0,1],2);
Wfd = fd([0;5],wbasis);

% EPS = 1e-4; 
% odeoptions = odeset('RelTol', EPS, 'AbsTol', EPS, ...
%                     'MaxStep',     5e-1,        ...
%                     'InitialStep', 1e-1);
% 
% [xval, hval] = ode23(@Dyfn, [0,1], 0, odeoptions, Wfd);

tic
y1 = monfn_ODE(x, Wfd);
toc
y1 = y1*y0(n)/y1(n);
max(abs(y1-y0))

tic
y2 = monfn(x, Wfd)/5;
toc
y2 = y2*y0(n)/y2(n);
max(abs(y2-y0))

plot(x, z, '-o')

max(abs(z1-z2))

tic
y1 = mongrad_ODE(x, Wfd);
toc

tic
y2 = mongrad(x, Wfd)/5;
toc

plot(x, z, '-o')

max(abs(z1-z2))

% global Wfd
% 
% t=0;   % initial time
% y=0;   % initial conditions
% fun=@Dyfn;        % create function handle to ODEFUN
% yp=feval(fun,t,y);  % compute initial slopes
% odess('Reset')      % reset parameters if new problem
% %                 odess('PName',Pvalue,...) % optionally set parameters
% test = 1;
% while test
%     disp(t)
%     odess('MaxStep', 1-t)
%     [t,y,yp]=odess(fun,t,y,yp);
%     % process data as needed
%     % modify ode parameters if needed
%     % do not change t,y,yp from call to call!
%     % break from loop when desired
%     test = t < 1;
% end
