function plot_pca_mult(pcastr, nx, harm, expand)
%  PLOT_PCA_MULT  Plots the harmonics for a principal components analysis.
%  It differs from PLOT_PCA in plotting these as multiple panels without
%    prompting.  It only accepts univariate functions.
%  Arguments:
%  PCASTR    ... Struct object returned by PCA_FD.
%  NX        ... Number of argument values for plotting. Default = 101.
%  HARM      ... If harm = 0 (the default) then all the computed harmonics
%                are plotted.   Otherwise those in HARM are plotted.
%  EXPAND    ... If expand =0 then effect of +/- 2 standard deviations of
%                each harmonic are given.
%                Otherwise the factor expand is used.

%  last modified 20 July 2006

  %  set up default argument values

  if nargin < 4
    expand = 0;
  end
  if nargin < 3
    harm = 0;
  end
  if nargin < 2
    nx = 51;
  end

  harmfd  = pcastr.harmfd;
  basis   = getbasis(harmfd);
  rangex  = getbasisrange(basis);
  x       = linspace(rangex(1), rangex(2), nx);
  fdmat   = eval_fd(harmfd, x);
  meanmat = squeeze(eval_fd(pcastr.meanfd, x));
  dimfd   = size(fdmat);
  nharm   = dimfd(2);
  if harm == 0
    harm = (1:nharm);
  end
  onesharm = ones(1,nharm);
  if expand == 0
      fac = ones(nx,1)*sqrt(pcastr.eigvals(1:nharm))';
  else
      fac = expand;
  end 
  meanplus  = meanmat*onesharm + fac.*fdmat;
  meanminus = meanmat*onesharm - fac.*fdmat;
  plottop = max(max([meanplus;meanminus]));
  plotbot = min(min([meanplus;meanminus]));
  if length(dimfd) == 2
    %  plotting for univariate functions
    for iharm = harm
      percentvar = round(100 * pcastr.varprop(iharm));
      if nharm < 4
          subplot(1,nharm,iharm)
      else
          subplot(2,2,iharm)
      end
      if plottop*plotbot < 0
          plot(x, meanmat, '-', [min(x),max(x)], [0,0], '--')
      else
          plot(x, meanmat, '-')
      end
      text(x-.2, meanplus(:,iharm),  '+')
      text(x-.2, meanminus(:,iharm), '-')
      axis([rangex(1), rangex(2), plotbot, plottop])
      if ~(nharm == 4), axis('square'); end
      title(['\fontsize{12} Component ', num2str(iharm), '  ', ...
              num2str(percentvar), '%'])
    end
  else
    error('Multivariate functions cannot be plotted.');
  end

