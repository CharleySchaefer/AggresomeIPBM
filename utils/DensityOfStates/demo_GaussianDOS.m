function demo_GaussianDOS()
  clc; close all;
  %[filepath,scriptname,ext] = fileparts(mfilename('fullpath'));
  %addpath(sprintf('%s/DensityOfStates',filepath));

  nA=150; % number of proteins
  NX=50;
  NY=50;
  Emean=2.0;
  sigma=0.4;
  NE=100;


  % Sample Gaussian Density of States (DOS)
  pmat=rand(NX,NY);
  Emat=dos_gaussian_cumm_inv( pmat, Emean, sigma ) ;

  % Verify: compare distribution to expectation values 
  Nbins=200; % Number of bins
  [Erow, prow]=get_E_hist(Emat, Nbins);
      % Exact distribution
  E=linspace( Emean-4*sigma, Emean+4*sigma, NE );
  p=dos_gaussian(E, Emean, sigma );

  % fill low-energy states 
  maxiter=1000;
  EF=getFermiLevel(Emat, nA, maxiter); % Fermi level
  cnf=(Emat<EF);
  Edoos=Emat(find(cnf==1));
  [Erowdoos, prowdoos]=get_E_hist(Edoos, Nbins);
  Conc=nA/(NX*NY);

  figure
  subplot(2,2,1)
  plot(Erow,prow); hold on
  plot(E,p, 'LineWidth', 2)
  xlabel('E/kT')
  ylabel('DOS(E)')
  title('Density of States')
  subplot(2,2,2)
  plot(Erowdoos,prowdoos); hold on
  plot(E,p/Conc, 'LineWidth', 2)
  xlabel('E/kT')
  ylabel('DOOS(E)')
  title('Density of occupied States')
  subplot(2,2,3)
  imagesc(Emat)
  title('Energy landscape')
  subplot(2,2,4)
  imagesc(cnf+(Emat-Emean)*0.1)
  title('Protein configuration')
  


end

% Uses bisection method
function EF=getFermiLevel(Emat, nAtarget, maxiter)
  EFl=min(Emat(:));
  EFu=max(Emat(:));
  iter=0;
  while iter<maxiter
    iter=iter+1;
    EF=0.5*(EFl + EFu);
    nA=sum( Emat(:)<EF );
    if nA>nAtarget
      EFu=EF;
    elseif nA<nAtarget
      EFl=EF;
    else
     break;
    end
  end
  if iter==maxiter
    fprintf('Warning: getFermiLevel() has not converged.\n');
  end
end


