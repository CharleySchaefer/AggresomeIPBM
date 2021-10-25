function createDensityOfStates()
  clc; close all;
  [filepath,scriptname,ext] = fileparts(mfilename('fullpath'));
  addpath(sprintf('%s/DensityOfStates',filepath));

  CellTopology=importdata('TopologyTwoAggresomesNX50NY150.dat');
  [NX, NY]=size(CellTopology);
  NZ=1;
  nA=150; % number of proteins

  % Gaussian density of states in cytosol
  Emean0=2.0;
  sigma0=0.4;

  % Gaussian density of states in aggresome
  BindingEnergy=2.0;
  Emean1=Emean0-BindingEnergy;
  sigma1=0.4;

  % Estimated concentration difference (fraction)
  cfrac=exp(BindingEnergy) 

  % Number of proteins in phase 0 (cytosol) and 1 (aggresome) 
    % cfrac = c1/c0 = (n1/V1)/(n0/V0), nA=n0+n1
    % --> cfrac = (n1/V1)/(( nA-n1)/V0)
    % --> cfrac*(( nA-n1)/V0) = n1/V1
    % --> n1*(1.0/V1+cfrac/V0) = cfrac*nA/V0
  ctot=nA./(NX*NY*NZ);
  V0=sum(CellTopology(:)==0) % Cytosol volume
  V1=sum(CellTopology(:)==1) % Aggresome volume
  n1=round( cfrac*nA./(V0*1.0/V1+cfrac) );
  n0=nA-n1;

  % Create energy landscape
  Emat=zeros(NX,NY,NZ);
  for i=1:NX
    for j=1:NY
      for k=1:NZ
        if CellTopology(i,j,k)==1 % Aggresome
          Emat(i,j,k) = dos_gaussian_cumm_inv( rand(), Emean1, sigma1 ) ;
        elseif CellTopology(i,j,k)==0 % cytosol
          Emat(i,j,k) = dos_gaussian_cumm_inv( rand(), Emean0, sigma0 ) ;
        end
      end
    end
  end

  % Get Fermi levels in cytosol and aggresome
  maxiter=1000;
  cnf=zeros(NX,NY,NZ);
  EF0=getFermiLevel( Emat(CellTopology==0) , n0, maxiter)
  EF1=getFermiLevel( Emat(CellTopology==1) , n1, maxiter);

  cnf=( ((CellTopology==0).*(Emat<=EF0)) + ...
            ((CellTopology==1).*(Emat<=EF1))  );
  
  % Histogram DOS
  Nbins=200; % Number of bins
  [Erow, prow]=get_E_hist(Emat, Nbins);
  % Histogram DOOS
  Edoos=Emat(find(cnf==1));
  [Erowdoos, prowdoos]=get_E_hist(Edoos, Nbins);
  Conc=nA/(NX*NY);

  figure  
  subplot(2,2,1)
  plot(Erow,prow); hold on
 % plot(E,p, 'LineWidth', 2)
  xlabel('E/kT')
  ylabel('DOS(E)')
  title('Density of States')
  title('Density of occupied States')
  subplot(2,2,2)
  plot(Erowdoos,prowdoos); hold on
  plot(Erow,exp(-Erow).*prow, 'LineWidth', 2); hold on
  plot(Erow,prow./(exp(Erow-2.0)+1), 'LineWidth', 2); hold on
  %plot(E,p/Conc, 'LineWidth', 2)
  xlabel('E/kT')
  ylabel('DOOS(E)')
  title('Density of occupied States')
  subplot(2,2,3)
  imagesc(Emat)
  title('Energy landscape')
  subplot(2,2,4)
  imagesc(cnf)
  title('Protein configuration')

if 0
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
  EF=getFermiLevel(Emat, nA, maxiter)
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
  imagesc(cnf)
  title('Protein configuration')
  
end

end



