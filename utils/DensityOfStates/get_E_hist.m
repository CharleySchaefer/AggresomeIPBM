% Get energy histogram.
% Emat can be a 1,2 or 3D array
function [Erow, prow]=get_E_hist(Emat, Nbins)
  [NX,NY,NZ]=size(Emat);

  Emin=min(Emat(:));
  Emax=max(Emat(:));
  % Get margins
  LE=Emax-Emin;
  Emin=Emin-0.1*LE;
  Emax=Emax+0.1*LE;
  Erow=linspace(Emin, Emax, Nbins);
  dE=Erow(2)-Erow(1);
  prow=zeros(size(Erow));
  for i=1:NX
    for j=1:NY
      for k=1:NZ
        ind=round( (Emat(i,j,k)-Emin)/dE  );
        prow(ind)=prow(ind)+1; 
      end
    end
  end
  prow=prow/sum(prow(:))/dE;
end
