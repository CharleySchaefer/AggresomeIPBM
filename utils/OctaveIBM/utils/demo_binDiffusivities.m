function demo_binDiffusivities()
  clc; close all;

  %-------------------------------
  % Generate artifical data set
  Ndata=1000;
  Diffusivities=zeros(1,Ndata);
  for i=1:Ndata
      r= exprnd(2.9);
      Diffusivities(i)=exp(-r)*0.23;
  end

  %-------------------------------
  % Do binning in log space
  Nbins=200;
 [logD_200, counts_200]=binDiffusivities(Diffusivities, Nbins);
  Nbins=20;
 [logD_20, counts_20]=binDiffusivities(Diffusivities, Nbins);

  %-------------------------------
  % Plot result
 figure
 plot(logD_200, counts_200); hold on;
 plot(logD_20, counts_20/10); hold on;
 xlabel('log10(D)')
 ylabel('counts')
 legend('200 bins', '20 bins (counts divided by 10)', 'Location', 'northwest')
end

% Input:  Diffusivities (can be 1D or 2D array)
%         Nbins     (integer value): number of bins
% Output: logDspace (1D array):      log10 of the diffusivity
%         counts    (1D array):      Number of diffusivities
function [logDspace, counts]=binDiffusivities(Diffusivities, Nbins)
  [NX,NY]=size(Diffusivities);

  % Create bins
  logDL=log10(min(Diffusivities(:)))-1;
  logDU=log10(max(Diffusivities(:)))+1;
  logDspace=linspace(logDL, logDU,Nbins );
  dd=logDspace(2)-logDspace(1); % bin spacing in log space

  % Counter
  counts=zeros(1,Nbins);
  for i=1:NX ;    for j=1:NY
    logD=log10(Diffusivities(i,j));
    ind=round((logD-logDL)/dd);
    counts(ind)=counts(ind)+1;
  end; 
end;
end


