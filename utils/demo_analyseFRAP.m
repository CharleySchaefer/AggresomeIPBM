function demo_analyseFRAP()
  clc; close all;
  dname="demo2D_topology";

  analyseFRAP(dname); 

  frapdata=importdata(sprintf('%s/frap.out', dname), ' ', 1);
  texp=frapdata.data(:,1);
  n1=frapdata.data(:,2);
  n2=frapdata.data(:,3);

  ntotmean=1;
  texp=texp(3:end);
  n1=n1(3:end);
  n2=n2(3:end);
  ntot=n1+n2;

  figure
  subplot(2,1,1)
  plot(texp, n1, 'LineWidth', 2, 'k'); hold on
  plot(texp, n2, 'LineWidth', 2, '--k'); hold on
  plot(texp, ntot, 'LineWidth', 1, '.-k');
  plot(texp, ones(1,length(texp)), 'LineWidth', 1, '--r');
  xlabel('time[kmc units]')
  ylabel('FRAP')
  subplot(2,1,2)
  loglog(texp(n1>0), n1(n1>0), 'LineWidth', 2, 'k'); hold on
  loglog(texp(n1>0), (texp(n1>0)).^0.5/2000, 'LineWidth', 2, 'r'); hold on
  xlabel('time[kmc units]')
  ylabel('FRAP')

end 
