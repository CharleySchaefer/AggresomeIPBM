function getRMSDFromParticleTracking(pth, dx, DA, label)
  clc; close all;
  Nheader=1;
%  dx=20; % nm
%  DA=0.3; % micron^2/s
  delim='\t';
  data=importdata(sprintf('%s/ParticleTracking%c.out', pth, label), delim, Nheader);
  data=data.data;
  [Ntime, nA3]=size(data);
  nA=(nA3-1)/3; % Number of particles

  %trow=data(:,1);
  pos0=data(1,2:end);     % Initial position
  dpos=zeros(size(pos0)); % displacements
  %MSD=zeros(1,Ntime);
  %msd;
  timefac=dx*dx/DA/1000000; % Conversion kMC time to real time
  fp=fopen(sprintf('%s/rmsd.out', pth), 'w');
    fprintf(fp, '%12s %12s\n', 'time(s)', 'rmsd(micron)');
  for i=2:Ntime
    tt=data(i,1);            % time
    dpos=data(i,2:end)-pos0; % displacement
    msd=dpos*dpos'/nA;       % mean square displacement (averaged over all particles)
    %MSD(i)=dpos*dpos'/(nA3/3);
    fprintf(fp, '%12e %12e\n', timefac*tt, (dx/1000)*sqrt(msd));
  end
  fclose(fp);


#  figure
#  loglog(timefac*trow(2:end), (dx/1000)*sqrt(MSD(2:end))); hold on
#  loglog(timefac*trow(2:end), sqrt(4*DA*timefac*trow(2:end)))
#  xlabel('time [s]')
#  ylabel('RMSD [micron]')
#  legend('simulation (with interaction)', 'theory (without interactions)', 'Location', 'northwest')

end
