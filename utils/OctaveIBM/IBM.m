function IBM()
  % Charley Schaefer, University of York, 2020.
  % Lattice simulation of A and B proteins in a cell.
  % Interactions:
  %   Excluded volume: 
  %     > Each site can only contain one A protein
  %       and one B protein
  %     > A-A nearest neighbours have an interaction
  %       energy epsAA
  %     > When an A and B protein occupy the same site
  %       they have an interaction energy epsAB
  %   The A-A interaction driven phase separation
  %   The A-B interaction gives protein B an affinity
  %   to diffuse to A-rich regions. 
  %   
 close all; clc
 addpath('src');
  %--------------------------------------------
 % % USER INPUT
  outdir='data_wholefrapB'
  Ntime=4096000;  % Number of time steps
  Nprint=40960;  % Number of time steps
  tstart=0.0;     % time t=0;

  % Box size
  Lx=1000;   % cell width in nm
  Ly=3000;   % cell length in nm
  dx=50; % resolution of the simulation in nm
  NX=Lx/dx; NY=Ly/dx; % Number of lattice sites in X and Y direction

  % Protein properties
  DA=0*0.23;     % Maximum diffusivity [micron^2/s] of A 
  DBin =0.00010;   % Maximum diffusivity [micron^2/s] of B in A phase
  DBout=0.0004;  % Maximum diffusivity [micron^2/s] of B outside A phase
  epsAA=-2.0;  % Interaction energy  [kT] between A-A nearest neighbours 
  epsAB=-2.5;  % Interaction energy [kT] between A and B at the same site
  nA=400;
  nB=320 ;    

  % Photobleaching
  tbleach=600.0;       % time [seconds] at which bleaching takes place
                     % (time needed either to develop morphology
                     %  or to equilibrate the system)  
   Xbleach=NX/2;     % center of laser focus [units of pixel size]
   Ybleach=NX/2; 
  bleach_radius=0.37; % radius [micron] of laser focus
  bleach_radius=bleach_radius/(dx*1e-3); % units of pixel size
  

  % OUTPUT
  %colourmap=[0 0 0; ;  0.8 0.1 0.1; 0.1 0.8 0.1;]; 
  colourmap=[0 0 0;  ; 0.22 0.86 0.84; 0.8 0.1 0.1];
  itarget=1; % first iteration after t=tbleach to take statistics
  itarget_multiplier=1.5; % Increase increment in time steps to take statistics
  %--------------------------------------------
  system(sprintf('mkdir -p %s', outdir));
  system(sprintf('mkdir -p %s/config', outdir));
  system(sprintf('mkdir -p %s/img', outdir));
  fname=sprintf('%s/timeprogress.dat', outdir);
  ifp=fopen(fname, 'w');
  fprintf(ifp, '#Settings\n');
  fprintf(ifp, '#Ntime %d\n', Ntime);
  fprintf(ifp, '#Nprint %d\n', Nprint);
  fprintf(ifp, '#tstart %e\n', tstart);
  fprintf(ifp, '#LX %e\n', Lx);
  fprintf(ifp, '#LY %e\n', Ly);
  fprintf(ifp, '#dx %e\n', dx);
  fprintf(ifp, '#NX %d\n', NX);
  fprintf(ifp, '#NY %d\n', NY);
  fprintf(ifp, '#DA %e\n', DA);
  fprintf(ifp, '#DBin %e\n', DBin);
  fprintf(ifp, '#DBout %e\n', DBout);
  fprintf(ifp, '#epsAA %e\n', epsAA);
  fprintf(ifp, '#epsAB %e\n', epsAB);
  fprintf(ifp, '#DBout %e\n', DBout);
  fprintf(ifp, '#nA %e\n', nA);
  fprintf(ifp, '#nB %e\n', nB);
  fprintf(ifp, '#tbleach %e\n', tbleach);
  fprintf(ifp, '#Xbleach %e\n', Xbleach);
  fprintf(ifp, '#Ybleach %e\n', Ybleach);
  fprintf(ifp, '#bleach_radius[micron] %e\n', bleach_radius*dx*1e-3);
  fprintf(ifp, '#==\n');
  fclose(ifp); %=fopen(fname, 'w');
  Nin=0;
tic
  %--------------------------------------------
 % DEFINE MOBILITIES
  nuA   =  DA/(dx*dx*1e-6);     % Hop rate of A
  rAmax=nuA;
  nuBin_max =DBin/(dx*dx*1e-6); % Hop rate of B hops inside A phase
  nuBout_max=DBout/(dx*dx*1e-6);% Hop rate of B hops outside A phase
  nuBin=zeros(NX,NY);  % Fixed position dependence
  nuBout=zeros(NX,NY); % Fixed position dependence
  rBmax=max(nuBin_max,nuBout_max);
  for i=1:NX
    for j=1:NY
      r = 0; %1; %exprnd(2);
      nuBin(i,j)=exp(-r)*nuBin_max;

      r= exprnd(2.9);
      nuBout(i,j)=exp(-r)*nuBout_max;
    end
  end
  Nbins=200;

  [DBoutspace, DBoutbins]=binMobilities(nuBout, Nbins);
  [DBinspace, DBinbins]=binMobilities(nuBin, Nbins);

  ifp=fopen(sprintf('%s/mobilities.dat', outdir), 'w');
  fprintf(ifp, '%16s %16s %16s %16s', 'Din', 'counts', 'Dout', 'counts');
  for i=1:Nbins
    fprintf(ifp,'%12e %12d %12e %12d\n', DBoutspace(i), DBoutbins(i), DBinspace(i), DBinbins(i));
  end
  fclose(ifp);

 % [DAspace, DAbins]=binMobilities(rAmax.*ones(NX,NY), Nbins);
  figure
  plot(DBoutspace, DBoutbins); hold on
  plot(DBinspace, DBinbins); hold on
 % plot(DAspace, DAbins); hold on
  %--------------------------------------------

  %--------------------------------------------
  % INITIAL CONFIGURATION
  % cnfA(x,y):
  %   -1: no A present at location x,y
  %    1:    A present at location x,y
  % cnfB(x,y):
  %   -1:         no B present at location x,y
  %    0: unbleached B present at location x,y
  %    1:   bleached B present at location x,y
  if 0  % Initially random uniform
    [cnfA, cnfB] =initialise_cnf_uniform(NX,NY,nA, nB);
  else  % For testing: A initially in 2 disk-like aggresome
    [cnfA, cnfB]=initialise_cnf_TwoAggresomes(NX,NY,epsAB);    
    nA=sum(cnfB(:)>-1);  
    nB=sum(cnfB(:)>-1);

    ifp=fopen(fname, 'a');
    fprintf(ifp, 'Fixed Aggresome Mode - Number of proteins override:\n');
    fprintf(ifp, '#nA %e\n', nA);
      fprintf(ifp, '#nB %e\n', nB);
    fclose(ifp);
  end



  % Print initial configuration  
  figure
  subplot(3,2,1)
  imagesc(linspace(0,3,NY),linspace(0,1,NX), (cnfA+1)); hold on
  daspect([1 1 1])
  colormap(colourmap)
  title('Protein A - Initial')
  subplot(3,2,2)
  imagesc(linspace(0,3,NY),linspace(0,1,NX), (2+cnfB)); hold on
  daspect([1 1 1])
  title('HSlU - Initial')


toc
fprintf('Initialisation done.\n==\n');
  %=================================================

tic
  %=================================================
  % TIME STEPS
  % Select random site
  % A: Hop in 4 directions from NXNY sites
  %    with maximum rate rAmax; --> Total A rate: RA=4*NX*NY*rAmax
  % B: Hop in 4 directions from NXNY sites
  %    with maximum rate rBmax; --> Total B rate: RB=4*NX*NY*rBmax
  % Total rate: Rtot=4*NX*NY*(rAmax+rBmax)
  % Probability of A process happing: PA=RA/Rtot
  RA=4*NX*NY*rAmax;
  RB=4*NX*NY*rBmax*2; % factor 2: diffusion A->V and V->A
  Rtot=RA+RB;
  AFRAC=RA/Rtot;
  TFAC=1.0/(Rtot);
  fprintf('Volume fraction A: %.3f\n', nA*1.0/(NX*NY));
  fprintf('Volume fraction B: %.3f\n', nB*1.0/(NX*NY));
  fprintf('Fraction of A attempts: %.3f\n',AFRAC);
  fprintf('Start time stepper.\n');
  N_A_hops_processed=0; % count number of executed A hops
  N_B_hops_processed=0; % count number of executed B hops
  tt=tstart; istart=-1; k=0; % Count statistics

  % Print header to timeprogress file
  ifp=fopen(fname, 'a');
  fprintf(ifp, '%12s %16s %12s %12s %12s %12s\n', '#time[s]', 'frap[#proteins]', 'Ahops', 'Bhops', 'accepted[%]', 'cin/cout');
  fclose(ifp);

  for iter=1:Ntime
    % BLEACH
    if istart==-1 && tt>tbleach
      istart=iter;
      % Xbleach, Xbleach,bleach_radius: units of pixels;
      [cnfB, mask, Nmask]=bleach_B_proteins(cnfB, Xbleach, Ybleach, bleach_radius);
      n_B_bleached=sum(cnfB(:)==1); % Number of B proteins
      FRAP_plateau=1.0-n_B_bleached/nB;
      
fprintf('Bleached %d B proteins; frap plateau: %.2f\n', n_B_bleached, FRAP_plateau);
      ifp=fopen(sprintf('%s/frapplateau.out', outdir), 'w');
      fprintf(ifp, '%d %f', n_B_bleached, FRAP_plateau);
      fclose(ifp);
      

  subplot(3,2,3)
  imagesc(linspace(0,3,NY),linspace(0,1,NX), (cnfA)); hold on
  daspect([1 1 1])
  title('t_{bleach}')
  subplot(3,2,4)
  imagesc(linspace(0,3,NY),linspace(0,1,NX), (2*(cnfB==0)+(cnfB==1)));
  daspect([1 1 1])
  title('t_{bleach}')
    end

    % GET STATISTICS (placed early in the loop because 
    %                 else it will be skipped by the 'continue' 
    %                 command used) 
    if iter==Ntime || mod(iter,Nprint)==0 || (istart~=-1 && mod(iter-istart,itarget)==0 )

      if istart==-1 % Still in equilibration stage
         N_Bin = sum((cnfB(:)>-1).*(cnfA(:)==1) ); % sites containing both proteins
        cin  =N_Bin./nA;
        cout =(nB-N_Bin)/(NX*NY-nA);
       ifp=fopen(fname, 'a');
        fprintf(ifp,'%12e %16s %12d %12d %12.1f %12.2e\n', tt, '0',N_A_hops_processed, N_B_hops_processed, 100.0*(N_A_hops_processed+N_B_hops_processed)/iter, cin/cout);
        fclose(ifp);
      else    % FRAP recovery
        itarget=round(itarget*itarget_multiplier);    
        k=k+1;  % Count number of data points
        trow(k)=tt;
        tmpf=0; 
        for j=1:Nmask
          tmpf=tmpf+(0==cnfB( mask(j,1),mask(j,2) ));
          Nin=Nin+(-1<cnfB( mask(j,1),mask(j,2) ));
        end
        frap(k)=tmpf;
        N_Bin = sum((cnfB(:)>-1).*(cnfA(:)==1) ); % sites containing both proteins
        cin  =N_Bin./nA;
        cout =(nB-N_Bin)/(NX*NY-nA);

        ifp=fopen(fname, 'a');
        fprintf(ifp,'%12e %16e %12d %12d %12.1f %12.2e\n', trow(k), tmpf/(Nin/k),N_A_hops_processed, N_B_hops_processed, 100.0*(N_A_hops_processed+N_B_hops_processed)/iter, cin/cout)
        fclose(ifp);

        % Export cnf
        export_cnf(outdir, k, cnfA, cnfB);
        imwrite(cnfA+1, colourmap, sprintf('%s/img/figA%04d.png', outdir, k));
        imwrite(cnfB+2, colourmap, sprintf('%s/img/figB%04d.png', outdir, k));
      end 
    end % End snapshot


    % Get time step and update time regardless 
    % of whether a process is accepted or rejected
    % according to the 'random selection method'
    dt=-TFAC*log(rand());
    tt=tt+dt; 
  
    x=1+floor( rand()*(NX) ); % ten-fold faster than randi(NX)
    y=1+floor( rand()*(NY) );
    if rand()<AFRAC 
      Aid=cnfA(x,y);
      if Aid>-1  % A at site (x,y)    --> do attempt to do an A hop
        Attempt_A_hop=1;Attempt_B_hop=0; Bid=cnfB(x,y); 
      else       % no A at site (x,y) --> no attempt
        continue; % end iteration
      end
    else
      Bid=cnfB(x,y);
  %    if Bid>-1  % B at site (x,y)    --> do attempt to do a B hop
        Attempt_B_hop=1;Attempt_A_hop=0; Aid=cnfA(x,y); 
  %    else       % no B at site (x,y) --> no attempt
  %      continue; % end iteration
  %    end
    end

    % Select a random nearest neighbour for the particle to hop to
    rnn=1+floor( rand()*4 ); %randi(4);
    switch rnn  % slightly faster than if else then
      case 1
        if x==NX ;  continue;  else
          xnn=x+1; ynn=y;
        end
      case 2
        if x==1 ;  continue;  else
          xnn=x-1; ynn=y;
       end
      case 3
        if y==NY ;  continue;  else
          xnn=x; ynn=y+1;
        end
      case 4
        if y==1 ;   continue; else
          xnn=x; ynn=y-1;
        end
    end
    %===========================================================
    % ATTEMPT TO HOP PARTICLE A TO NEAREST NEIGHBOUR
    if Attempt_A_hop % Attempt a hop of particle A
      Aid2=cnfA(xnn,ynn);
      if (Aid==-1&&Aid2==-1)||(Aid>-1&&Aid2>-1);   % neighbour already occupied -> move not allowed
        continue; % reject -> end iteration
      end
    %  if Aid2~=-1;   % neighbour already occupied -> move not allowed
    %    continue;
    %  end
      Bid2= cnfB(xnn,ynn); 

      %-----------------------------
      % Calculate rate      M=M0*exp(-DeltaE) if DeltaE>0
      %                 and M=M0              if DeltaE<=0
      %                 DeltaE=epsAA*DeltaNA + DeltaNB*epsAB
      %  Below:
      %   1. calculate DeltaNA (difference in number of nearest neighbours)
      %   2. calculate DeltaE
      %   3. get attempt frequence M0
      % 1. Calculate DeltaNA
      DeltaNA=0; % initialise counter
      switch rnn  % slightly faster than if else then
      case 1   % Increase in x
        if x>1
          DeltaNA=DeltaNA-(cnfA(x-1,y)==1);
        end
        if xnn<NX
          DeltaNA=DeltaNA+(cnfA(xnn+1,y)==1);
        end
        if y>1
          DeltaNA=DeltaNA+(cnfA(xnn,y-1)==1);
          DeltaNA=DeltaNA-(cnfA(x , y-1)==1);
        end
        if y<NY
          DeltaNA=DeltaNA+(cnfA(xnn,y+1)==1);
          DeltaNA=DeltaNA-(cnfA(x,  y+1)==1);
        end
      case 2   % Decrease in x
        if xnn>1
          DeltaNA=DeltaNA+(cnfA(xnn-1,y)==1);
        end
        if x<NX
          DeltaNA=DeltaNA-(cnfA(x+1,y)==1);
        end
        if y>1
          DeltaNA=DeltaNA+(cnfA(xnn,y-1)==1);
          DeltaNA=DeltaNA-(cnfA(x,  y-1)==1);
        end
        if y<NY
          DeltaNA=DeltaNA+(cnfA(xnn,y+1)==1);
          DeltaNA=DeltaNA-(cnfA(x,  y+1)==1);
        end
      case 3   % Increase in y
        if y>1
          DeltaNA=DeltaNA-(cnfA(x,y-1)==1);
        end
        if ynn<NY
          DeltaNA=DeltaNA+(cnfA(x,ynn+1)==1);
        end
        if x>1
          DeltaNA=DeltaNA+(cnfA(x-1,ynn)==1);
          DeltaNA=DeltaNA-(cnfA(x-1, y)==1);
        end
        if x<NX
          DeltaNA=DeltaNA+(cnfA(x+1,ynn)==1);
          DeltaNA=DeltaNA-(cnfA(x+1, y)==1);
        end
      case 4   % Decrease in y
        if ynn>1
          DeltaNA=DeltaNA+(cnfA(x,ynn-1)==1);
        end
        if y<NY
          DeltaNA=DeltaNA-(cnfA(x,y+1)==1);
        end
        if x>1
          DeltaNA=DeltaNA+(cnfA(x-1,ynn)==1);
          DeltaNA=DeltaNA-(cnfA(x-1, y)==1);
        end
        if x<NX
          DeltaNA=DeltaNA+(cnfA(x+1,ynn)==1);
          DeltaNA=DeltaNA-(cnfA(x+1, y)==1);
        end
      end % End switch -> DeltaNA calculated

      % 2. Calculate DeltaNA
        DeltaNB=(Bid2>-1)-(Bid>-1);
      if Aid>1
        DeltaE= epsAA*DeltaNA + epsAB*DeltaNB; % binding energy
      else
        DeltaE=-epsAA*DeltaNA - epsAB*DeltaNB; % binding energy
      end

      % 3. Calculate rate
      if (DeltaE<0)
        M=rAmax;
      else
        M=rAmax*exp(-DeltaE);
      end

      % accept/reject process
      if(  rand()<M/rAmax  ) % exchange particles
        N_A_hops_processed=N_A_hops_processed+1;
        cnfA(x,y)=Aid2;
        cnfA(xnn,ynn)=Aid;
      end

    %===========================================================
    % ATTEMPT TO HOP PARTICLE B TO NEAREST NEIGHBOUR
    else % Attempt a hop of particle B
      % Get particle types at target site
      Bid2=cnfB(xnn,ynn);
      if (Bid==-1&&Bid2==-1)||(Bid>-1&&Bid2>-1);   % neighbour already occupied -> move not allowed
        continue; % reject -> end iteration
      end
      Aid2= cnfA(xnn,ynn); 
      
      % calculate rate M

      if Bid>-1
        if (Aid==1)      % Particle B is initially inside aggresome
          if Aid2==1     % Particle B remains inside aggresome
            M=nuBin(x,y);
          else           % Particle B leaves aggresome
            M=nuBout(x,y)*exp(epsAB);
          end
        else             % Particle B is initially  outside aggresome
          M=nuBout(x,y);
        end
      else
        if (Aid2==1)      % Particle B is initially inside aggresome
          if Aid==1     % Particle B remains inside aggresome
            M=nuBin(xnn,ynn);
          else           % Particle B leaves aggresome
            M=nuBout(xnn,ynn)*exp(epsAB);
          end
        else             % Particle B is initially  outside aggresome
          M=nuBout(xnn,ynn);
        end
      end
      

      % accept/reject process
      if(  rand()<M/rBmax  ) % exchange particles
        N_B_hops_processed=N_B_hops_processed+1;
        cnfB(x,y)=Bid2;
        cnfB(xnn,ynn)=Bid;
      end
    end % Processing of B particle done
  end % TIME ITERATIONS FINISHED


toc
  %=================================================
  fprintf('==\nTime stepper completed\n');
  fprintf('Time: %e\n', tt);
  fprintf('Number of time steps: %d\n', Ntime);
  fprintf('Number of A hops: %d (%.1f%% of time steps)\n', N_A_hops_processed, N_A_hops_processed*100.0/Ntime);
  fprintf('Number of B hops: %d (%.1f%% of time steps)\n', N_B_hops_processed, N_B_hops_processed*100.0/Ntime);
  fprintf('Accepted processes:  %.1f%%\n', (N_A_hops_processed+ N_B_hops_processed)*100.0/Ntime);

  cmean=  sum( cnfB(:)~=-1 )/(NX*NY);

  n_A=sum(cnfA(:)==1);     % Number of A proteins
  n_B=sum(cnfB(:)>-1); % Number of B proteins
  n_B_bleached=sum(cnfB(:)==0); % Number of B proteins
  N_Bin = sum((cnfB(:)>-1).*(cnfA(:)==1) ); % Number of A proteins that share the same site with a B protein.
  cin  =N_Bin./n_A;
  cout =(n_B-N_Bin)/(NX*NY-n_A);
 % cout = ( cmean*NX*NY-cin*sum(cnfA(:)) )/(NX*NY-sum(cnfA(:)));
  fprintf('Theoretical cout/cin=%e\n', exp(epsAB));
  fprintf('Simulated   cout/cin=%e\n', cout/cin);
  
  
  % Display final configuration  
  subplot(3,2,5) % colourmap
  imagesc(linspace(0,3,NY),linspace(0,1,NX), 1+(cnfA) ); hold on
  daspect([1 1 1])
  title('Final')
  subplot(3,2,6)
  imagesc(linspace(0,3,NY),linspace(0,1,NX), (2*(cnfB==1)+(cnfB==0))); hold on
  daspect([1 1 1])
  title('Final')






  % Normalise FRAP using mean number of proteins B in focus area
  frap_normalised=frap/(Nin/k);
  ifp=fopen(sprintf('%s/frap.out', outdir), 'w');
  for i=1:k
    fprintf(ifp, '%12e %12e\n', trow(i)-tbleach, frap_normalised(i));
  end
  fclose(ifp);

  figure
  subplot(1,2,1)
  loglog(trow(frap>0)-tbleach,  frap(frap>0)/(Nin/k)); hold on
  loglog(trow(frap>0)-tbleach, 3.5*[trow(frap>0)-tbleach].^0.25/50, 'LineWidth', 3)
  loglog(trow(frap>0)-tbleach, 8*1.5*[trow(frap>0)-tbleach].^0.5/100, 'LineWidth', 3)
  xlabel('time [s]')
  ylabel('FRAP')

  subplot(1,2,2)
  plot(trow(frap>0)-tbleach,  frap(frap>0)/(Nin/k)); hold on
  plot(trow(frap>0)-tbleach, 3.5*[trow(frap>0)-tbleach].^0.25/50, 'LineWidth', 3)
  plot(trow(frap>0)-tbleach, 8*1.5*[trow(frap>0)-tbleach].^0.5/100, 'LineWidth', 3)
  xlabel('time [s]')
  ylabel('FRAP')

end

% Input:  Diffusivities (can be 1D or 2D array)
%         Nbins     (integer value): number of bins
% Output: logDspace (1D array):      log10 of the diffusivity
%         counts    (1D array):      Number of diffusivities
function [logDspace, counts]=binMobilities(Diffusivities, Nbins)
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
  end; end;
end








