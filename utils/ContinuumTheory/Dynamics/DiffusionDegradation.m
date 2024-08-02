function DiffusionDegradation()
  clc; close all;
  
  %=============================================
  % USER SETTINGS
  % Physical parameters
  Kdeg    = 1e-6; % Degradation rate
  Kdiff=1.0;
  
  
  DeltaH_row=[1e-3, 2.e-3, 4e-3, 8e-3] ;
  nudiff=0.5;
  DT=2.0;
  B=0.2457;
  
  % Initial distribution
  Nmean0=600;
  Nsigma=200;
  Intensity=1e8;
  
  % Numerical parameters
  Nmax=(Nmean0+5*Nsigma);
  Ndiscr=100;
  Niter=2^17; % 8=2^3; 64=2^6; 512=2^9 4096=2^12
  % END USER SETTINGS
  %=============================================
  Nrow=[1:Ndiscr].*Nmax/Ndiscr;
 % function p=TotConc(x)
 %   p=( log(250)**2*300 - (log(x) - log(250)).**2*300 ).*exp(-x/1000)*1.15;
 % end
  function out=TotConc(Nr)
  
     out= Nr.*exp(-1e-6*(Nr+600).^2) ;
  end
 % fnc(Nrow,2000,150,100,2,1)
  %PC=exp( -(Nrow-Nmean0).^2/(2*Nsigma^2));
  PC=TotConc(Nrow,2000,150,100,2,1);
   % PCtot=sum(Nrow.*PC);
   % PC=PC/PCtot;
  PC=PC-min(PC);
   PCtot=sum(Nrow.*PC);
  PC=PC/PCtot; % 1e8 ~matches this distribution to 
                   % experimental intensity peak (IA + B*IC) 
  PA=zeros(size(PC));
  meanNr= sum(Nrow.*PC)/sum(PC);
  fprintf('Mean RNA length at t=0:  %f\n', meanNr);
 
  figure
  for iter=1:length(DeltaH_row)
    DeltaH=DeltaH_row(iter);
    KdiffA0 = B*Kdiff*Nrow.^(-nudiff);
    KdiffB0 = Kdiff.*Nrow.^(-nudiff).*exp(-Nrow.*DeltaH);
  
      plot(Nrow,      PC, 'k'); hold on
    plot(Nrow,      KdiffB0./(KdiffB0+KdiffA0).*PC, 'b'); hold on
    plot(Nrow,      KdiffA0./(KdiffB0+KdiffA0).*PC, 'r')
    title('Steady State')
    set(gca, 'XScale', 'log')
    legend('tot','Cytoplasm', 'Aggresome')
  end
  
  % return;
  DeltaH=5.5e-3;
  KdiffB0 =   Kdiff.*Nrow.^(-nudiff).*exp(-Nrow.*DeltaH);
  %figure
  %plot(Nrow, KdiffA0./KdiffB0); hold on
  %plot(Nrow, KdiffB0)
  
  
  
  figure
  subplot(2,2,1)
  plot(Nrow, PC, 'b'); hold on
  set(gca, 'XScale', 'log')
  %axis([1 1e4 0 0.08])
  xlabel('Nucleotides')
  ylabel('Cytoplasm')
  subplot(2,2,3)
  plot(Nrow, PA, 'r'); hold on
  set(gca, 'XScale', 'Log')
  %axis([1 1e4 0 0.08])
  xlabel('Nucleotides')
  ylabel('Aggresome')
  
  iter_print=1;
  ifp=fopen(sprintf('SizeVsTime_KdiffByKdeg%.2f.txt', Kdiff/(meanNr*Kdeg)),'w');
  fprintf(ifp, '#DeltaH=%e\n', DeltaH);
  fprintf(ifp, '#Kdeg=%e\n', Kdeg);
  fclose(ifp);
  ifp2=fopen(sprintf('Distributions_KdiffByKdeg%.2f.txt', Kdiff/(meanNr*Kdeg)),'w');
  for iter=1:Niter
  
    if iter>1
      KITER=1;
    elseif iter==128
      KITER=1;
     % figure
    else
      KITER=1;
    end
    for k=1:KITER
    Wgain_deg=zeros(size(PC));
    Wloss_deg=Kdeg*[Nrow-Nrow(1)].*PC;
    Wloss_deg(1)=0;

    for i=1:Ndiscr
      sumf=0;
        for j=i+1:Ndiscr
          sumf=sumf+2*Kdeg*PC(j);
        end
        Wgain_deg(i)=sumf;
        end
      Wgain_deg=Wgain_deg*Nmax/(Ndiscr);
      PC = PC + DT*(Wgain_deg - Wloss_deg);
    end
    
    
    W_diff = -KdiffA0.*PC + KdiffB0.*PA;
    PC = PC + DT*W_diff;
    PA = PA - DT*W_diff;
    if iter==iter_print
      iter_print=ceil(iter_print*1.4);
      fprintf('%12e %12e %12e\n', DT*iter, sum(Nrow.*(PA+PC)), sum(Nrow.*Nrow.*(PA+PC))/sum(Nrow.*(PA+PC)))
      
      ifp=fopen(sprintf('SizeVsTime_KdiffByKdeg%.2f.txt', Kdiff/(meanNr*Kdeg)),'a');
  
      fprintf(ifp,'%e %e\n', DT*iter, sum(Nrow.*(PA+PC))/sum((PA+PC)));
      fclose(ifp);
        
      fprintf(ifp2,'%12e ', DT*iter);
      for i=1:Ndiscr
        fprintf(ifp2,'%12e ',PC(i));
      end
      fprintf(ifp2,'\n');
      fprintf(ifp2,'%12e ', DT*iter);
      for i=1:Ndiscr
        fprintf(ifp2,'%12e ',PA(i));
      end
      fprintf(ifp2,'\n');
      
        
      subplot(2,2,1)
      plot(Nrow, Intensity*PC, 'b'); hold on
      set(gca, 'XScale', 'log')
      subplot(2,2,3)
      plot(Nrow, Intensity*PA, 'r'); hold on
      set(gca, 'XScale', 'Log')
  
      subplot(2,2,[2,4])
      loglog(meanNr*DT*iter, sum(Nrow.*(PA+PC))/sum((PA+PC)), '.k', 'MarkerSize', 15); hold on;
      if sum((PA))>0
        loglog(meanNr*DT*iter, sum(Nrow.*(PA))/sum((PA)), '.r', 'MarkerSize', 15); hold on;
      end
      loglog(meanNr*DT*iter, sum(Nrow.*(PC))/sum((PC)), '.b', 'MarkerSize', 15); hold on;
      xlabel('time * K_{def}*N_{R}^{0}')
      ylabel('< N_R >')
    end
  end
  %fclose(ifp);
 
  
%  plot(Nrow, Nrow.*PC, 'k', 'LineWidth', 3); hold on
%  plot(Nrow, Nrow.*PA, 'g', 'LineWidth', 3); hold on
 % subplot(1,2,2)
 % loglog(Nrow, PA./PC, 'r', 'LineWidth', 3); hold on
  %legend('Cytoplasm', 'Aggresome')

end
