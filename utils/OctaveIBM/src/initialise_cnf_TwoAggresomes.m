function [cnfA, cnfB]=initialise_cnf_TwoAggresomes(NX,NY,epsAB)
    % Generated morphology:
    cin=0.7;       % inside concentration
    cout=1.5*cin*exp(epsAB);
  N_agg  =2; % Number of aggresomes
  phi_agg=0.33; % volume fraction occupied by aggresome
  V_agg = phi_agg*NX*NY/N_agg;
  Radius = sqrt(V_agg/pi)

  % DEFINE CONFIGURATION
  X1=NX/2; X2=NX/2  
  Y1=NX/2; Y2=NY-NX/2
  N_agg  =2; % Number of aggresomes
  phi_agg=0.33; % volume fraction occupied by aggresome
  V_agg = phi_agg*NX*NY/N_agg;
  Radius = sqrt(V_agg/pi)
  X1=NX/2; X2=NX/2 ;
  Y1=NX/2; Y2=NY-NX/2;
  cnfA=-1*ones(NX,NY);  cnfB=zeros(NX,NY);
  for i=1:NX
    for j=1:NY
      R1=sqrt((i-X1)^2+(j-Y1)^2);
      R2=sqrt((i-X2)^2+(j-Y2)^2);
      if R1<Radius || R2 < Radius %|| R3 < Radius
        cnfA(i,j)=1;

        if R1<Radius
          if rand()<cin
            cnfB(i,j)= 0;       % Bleached
          else
            cnfB(i,j)=-1;       % Remove 
          end
        elseif  R2<Radius
          if rand()<cin
            cnfB(i,j)= 0;       %unBleached
          else
            cnfB(i,j)=-1;       % Remove
          end
        end
      elseif rand()>cout
          cnfB(i,j)=-1;         % Remove
      end
    end
  end
end
