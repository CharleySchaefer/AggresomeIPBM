function [cnfA, cnfB] =initialise_cnf_uniform(NX,NY,nA, nB)
    cnfA=-1*ones(NX,NY);    %-1: no A present
    if nA>NX*NY
      fprintf('Error: nA should be smaller than NX*NY.\n');
    elseif nA==NX*NY
      cnfA=ones(NX,NY);     %1:  A present everywhere
    else
      i=0;
      while i<nA
        j=1+floor(rand()*NX*NY);
        if(cnfA(j)==-1)
          cnfA(j)=1; i=i+1; %1:  A present at site j;
        end
      end
    end
    cnfB=-1*ones(NX,NY);    %-1: no B present
    if nB>NX*NY
      fprintf('Error: nA should be smaller than NX*NY.\n');
    elseif nB==NX*NY
      cnfB=zeros(NX,NY);     %0:  unbleached B present everywhere
    else
      i=0;
      while i<nB
        j=1+floor(rand()*NX*NY);
        if(cnfB(j)==-1)
          cnfB(j)=0; i=i+1;  %0:  unbleached B present at site j
        end
      end
    end
end
