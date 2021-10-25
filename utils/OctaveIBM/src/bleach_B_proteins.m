%================================================================
% PHOTOBLEACHING
% mask:  pixels corresponding to bleaching area
% Nmask: number of pixels in bleaching area
% cnfB:  NXxNY array with entry values cnf(x,y)
%         -1: no B present at site x,y; 
%          0: unbleached B at site x,y; 
%          1:   bleached B at site x,y; 
% This function 'bleaches' cnfB at focus point X1,X1, with given radius. 
function [cnfB, mask, Nmask]=bleach_B_proteins(cnfB, X1, Y1, radius)
  [NX,NY]=size(cnfB);
  % BLEACH: Label proteins B in a circular laser focus
  mask=[];
  for i=1:NX
    for j=1:NY
      R1=sqrt((i-X1)^2+(j-Y1)^2);
        if R1<radius
          if( cnfB(i,j)== 0);       % unbleached present
            cnfB(i,j)=1;             % bleach
          end
          mask=[mask; i,j];
        end
    end
  end
  [Nmask, b]=size(mask); % Number of bleached pixels
end
