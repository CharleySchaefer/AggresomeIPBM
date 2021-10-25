function createCellTopology()
  clc; close all;
  % -1: Excluded volume
  %  0: Cytosol
  %  1: Region of attraction (aggresome/stress granule)
  NX=50;
  NY=150;
  fname=sprintf("TopologyTwoAggresomesNX%dNY%d.dat", NX, NY);

  % Initialise 
  cnf=zeros(NX,NY);

  % Two aggresomes
  R=round(NX*0.8/2);
  X1=NX/2; Y1=X1;
  X2=NX/2 ;
  Y2=NY-X1;
  R2target=R*R;
  for i=1:NX;
    RX2=(X2-i)^2 ;
    RX1=(X1-i)^2 ;
    for j=1:NY
      RY2=(Y2-j)^2 ;
      RY1=(Y1-j)^2 ;
      R1=RX1+RY1;
      R2=RX2+RY2;
      if  (  R1<=R2target || R2<=R2target  )
        cnf(i,j)=1;
      end
    end ;
  end

  % One nucleoid
  if 1
  R=round(NX*0.8/2);
  X1=NX/2; Y1=NY/2;
  R2target=R*R;
  for i=1:NX;
    RX1=(X1-i)^2 ;
    for j=1:NY
      RY1=(Y1-j)^2 ;
      R1=RX1+RY1;
      if  (  R1<=R2target   )
        cnf(i,j)=-1;
      end
    end ;
  end
  end

  % Export cnf
  fp=fopen(fname, "w");
  for i=1:NX;
    for j=1:NY
      fprintf(fp, "%3d", cnf(i,j));
    end
    fprintf(fp, "\n"); % end of line
  end
  fclose(fp);

  figure
  imagesc(cnf)


end


