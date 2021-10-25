function export_cnf(outdir, k, cnfA, cnfB)
  [NX, NY]=size(cnfA);
  ifp=fopen(sprintf('%s/config/cnf%04d.dat', outdir, k), 'w');
      for xn=1:NX
        for ny=1:NY
          if     cnfA(xn,ny)==-1 && cnfB(xn,ny)==-1
            val=0;
          elseif cnfA(xn,ny)==-1 && cnfB(xn,ny)==0
            val=2;
          elseif cnfA(xn,ny)==-1 && cnfB(xn,ny)==1
            val=3;
          elseif cnfA(xn,ny)==1 && cnfB(xn,ny)==-1
            val=1;
          elseif cnfA(xn,ny)==1 && cnfB(xn,ny)==0
            val=4; 
          elseif cnfA(xn,ny)==1 && cnfB(xn,ny)==1
            val=5;
          end
          fprintf(ifp, '%2d', val);
        end
        fprintf(ifp, '\n');
      end
      fclose(ifp);
end
