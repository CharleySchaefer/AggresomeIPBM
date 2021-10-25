function cnf2png_frap()
  map=[0 0 0; 0.1 0.1 0.2; 0.1 0.9 0.0;  0.9 0.1 0.1];
  pth='cnf';

  contents=dir(pth);
  Nfiles=length(contents);

  Ncnf=(Nfiles-2)/2; % Subtract '.' and '..'
  for i=1:Nfiles
    fcnfA=sprintf('%s/cnfA%05d.out', pth, i);
    fcnfB=sprintf('%s/cnfB%05d.out', pth, i);

   
    cnfA=importdata(fcnfA, ' ');
    cnfB=importdata(fcnfB, ' ');

    cnf=( (cnfA==0).*(cnfB~=0) ) + 2*(cnfA==1) + 3*(cnfA==2);
    
    figname=sprintf('img/cnfall%05d.png', i);
       fprintf('creating %s\n', figname);
       imwrite(1+cnf, map, figname);
    

  end
end
