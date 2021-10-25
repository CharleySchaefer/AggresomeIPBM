function cnf2png()
  map=[0 0 0; 0.1 0.8 0.1;;  0.8 0.1 0.1];
  pth='cnf';

  contents=dir(pth);
  Nfiles=length(contents)
  for i=1:Nfiles
    fname=contents(i).name;
    if(strncmp(fname, 'cnf', 3) && strcmp(fname(end-3:end), '.out'))
       prefix=fname(1:end-4);
       matrix=importdata(sprintf('%s/%s', pth, fname), ' ');
       figname=sprintf('img/%s.png', prefix);
       fprintf('creating %s\n', figname);
       imwrite(1+matrix, map, figname);
    end
  end
end
