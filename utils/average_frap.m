function average_frap(directory)

  #directory='.';
  NSEEDS=0;
  contents=dir(directory);
  Nfiles=length(contents);
  for i=1:Nfiles
    fname=contents(i).name;
    if(strncmp(fname, 'seed', 4)) #&& strcmp(fname(end-3:end), '.out'))
      NSEEDS=NSEEDS+1; % Count
    end
  end

% MEAN
trow=0;
frow=0;
for i=1:NSEEDS

  data=importdata(sprintf('%s/seed%d/frap.out', directory, i), ' ', 2);
   tmprow1=data.data(:,1);
   tmprow2=data.data(:,2);
  if i>1&&length(tmprow1)>length(trow)
    tmprow1=tmprow1(1:length(trow));
    tmprow2=tmprow2(1:length(trow));
  elseif i>1&& length(tmprow1)<length(trow)
    trow=trow(1:length(tmprow1));
    frow=frow(1:length(tmprow1));
  end
  trow=trow+tmprow1;
  frow=frow+tmprow2;
end
trowmean=trow/NSEEDS;
frowmean=frow/NSEEDS;

% SE
trow=0;
frow=0;
for i=1:NSEEDS

  data=importdata(sprintf('%s/seed%d/frap.out', directory, i), ' ', 2);

   tmprow1=data.data(:,1);
   tmprow2=data.data(:,2);
  if length(tmprow1)>length(trowmean)
    tmprow1=tmprow1(1:length(trowmean));
    tmprow2=tmprow2(1:length(frowmean));
  elseif  length(tmprow1)<length(trowmean)
    trowmean=trowmean(1:length(tmprow1));
    frowmean=frowmean(1:length(tmprow1));
  end


   tmprow1=(tmprow1-trowmean).^2;
   tmprow2=(tmprow2-frowmean).^2;



  trow=trow+tmprow1;
  frow=frow+tmprow2;
end
trowse=sqrt(trow/NSEEDS)/sqrt(NSEEDS-1);
frowse=sqrt(frow/NSEEDS)/sqrt(NSEEDS-1);

fp=fopen(sprintf('%s/averagefrap.out', directory), 'w')
for i=1:length(trow)
  fprintf(fp,'%12e %12e %12e %12e\n', trowmean(i), frowmean(i), trowse(i), frowse(i));
end
fclose(fp);
end
