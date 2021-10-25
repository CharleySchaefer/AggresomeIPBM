function average_lengthscale(directory)
  #  TIME AND LENGTH DATA IN [s] and [micron], respectively

# TODO: script assumed the directory names are seed1 seed2 seed3, etc starting always from seed1
#       more flexibility is to be implemented.

NSEEDS=0;
  contents=dir(directory);
  Nfiles=length(contents);
  for i=1:Nfiles
    fname=contents(i).name;
    if(strncmp(fname, 'seed', 4)) #&& strcmp(fname(end-3:end), '.out'))
      NSEEDS=NSEEDS+1; % Count
    end
  end


 % INPUT; TODO: USE AS INPUT ARGUMENT TO SCRIPT
% GET AVERAGE TIME AND LENGTHSCALE DATA
trowmean=0;
frowmean=0;
for i=1:NSEEDS
  data=importdata(sprintf('%s/seed%d/AnalysedLengthScales.out', directory, i), ' ', 2);
  trowmean=trowmean+data.data(:,1); % Time 
  frowmean=frowmean+data.data(:,4); % Length scale
end
trowmean=trowmean/NSEEDS;
frowmean=frowmean/NSEEDS;

% GET STANDARD ERROR
trow=0;
yrow=0;
for i=1:NSEEDS
  data=importdata(sprintf('%s/seed%d/AnalysedLengthScales.out', directory, i), ' ', 2);
  trow=trow+(data.data(:,1)-trowmean).^2;
  yrow=yrow+(data.data(:,4)-yrow).^2;
end
trowse=sqrt(trow/NSEEDS)/sqrt(NSEEDS-1);
yrowse=sqrt(yrow/NSEEDS)/sqrt(NSEEDS-1);


% EXPORT 
fp=fopen(sprintf('%s/AverageLengthScale.out', directory), 'w');
  fprintf(fp, '%12s %12s %12s %12s\n', '#time_mean', 'R_mean', 'time_SE', 'R_SE');
for i=1:length(trow)
  fprintf(fp, '%12e %12e %12e %12e\n', trowmean(i), frowmean(i), trowse(i), yrowse(i));
end
fclose(fp);
end
