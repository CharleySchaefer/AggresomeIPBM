function analyseFRAP(dname, label)
  # label: 'A' or 'B'
  # FRAP CURVE WITH TIME DATA IN MONTE CARLO UNITS

  % Get bleach profile
  bleach_profile=importdata(sprintf("%s/BleachProfile.out", dname));

  % Get bleach time
  fsettings=(sprintf("%s/settings.out", dname));
fid = fopen(fsettings);
tline = fgetl(fid);
while ischar(tline)
    tline = fgetl(fid);
    if( strncmp(tline, '  tbleach', 9) ==1)
      tbleach=str2num(tline(10:end));
    end
end
fclose(fid);
  Nheader=1;

  % Get time data
  timedata=importdata(sprintf("%s/timeprogress.out", dname), ' ', Nheader);
  trow=timedata.data(:,2);


  % SWEEP CONFIGURATION FILES
  time_indices=find(trow>=tbleach);      % times after t_bleach
  bleach_indices=find(bleach_profile>0); % only analyse bleach area
  for i=1:length(time_indices)
    cnf=importdata(sprintf("%s/cnf/cnf%c%05d.out", dname, label, time_indices(i)));
    n1(i)=sum(cnf(bleach_indices)==1);
    n2(i)=sum(cnf(bleach_indices)==2);
  end
  ntot=n1+n2;
  ntotmean=mean(ntot);

  % EXPORT TO FILE
  fp=fopen(sprintf('%s/frap.out', dname), 'w');
  fprintf(fp, '#mean number of proteins in bleach area: %d\n', ntotmean);
  fprintf(fp, '%12s %12s %12s\n', '#kmctime', 'unbleached', 'bleached');
  for i=1:length(time_indices)
    fprintf(fp, '%12e %12e %12e\n', trow(time_indices(i))-tbleach, n1(i)/ntotmean, n2(i)/ntotmean);
  end
  fclose(fp);

end 
