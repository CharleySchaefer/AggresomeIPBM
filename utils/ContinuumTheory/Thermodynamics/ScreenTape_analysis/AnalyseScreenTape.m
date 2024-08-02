function AnalyseScreenTape()
close all; clc
  % Charley Schaefer, University of York, 2023/04/27
  %
  %   - Import screentape data
  %   - User control data to determine baseline + error margin
  %   - Filter intensity data using the baseline -> data whose error margin overlaps with that of the baseline is excluded.
  %   - use lsqnonlin to fit the model
  % 
  %             IA/IC = Prefac*exp( RNAlength * DeltaH/kT)
  %
  %     to the data.
  %
  %===========================================================
  % USER INPUT
  mprog='octave';                   % octave | matlab
  if 0
    f_in='Figure_2c_screentape.csv';  % FILE NAME
    delim=' ';                        % DELIMITER BETWEEN DATA
    Nheader=4;                        % NUMBER OF HEADER LINES
    col_RNA_length=2;                 % COLUMN NUMBERS
    col_CR_mean=6;
    col_CR_err=7;
    col_control_mean=11;
    col_control_err=12;
    col_AR_mean=16;
    col_AR_err=17;
  param0=[0.4,0.006]
    baseline_mean=-1; % Automatic detection baseline
  %   Summary results:
  %     R2=0.989244 (goodness of fit)
  %     Prefac   =0.315376 +/- 0.083348
  %     DeltaH/kT=0.005486 +/- 0.000699
  %     Note: kT=2.5kJ/mol at room temperature
  %     so, DeltaH=14+/-2 J/mol
  elseif 1
    f_in='20230610_screentape.csv';  % FILE NAME
    delim=' ';                        % DELIMITER BETWEEN DATA
    Nheader=3;                        % NUMBER OF HEADER LINES
    col_RNA_length=1;                 % COLUMN NUMBERS
    col_CR_mean=2;
    col_CR_err=3;
    col_control_mean=4;
    col_control_err=5;
    col_AR_mean=6;
    col_AR_err=7;
  param0=[0.1,0.0006]
 % param0=[0.1,0.006,0.65]
    baseline_mean=70;
    RNA_cutoff=25;
  %   Summary results:
%  R2=0.989378
%Prefac   =0.316275 +/- 0.083062
%DeltaH/kT=0.005481 +/- 0.000694
%R2=0.989378

  end
  exponent=1.0
  %===========================================================
  
  
  clc; close all;

  %===========================================================
  % IMPORT DATA
  data=importdata(f_in, delim,Nheader);
  try data=data.data;
  end
  RNA_length=data(:,col_RNA_length);
  CR_mean     =data(:,col_CR_mean);
  CR_err      =data(:,col_CR_err);
  control_mean=data(:,col_control_mean);
  control_err =data(:,col_control_err);
  AR_mean     =data(:,col_AR_mean);
  AR_err      =data(:,col_AR_err);
  hfig1=figure
  subplot(2,1,1)
  errorbar(RNA_length, CR_mean, CR_err, '.k'); hold on
  errorbar(RNA_length, control_mean, control_err, '.k'); hold on
  errorbar(RNA_length, AR_mean, AR_err, '.k'); hold on
  set(gca,'XScale','log');

  %===========================================================
  % DETERMINE BASELINE USING THE CONTROL
  tolerance=mean(control_err./control_mean);
  Ndata=length(RNA_length)
  cumm=zeros(Ndata,1);
  cumm(1)=control_mean(end);
  for j=2:Ndata
    cumm(j)=cumm(j-1)+control_mean(Ndata+1-j);
  end
  cumm=cumm./[1:Ndata]';
  cumm_err=zeros(Ndata,1);
  err=0; j=1;
  while err<tolerance
    j=j+1;
    cumm_err(j)=sqrt( sum((control_mean(Ndata+1-j:end)-cumm(j)).^2)/j  );
    err=cumm_err(j)/cumm(j);
  end
  j=j-1;
  RNA_length(Ndata+1-j)
  if baseline_mean==-1
  baseline_mean=cumm(j)
  baseline_err=cumm_err(j)
  else
  baseline_err=0; %cumm_err(j)
  end
  
  %===========================================================
  % FILTER DATA
  for j=1:Ndata
    ind=Ndata+1-j;
    if AR_mean(ind)==0 ||CR_mean(ind)==0||control_mean(ind)==0 || CR_mean(ind)-CR_err(ind)<baseline_mean+baseline_err || AR_mean(ind)-AR_err(ind)<baseline_mean+baseline_err || RNA_length(ind) <=RNA_cutoff
    RNA_length(ind)=[]; %data(:,col_RNA_length);
  CR_mean(ind)     =[]; %data(:,col_CR_mean);
  CR_err(ind)      =[]; %data(:,col_CR_err);
  control_mean(ind)=[]; %data(:,col_control_mean);
  control_err(ind) =[]; %data(:,col_control_err);
  AR_mean(ind)     =[]; %data(:,col_AR_mean);
  AR_err(ind)      =[]; %data(:,col_AR_err);
    end
  end
  Iratio_mean=AR_mean./CR_mean;
  Iratio_err=Iratio_mean.*sqrt((AR_err./AR_mean).^2+(CR_err./CR_mean).^2);
  Iratio_mean;
  Iratio_err;
  figure(hfig1)
  subplot(2,1,1)
  errorbar(RNA_length, CR_mean, CR_err); hold on
  errorbar(RNA_length, control_mean, control_err); hold on
  errorbar(RNA_length, AR_mean, AR_err); hold on
  xlabel('RNA length')
  ylabel('Intensity')
  subplot(2,1,2)
  log(RNA_length)
  
  errorbar(RNA_length, Iratio_mean, Iratio_err, '.k'); hold on
  %plot(RNA_length,0.4*exp(RNA_length*0.006));
  set(gca,'XScale','log');
  set(gca,'YScale','log');
  xlabel('RNA length')
  ylabel('Intensity ratio')
  
  
  %===========================================================
  % CURVE FIT
  include_optimisation_pkg(mprog);
  function cost=fit_fnc(param, xdata,ydata, exponent)
    %cost=(log(param(1)*exp(param(2)*xdata))-log(ydata));
    cost=(log(param(1))+(param(2)*xdata.^(exponent))-log(ydata));
  end
  function cost=fit_fnc2(param, xdata,ydata)
    %cost=(log(param(1)*exp(param(2)*xdata))-log(ydata));
    cost=(log(param(1))+(param(2)*xdata.^(param(3)))-log(ydata));
  end
  paramL=1e-2*param0; %paramL(3)=2/3;
  paramU=1e2*param0;  %paramU(3)=3/3;
  options = optimset(...   
            'MaxIter',1000,...
            'Display','off',...
            'MaxFunEvals',100000,...
            'TolX',1e-10,...
            'TolFun',1e-10);
  [param, resnorm, residual, qqq1, qqq2, qqq3, jacobian]=lsqnonlin(@(param)fit_fnc(param, RNA_length,Iratio_mean,exponent), param0, paramL, paramU);
%  [param, resnorm, residual, qqq1, qqq2, qqq3, jacobian]=lsqnonlin(@(param)fit_fnc2(param, RNA_length,Iratio_mean), param0, paramL, paramU);
  chisquare= resnorm; 
  Npar = 2;                              % Number of fit parameters
  dof  = length(Iratio_mean) - Npar;            % Degrees of freedom
p=chi2cdf(chisquare, dof);
  C          = inv(full(jacobian)'*full(jacobian));
  covarpar   = chisquare/dof*C;            % Variance-covariance matrix
  param_std  = sqrt(diag(covarpar))';    % Standard deviation (STD)

  corr_mat = C./sqrt(diag(C)*diag(C)');  % Correlation matrix 

  Rsquared=1-resnorm/((Iratio_mean-mean(Iratio_mean))'*(Iratio_mean-mean(Iratio_mean)));
  %===========================================================
  fprintf('Results\n')
  fprintf('R2=%f\n',Rsquared );
  fprintf('Prefac   =%f +/- %f\n',param(1),param_std(1) );
  fprintf('DeltaH/kT=%f +/- %f\n',param(2),param_std(2) );
  %fprintf('exponent =%f +/- %f\n',param(3),param_std(3) );
  fprintf('R2=%f\n',Rsquared );
  
  fprintf('summary: %f %f %f %f %f %f\n', exponent, param(1), param_std(1), param(2), param_std(2), Rsquared)
  
  ifp=fopen('processed_data.txt', 'w');
    fprintf(ifp, '%12s %12s %12s\n', '#RNA_length', 'Iratio', 'Iratio_err')
  for i=1:length(Iratio_mean)
    fprintf(ifp,'%12e %12e %12e\n', RNA_length(i),Iratio_mean(i), Iratio_err(i));
  end
  fclose(ifp);
  
  plot(RNA_length,param(1)*exp(RNA_length.^(exponent)*param(2)), '-r', 'LineWidth', 2);
 % plot(RNA_length,param(1)*exp(RNA_length.^(param(3))*param(2)), '-r', 'LineWidth', 2);
%  legend('AR/CR data', 'initial guess', 'final fit')
end

function include_optimisation_pkg(mprog)
  % Charley Schaefer, University of York 2021
  switch mprog % RUN IN MATLAB
    case 'matlab'
      % Check licences
      if ~license('test','optimization_toolbox')
        error('ERROR: this code uses MATLAB''s optimization toolbox!')
        return
      end
      if ~license('test','statistics_toolbox')
         error('ERROR: this code uses MATLAB''s statistics toolbox!')
       return
      end
    case 'octave' % RUN IN OCTAVE
     % install octave packages: run "pkg install -forge struct optim statistics"
      %                           in command window
      pkg load struct
      pkg load optim       % lsqnonlin; dependency: struct
      pkg load statistics  % lhsdesign
    
  end
end
