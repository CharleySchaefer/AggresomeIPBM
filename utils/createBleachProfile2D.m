function BleachProfile=createBleachProfile2D(varargin)
  clc; close all;
  [filepath,scriptname,ext] = fileparts(mfilename('fullpath'));
  addpath(sprintf('%s',filepath));

  if nargin ~= 8
    fprintf('usage: createBleachProfile2D(NX,NY,LX,LY,xbleach, ybleach, Rbleach, outfile)\n');
    fprintf('       NX,NY:           resolution (number of pixels).\n');
    fprintf('       LX,LY:           size dimensions.\n');
    fprintf('       xbleach,ybleach: center of laser focus.\n');
    fprintf('       Rbleach:         radius of laser focus.\n');
    fprintf('       See demo_createBleachProfile2D.m\n');
    BleachProfile=0;
  else

  %=========================================================
  % Future development: unsharp boundary; z-dependence
  % 0: unbleached
  % 1: bleached
  NX=varargin{1};  % Size simulation volume
  NY=varargin{2}; 
 % NZ=varargin{3};
  LX=varargin{3};   % micron
  LY=varargin{4};   % micron
 % LZ=varargin{6};   % micron - currently unused!
  Xbleach=varargin{5}; % Bleach center of focus (micron)
  Ybleach=varargin{6}; % Bleach center of focus (micron)
  Rbleach=varargin{7}; % Radius of bleach focus
  outfile=varargin{8}; % Radius of bleach focus
  %=========================================================

  %========================================================= 
  % Check for errors
  dy=LY/NY;
  dx=LX/NX;
  dz=dx; % dz currently unused
  if(dx~=dy)
    error(sprintf('grid spacing dx=%f and dy=%f should be equal\n', dx,dy));
  end
  %if( NZ~=1 )
  %  error('NZ=1 not implemented yet.\n');
  %end
  %=========================================================

  %========================================================= 
  % CORE 
  BleachProfile=zeros(NX,NY);  
  Rbleachsq=Rbleach*Rbleach;
  for i=1:NX
    xsq= (i*dx-Xbleach)^2;    % x distance from laser focus squared
    for j=1:NY
      ysq= (j*dy-Ybleach)^2;  % y distance from laser focus squared
      Rsq=xsq+ysq; % Distance from laser focus squared
      if(Rsq<Rbleachsq)
        BleachProfile(i,j)=1; % Bleach site!
      end
    end
  end
  end
  %=========================================================


  %=========================================================
  % EXPORT
  fp=fopen(outfile, 'w');
  for i=1:NX
    for j=1:NY
      fprintf(fp, '%2d', BleachProfile(i,j) );
    end
    fprintf(fp, '\n');
  end
  fclose(fp);
end

