function demo_createBleachProfile2D()

  NX=50;  % Size simulation volume
  NY=150;
  LX=1;   % micron
  LY=3;   % micron

  Xbleach=0.5; % Bleach center of focus (micron)
  Ybleach=0.2; % Bleach center of focus (micron)
  Rbleach=0.4; % Radius of bleach focus
  outfile='bleach_profile.dat'


  BleachProfile=createBleachProfile2D(NX,NY,LX,LY,Xbleach, Ybleach, Rbleach, outfile);

  figure
  imagesc(BleachProfile);
  daspect([1 1 1])
end

