function p=dos_gaussian_cumm( E, Emean, sigma )
  p=0.5*erfc( (Emean-E)/(sqrt(2)*sigma) );
end
