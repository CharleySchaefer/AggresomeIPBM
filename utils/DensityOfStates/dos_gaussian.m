function p=dos_gaussian( E, Emean, sigma )
  two_sigma_square=2*sigma*sigma;
  p=exp(-(Emean-E).^2/(two_sigma_square))/sqrt( pi*two_sigma_square  );
end
