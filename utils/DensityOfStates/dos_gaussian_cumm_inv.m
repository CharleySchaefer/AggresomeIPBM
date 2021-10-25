% p: value between 0 and 1
% E: energy
function E=dos_gaussian_cumm_inv( p, Emean, sigma )
  E=Emean-sqrt(2)*sigma*erfcinv(2*p);
end
