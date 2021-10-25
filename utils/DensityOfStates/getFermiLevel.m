% Uses bisection method
function EF=getFermiLevel(Emat, nAtarget, maxiter)
  EFl=min(Emat(:));
  EFu=max(Emat(:));
  iter=0;
  while iter<maxiter
    iter=iter+1;
    EF=0.5*(EFl + EFu);
    nA=sum( Emat(:)<=EF );
    if nA>nAtarget
      EFu=EF;
    elseif nA<nAtarget
      EFl=EF;
    else
     break;
    end
  end
  if iter==maxiter
    fprintf('Warning: getFermiLevel() has not converged.\n');
  end
end
