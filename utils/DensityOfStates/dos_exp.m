function p=dos_exp(E,Emax,width )
  p=zeros(size(E));
  for i=1:length(E)
  %if E(i)<=Emax
    p(i)=exp( (E(i)-Emax)/width  )/width;
  %else
  %  p(i)=0;
  %end
  end
end
