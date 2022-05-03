%  rh_plus(eps,q)
%
%  Returns Xiao-Catto result for Rosenbluth Hinton GAM residual given q and eps
%  use with (e.g.) fplot(@(x) rh_plus(x,1.5),[0 1])  
%
function [residual] = rh_plus(eps,q)
 
  Odum= 1.6*eps^(3/2)+0.5*eps*eps+0.36*eps^(5/2);
  
  residual= 1/(1+q^2*Odum/(eps*eps));

end