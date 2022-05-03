
profiles = read_prof; 

ns = size(profiles.de,1); 

figure; 
for i = 1: ns 

  subplot(ns,5,5*(i-1)+1); 
  plot(profiles.xgr,profiles.de(i,:),'bx-'); 
  xlabel('\psi')
  ylabel(['density species ' int2str(i)]); 
  title('density'); 

end; 

for i = 1: ns 

  subplot(ns,5,5*(i-1)+2); 
  plot(profiles.xgr,profiles.tmp(i,:),'rx-'); 
  xlabel('\psi')
  ylabel(['temperature species ' int2str(i)]); 
  title('Temperature')

end; 

for i = 1: ns 

  subplot(ns,5,5*(i-1)+3); 
  plot(profiles.xgr,profiles.rln(i,:),'bo-'); 
  xlabel('\psi')
  ylabel(['density gradient species ' int2str(i)]); 
  title('R/L_N')

end; 

for i = 1: ns 

  subplot(ns,5,5*(i-1)+4); 
  plot(profiles.xgr,profiles.rlt(i,:),'rx-'); 
  xlabel('\psi')
  ylabel(['temperature gradient species ' int2str(i)]); 
  title('R/L_T');

end; 

for i = 1: ns 

  subplot(ns,5,5*(i-1)+5); 
  plot(profiles.xgr,profiles.uprim(i,:),'mx-'); 
  xlabel('\psi')
  ylabel(['velocity gradient species ' int2str(i)]); 
  title('u^\prime'); 

end; 
