load cycgene ;
load cycgkw_highres;

gkw = cycgkw_highres;
gene = cycgene; 

figure(99);
hold on; 

plot(gkw(:,1),gkw(:,2),'rs')
plot(gkw(:,1),gene(:,1),'bo')

ns = size(gkw,1);
for i = 1: ns 
  x = [gkw(i,1) gkw(i,1)];
  y = [gkw(i,2)-gkw(i,3) gkw(i,2)+gkw(i,3)];
  plot(x,y,'r-')
end;

x = linspace(6,20,100);
y = 15.4*(1 - 6./x);
axis([0 20 0 10]);
plot(x,y,'k');

figure(98) 
hold on; 
y = y / 2. / sqrt(2.)
plot(x,y)
plot(gkw(:,1),gkw(:,2)/2/sqrt(2.),'ro');
axis([0 16 0 4])
xlabel('R/L_T')
ylabel('\chi_i (\rho_i^2 v_{th} / L_N)')

