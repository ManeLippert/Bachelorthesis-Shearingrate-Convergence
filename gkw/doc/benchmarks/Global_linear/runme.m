

krgene = linspace(1,6,6)/10; 
gam = [0.082 0.195 0.265 0.2625 0.19 0.08]; 
gamg = [0.082 0.185 0.256 0.258 0.213 0.14];

figure(1)
hold off 
plot(krgene,gam,'rx-'); 

om = [0.18 0.36 0.59 0.83 1.04 1.19]; 
omg = [0.18 0.4 0.605 0.83 1.04 1.13];

figure(2);
hold off 
plot(krgene,om,'rx-'); 

load gyg5 
figure(1)
hold on
plot(krgene,gamg,'kx-');
plot(krgene,sqrt(2.)*gyg5(:,2),'bo-')
axis([0 0.7 0 0.3])
box on 

figure(2)
hold on 
plot(krgene,omg,'kx-');
plot(krgene,sqrt(2.)*gyg5(:,3),'bo-')
axis([0 0.7 0 1.3])
box on 


