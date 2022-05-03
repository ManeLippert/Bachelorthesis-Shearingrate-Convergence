set(0,'defaulttextinterpreter','none')
figure

subplot(2,3,1)

growth_rates_benchmark('miller_scans','kappa*','GEOM.kappa',2)
hold all
growth_rates_benchmark('miller_scans','kappa*','GEOM.kappa',2,1)
growth_rates_benchmark('miller_scans','base_circ','GEOM.kappa',2)
ylim([0 0.5])
set(gca,'box','on','fontsize',11,'Xminortick','on','Yminortick','on')
set(gca,'TickLength',[0.015,0.07])

subplot(2,3,2)

growth_rates_benchmark('miller_scans','skappa*','GEOM.skappa',2)
hold all
growth_rates_benchmark('miller_scans','skappa*','GEOM.skappa',2,1)
growth_rates_benchmark('miller_scans','base_circ','GEOM.skappa',2)
ylim([0 0.5])
set(gca,'box','on','fontsize',11,'Xminortick','on','Yminortick','on')
set(gca,'TickLength',[0.015,0.07])

subplot(4,3,3)
 
growth_rates_benchmark('miller_scans','betaprime*','SPCGENERAL.betaprime_ref',2)
hold all
growth_rates_benchmark('miller_scans','betaprime*','SPCGENERAL.betaprime_ref',2,1)
growth_rates_benchmark('miller_scans','base_circ','SPCGENERAL.betaprime_ref',2)
growth_rates_benchmark('miller_scans','gradp_betaprime*','GEOM.gradp',2)
ylim([0 0.5])
set(gca,'box','on','fontsize',11,'Xminortick','on','Yminortick','on')
set(gca,'TickLength',[0.015,0.07])

subplot(4,3,6)
 
growth_rates_benchmark('miller_scans','mhd_alpha*','GEOM.gradp',2)
hold all
growth_rates_benchmark('miller_scans','mhd_alpha*','GEOM.gradp',2,1)
growth_rates_benchmark('miller_scans','base_circ','GEOM.gradp',2)
ylim([0 0.5])
set(gca,'box','on','fontsize',11,'Xminortick','on','Yminortick','on')
set(gca,'TickLength',[0.015,0.07])

subplot(2,3,4)

growth_rates_benchmark('miller_scans','delta*','GEOM.delta',2)
hold all
growth_rates_benchmark('miller_scans','delta*','GEOM.delta',2,1)
growth_rates_benchmark('miller_scans','base_circ','GEOM.delta',2)
ylim([0 0.5])
set(gca,'box','on','fontsize',11,'Xminortick','on','Yminortick','on')
set(gca,'TickLength',[0.015,0.07])

subplot(2,3,5)

growth_rates_benchmark('miller_scans','s2delta*','GEOM.sdelta',2)
hold all
growth_rates_benchmark('miller_scans','sdelta*','GEOM.sdelta',2,1)
growth_rates_benchmark('miller_scans','base_circ','GEOM.sdelta',2)
ylim([0 0.5])
set(gca,'box','on','fontsize',11,'Xminortick','on','Yminortick','on')
set(gca,'TickLength',[0.015,0.07])

subplot(2,3,6)

growth_rates_benchmark('miller_scans','*dRmil*','GEOM.drmil',2)
hold all
growth_rates_benchmark('miller_scans','*dRmil*','GEOM.drmil',2,1)
growth_rates_benchmark('miller_scans','base_circ','GEOM.drmil',2)
ylim([0.1 0.6])
set(gca,'box','on','fontsize',11,'Xminortick','on','Yminortick','on')
set(gca,'TickLength',[0.015,0.07])
