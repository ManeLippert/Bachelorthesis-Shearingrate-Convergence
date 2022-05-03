clear rhtest_in rhtest_spec

cd(gkwpath('top'))

rhtest_spec=load('./nonspec_nl/other/RH_spec_small_newinit/rhtest');
time=load('./nonspec_nl/time/RH_spec_small_newinit')
fid = fopen('./nonspec_nl/other/RH_nonspec_new3a_hr2a/rhtest', 'r');

frewind(fid);
%lines=size(time,1)-1
lines=1600 % before recurrence problem starts

rhtest_in=zeros(lines,800);

for i = [1:lines] 
  rhtest_in(i,:) = fscanf(fid, '%f',800);
end 

fclose(fid);

figure
set(0,'defaulttextinterpreter','Tex')
set(gca,'box','on','fontsize',12,'Xminortick','on','Yminortick','on')
set(gca,'TickLength',[0.015,0.07])

plot(time(1:lines),rhtest_in(1:lines,101)/rhtest_in(1,101),'DisplayName','Non Spectral')
hold all

%cmplx=complex(rhtest_spec(:,5),rhtest_spec(:,6))

plot(time(1:lines),rhtest_spec(1:lines,6)/rhtest_spec(1,6),'--','DisplayName','Spectral')

xlabel('t (v_{th} / R)')
ylabel('\phi (A.U.)')
title('q=1.3, \epsilon=0.05 GKW RH test benchmark')

nav=8*(66-3); %average 8 periods

% use arbitrary place in radial grid (real part is odd columns)
resid_nonspec=mean(rhtest_in(lines-nav:lines,101))/rhtest_in(1,101);

nav=8*(61);   % average 8 slightly shorter periods

% use 3rd (forward) mode complex part (real part is smaller and out of phase)
resid_spec=mean(rhtest_spec(lines-nav:lines,6))/rhtest_spec(1,6)

analytic=rh_plus(0.05,1.3)

text(30,0.45,['Xiao Catto residual:      ' num2str(analytic,3)],'Fontsize',12)
text(30,0.6,['Spectral residual:         ' num2str(resid_spec,3)],'Fontsize',12)
text(30,0.75,['Non-spectral residual: ' num2str(resid_nonspec,3)],'Fontsize',12)
set(gca,'YGrid','on')
