% to check GKW outputs: plots of growth rate, fluxes, potential and distribution function
%	function check_gkwoutput(flnm,proj,col,fdis)
% Inputs:
%	flnm		file name
%	proj:	project name (optional)
%		path for the input files is obtained from the gkwpath function with "proj" in argument
%	col:		color used for the 2D plots
%	fdis:	fdis=1 -> plots the distribution function

function []=check_gkwoutput(flnm,proj,col,fdis)


% default 
if ~exist('proj')
	proj=[];
end
if ~exist('col')
	col='b';
end
if ~exist('fdis')
	fdis=0;
end
flpth=gkwpath('input',proj);


% read standard outputs:
[time, anrate, freq, s, phi, apar, n, Tpar, Tper, wflow, ...
thetasq_av, kperp_av,  GKWin, pflux, eflux, vflux,phi_norm]=read_gkwoutput(flnm,proj,1);


% growth rate convergence
hf=findobj('Name','growth_rate');
if isempty(hf)
	hf=figure;
	set(gcf,'Name','growth_rate');
	hold on
	xlabel('time')
	ylabel('\gamma')
	title('growth rate convergence')
	set(gcf,'UserData',[0 max(time) 0.8*anrate(end) 1.2*anrate(end)])
end
figure(hf)
plot(time,anrate,col)
lim=get(gcf,'UserData');
lim(2)=max([lim(2) max(time)]);
lim(3)=min([lim(3) 0.8*anrate(end)]);
lim(4)=max([lim(4) 1.2*anrate(end)]);
axis(lim)
set(gcf,'UserData',lim)

% fluxes convergence
hf=findobj('Name','fluxes');
flux=pflux(:,1)./eflux(:,1);
%flux=1;
if isempty(hf)
	hf=figure;
	set(gcf,'Name','fluxes');
	hold on
	xlabel('time')
	ylabel('\Gamma_1/Q_1')
	title('fluxes convergence')
	set(gcf,'UserData',[0 max(time) 0.8*flux(end) 1.2*flux(end)])
end
figure(hf)
plot(time,flux,col)
lim=get(gcf,'UserData');
lim(2)=max([lim(2) max(time)]);
lim(3)=min([lim(3) min([0.8*flux(end) 1.2*flux(end)]) -1e-13]);
lim(4)=max([lim(4) max([0.8*flux(end) 1.2*flux(end)]) 1e-13]);
axis(lim)
set(gcf,'UserData',lim)


% potential, real part
hf=findobj('Name','real_imag');
if isempty(hf)
	hf=figure;
	set(gcf,'Name','real_imag');
	hold on
	xlabel('s')
	ylabel('\phi_R')
	title('potential')
	set(gcf,'UserData',[min(s) max(s) 0.95*min(phi(:,1)) 1.05*max(phi(:,1))])
end
figure(hf)
plot(s,phi(:,1),col,'+-')
lim=get(gcf,'UserData')
lim(3)=min([lim(3) min([0.95*min(phi(:,1)) 1.05*min(phi(:,1))])]);
lim(4)=max([lim(4) max([0.95*min(phi(:,1)) 1.05*max(phi(:,1))])]);
axis(lim)
set(gcf,'UserData',lim)

% potential, norm
hf=findobj('Name','phi_norm');
if isempty(hf)
	hf=figure;
	set(gcf,'Name','phi_norm');
	hold on
	xlabel('s')
	ylabel('\phi_norm')
	title('potential')
	set(gcf,'UserData',[min(s) max(s) 0.95*min(phi(:,2)) 1.05*max(phi(:,2))])
end
figure(hf)
plot(s,sqrt(phi(:,1).^2+phi(:,2).^2),col,'+-')
lim=get(gcf,'UserData');
lim(3)=min([lim(3) 0.95*min(phi(:,2))]);
lim(4)=max([lim(4) 1.05*max(phi(:,2))]);
axis(lim)
set(gcf,'UserData',lim)


% distribution function
if fdis==1
  vpargr=load([gkwpath('distr1',proj)  flnm]); % parallel velocity grid: vpar/vth (vth is species dependent)
  vperpgr=load([gkwpath('distr2',proj)  flnm]); % perpendicular velocity grid: vperp/vth 
  f_im=load([gkwpath('distr3',proj)  flnm]); % distribution function at maximum potential for the last species
  f_r=load([gkwpath('distr4',proj)  flnm]);
  f=sqrt(f_im.^2+f_r.^2);
  ph=angle((f_r+i*f_im)./phi_norm); % phase(f*phi)@s=0
  G=read_geom(flnm,proj);
  [phim, ii] = max(sqrt(phi(:,1).^2)+phi(:,2).^2);
  bratio=sqrt((max(G.bn)-G.bn(ii))./G.bn(ii));
  
  % amplitude
  figure
  surf(vpargr,vperpgr,f)
  hold on
  shading flat
  colorbar
  xlabel('vpar')
  ylabel('vperp')
  zlabel('abs(\deltaf)')

  vperpmax=max(max(vperpgr));
  h1=fill3([0 vperpmax*bratio vperpmax*bratio 0],[0 vperpmax vperpmax 0],[0 0 max(max(f)) max(max(f))],'b');
  set(h1,'FaceAlpha',0.5)
  h2=fill3(-[0 vperpmax*bratio vperpmax*bratio 0],[0 vperpmax vperpmax 0],[0 0 max(max(f)) max(max(f))],'b');
  set(h2,'FaceAlpha',0.5)

  % phase
  figure
  surf(vpargr,vperpgr,cos(ph))
  hold on
  shading flat
  colorbar
  xlabel('vpar')
  ylabel('vperp')  
  zlabel('phase(\deltaf\phi)')

end

