% read the outputs for a linear gkw run
%	function [time, anrate, freq, s, phi, apar, n, Tpar, Tper, wflow, thetasq_av, kperp_av,  GKWin, pflux, eflux, vflux]=read_gkwoutput(flnm,proj,full)
% Inputs:
%	flnm:	file name
%	proj:	project name (optional)
%		path for the input files is obtained from the gkwpath function with "proj" in argument
%	full:	full=1 -> load time arrays, full=0 -> load only the last time step
%Note1: the fluxes are still semi-normalised, i.e. the species dependent normalisation (factors n_s, n_s*T_s and m_s*n_s*vth_s for pflux, eflux, vflux, respectively) has been removed but the vthref*rhostar^2 factor is still there. 
%Note2: the electrostatic potential is rotated in the complex plane such as to have real(phi(s=0))=max(real(phi)) and imag(phi(s=0))=0. Same representation than GS2.
% The non rotated potential (i.e. the one given as a GKW output) is given by (phi(:,1) - i*phi(:,2))./phi_norm

function [time, anrate, freq, s, phi, apar, n, Tpar, Tper, wflow, thetasq_av, kperp_av,  GKWin, pflux, eflux, vflux,phi_norm]=read_gkwoutput(flnm,proj,full)

if ~exist('proj') 
	proj=[];
end
if exist('full')~=1 
	full=0;
end
if full==1,
	range='1:end';
else
	range='end';
end


tmp=load([gkwpath('time',proj)  flnm]);
eval(['time=tmp(' range ',1);']);
%Nonlinear runs don't have a growth rate
if size(tmp,2)>1
eval(['anrate=tmp(' range ',2);']);
else
  anrate = time;
end

if size(tmp,2)>2 
  eval(['freq=tmp(' range ',3);'])
else
  freq = NaN.*time;
end

tmp=load([gkwpath('parallel',proj)  flnm]);
s=tmp(:,1);
phi=tmp(:,2:3);
apar=tmp(:,4:5);
n=tmp(:,6:7);
Tpar=tmp(:,8:9);
Tper=tmp(:,10:11);
wflow=tmp(:,12:13);


GKWin=read_gkwinput(flnm,proj);
% normalize phi to be in the same configuration as GS2
n=GKWin.GRIDSIZE.number_of_species;
%L=length(s)./n;
L=GKWin.GRIDSIZE.n_s_grid;
s=s(1:L);
phi=phi(1:L,:);
%I=find(-1e-6<s&s<1e-6);
[tmp I] = max(abs(phi(:,1)+i*phi(:,2)));
phi_0=phi(I,1)+i*phi(I,2);
A=abs(phi_0);
B=phase(phi_0);
phi_tmp=(phi(:,1)+i*phi(:,2))*exp(-i*B)./A;
phi(:,1)=real(phi_tmp);
phi(:,2)=-imag(phi_tmp);
phi_norm=exp(-i*B)./A;

% computes k_perpav using theta=s*pi 
% from Clemente's GS2 routine
%  n=GKWin.GRIDSIZE.number_of_species;
%  ll=length(phi);
%  for ii=1:n,	
%  	II=1+(ii-1)*ll/n:ii*ll/n;
%  	phisq=sum(phi(II,:).^2,2);
%  	theta=s(II).*2*pi;
%  	ith = find(theta > -3.5*pi & theta < 3.5*pi );
%  	zdtheta = diff(theta(ith));
%  	dtheta = [zdtheta(1)./2; zdtheta(1:end-1); zdtheta(end)./2];
%  	at = trapz(phisq(ith).*dtheta);
%  	bt = trapz(theta(ith).^2.*phisq(ith).*dtheta);
%  	thetasq_av(ii) = bt./at;
%  	kperp_av(ii) = sqrt(GKWin.MODE.kthrho.^2.*(1+thetasq_av(ii).*GKWin.GEOM.shat^2) );
%  end
thetasq_av=0;
kperp_av=0;


tmp=load([gkwpath('fluxes',proj)  flnm]);
spcnumber=[1:length([GKWin.SPECIES.z])];
spcnumber=spcnumber([GKWin.SPECIES.z]~=-1|~GKWin.SPCGENERAL.adiabatic_electrons);
pnorm=[GKWin.SPECIES.dens];
enorm=[GKWin.SPECIES.temp].*[GKWin.SPECIES.dens];
vnorm=[GKWin.SPECIES.mass].*[GKWin.SPECIES.dens].*sqrt(2*[GKWin.SPECIES.temp]./[GKWin.SPECIES.mass]);
eval(['pflux=tmp(' range ',1:3:end).*repmat(pnorm(spcnumber),size(tmp(' range ',1),1),1);']);
eval(['eflux=tmp(' range ',2:3:end).*repmat(enorm(spcnumber),size(tmp(' range ',1),1),1);']);
vnorm=[GKWin.SPECIES.mass].*vnorm(spcnumber);
eval(['vflux=tmp(' range ',3:3:end).*repmat(vnorm(spcnumber),size(tmp(' range ',1),1),1);']);

