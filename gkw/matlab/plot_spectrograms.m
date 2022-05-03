% This function plots GKW contour spectrograms of spectra against time
% This is useful for understanding the evolution of a nonlinear run
% and identifying high-k pileup problems, or low-k blowup.
%
% usage: plot_spectrograms('proj,'file',tstart,tend)
%
% tstart and tend are optional integers for start and end range of plot
%
% Requires gkwpath.m and gkwnlin folder structure.
% 


function []=plot_spectrograms(proj,file,tstart,tend)

input=read_gkwinput(file,proj,0);
nx=input.GRIDSIZE.nx;
nmod=input.GRIDSIZE.nmod;
time=load([gkwpath('time',proj) file]);

if ~exist('tstart','var'); tstart=1
end
if ~exist('tend','var'); tend=size(time,1)
end

legendx='k_\psi \rho_i';
legendy='k_\theta \rho_i';

scale=['spectrum/kxspec/' file '.kxrh'];
if exist([gkwpath('root',proj) scale],'file')
    scx=importdata([gkwpath('root',proj) scale]);
elseif exist([gkwpath('kxrh',proj) file],'file')
    scx=importdata([gkwpath('kxrh',proj) file]);
else
    lx=input.GRIDSIZE.lx;
    kxmin=2*pi/lx
    kxmax=kxmin*nx/2
    scx(1:nx/2)=linspace(0,kxmax,nx/2);
    scx(nx/2+1:nx)=linspace(-kxmax,0,nx/2);
end
scx=scx';
scx=scx(:,1);

scale=['spectrum/kyspec/' file '.krho'];
if exist([gkwpath('root',proj) scale],'file')
    scy=importdata([gkwpath('root',proj) scale]);
else
    scy=importdata([gkwpath('krho',proj) file]);
end
scy=scy(:,1);

kyspec=importdata([gkwpath('spectrum',proj)  '/kyspec/' file]);
kyvort=importdata([gkwpath('spectrum',proj)  '/kyvort/' file]);
kyspec_em=importdata([gkwpath('spectrum',proj)  '/kyspec_em/' file]);
kxspec=importdata([gkwpath('spectrum',proj)  '/kxspec/' file]);

figure
a1=subplot(2,2,1);

contourf(time(tstart:tend,1),scy,log(kyspec(tstart:tend,:)'),20);
xlabel('t_N');
ylabel(legendy)
title('log(kyspec)')
shading flat

a2=subplot(2,2,2);
contourf(time(tstart:tend,1),scy,log(kyvort(tstart:tend,:)'),20);
xlabel('t_N');
ylabel(legendy)
title('log(kyvort))')
shading flat

a3=subplot(2,2,3);
contourf(time(tstart:tend,1),scy,log(kyspec_em(tstart:tend,:)'),20);
xlabel('t_N');
ylabel(legendy)
title('log(kyspec EM)')
shading flat

a4=subplot(2,2,4);

%Circular shift for FFT of nonspectral
%scx=circshift(scx,[size(scx,1)/2]); 
%kxspec=circshift(kxspec,[0 size(kxspec,2)/2]);
[scx ind]=sort(scx);
kxspec=kxspec(:,ind);

contourf(time(tstart:tend,1),scx,log(abs(kxspec(tstart:tend,:))'),20);
xlabel('t_N');
ylabel(legendx);
title('log(kxvort)')
shading flat

linkaxes([a1 a2 a3 a4],'x');
%linkaxes([a1 a2 a3],'xy');

