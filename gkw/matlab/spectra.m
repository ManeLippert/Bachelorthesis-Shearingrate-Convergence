% This file defines a matlab function which plots the averaged spectrum
% against time on a log plot, useful for nonlinear runs.
% 
% usage: spectra('proj,'file','spectrum',sp,tstart,tend,optline)
%
% Where spectrum is one of
%
% kyspec, eflux, pflux, vflux, kyspec_em, eflux_em, pflux_em, vflux_em
% kxspec, eflux_kx, pflux_kx, vflux_kx, kxspec_em, eflux_em_kx, pflux_em_kx, vflux_em_kx
%
% sp is the species number
%
% tstart and tend are optional data points to average from (default last 2)
%
% The spectrum in file is plotted against  file.krho values
% Hence this script assumes the naming conventions of the script gkwnlin and gkwpath.m


function []=spectra(proj,file,spec,sp,from,to,line)

if ~exist('spec','var')
   spectrum='kyspec'
else
   spectrum=spec
end

if ~exist('line','var')
   line='-';
end

input=read_gkwinput(file,proj,0);
nx=input.GRIDSIZE.nx;
nmod=input.GRIDSIZE.nmod;

switch spectrum
    case {'kxspec', 'eflux_kx', 'pflux_kx', 'vflux_kx','kxspec_em','eflux_em_kx','pflux_em_kx','vflux_em_kx','kxvort'}
    scale=['spectrum/kxspec/' file '.kxrh'];
    legend='k_\psi \rho_i';
    if exist([gkwpath('root',proj) scale],'file')
     sc=importdata([gkwpath('root',proj) scale])/sqrt(2);
    elseif exist([gkwpath('kxrh',proj) file],'file')
     sc=importdata([gkwpath('kxrh',proj) file])/sqrt(2);
    else
     lx=input.GRIDSIZE.lx;
     kxmin=2*pi/lx
     kxmax=kxmin*nx/2
     sc(1:nx/2)=linspace(0,kxmax,nx/2);
     sc(nx/2+1:nx)=linspace(-kxmax,0,nx/2);
    end
        
    sc=sc';
    semilog=1;
    sign=1
    %sign=-1
    
    otherwise
    scale=['spectrum/kyspec/' file '.krho'];
    if exist([gkwpath('root',proj) scale],'file')
     sc=importdata([gkwpath('root',proj) scale]);
    else
     sc=importdata([gkwpath('krho',proj) file]);
    end
     semilog=0;
     sign=1;
     legend='k_\theta \rho_i';
end

fac=1;
if isfield(input.CONTROL,'spectral_radius')
if (isequal(input.CONTROL.spectral_radius,'true'))
  fac=1   
else
  fac=nx*4
  lx=input.GRIDSIZE.lx;
  %fac=1
end
end

data=importdata([gkwpath('spectrum',proj)  spectrum '/' file]);

if ~exist('to')
   to=size(data,1)
end

if ~exist('from')
   from=to-2
end

%line = '-';
modes=size(sc);
modes=modes(1);

switch spectrum
    case {'kyspec','kxspec','kyvort','kxvort'}
        sp=1;
        %line='--';
end

start=modes*(sp-1)+1;
fin=modes*sp;

if ~exist('from')
   from=100
end

if ~exist('to')
   to = length(data(:,1)) 
end

if (semilog ==1)
  semilogy(sc(:,1),mean(abs(data(from:to,start:fin)))/fac,line,'DisplayName',[proj ' ' file ' ' spectrum ' ' int2str(sp)]);
elseif(semilog==2)  
  semilogx(sc(:,1),mean(abs(data(from:to,start:fin)))/fac,line,'DisplayName',[proj ' ' file ' ' spectrum ' ' int2str(sp)]);
else
  loglog(sign*sc(:,1),mean(abs(data(from:to,start:fin)))/fac,line,'DisplayName',[proj ' ' file ' ' spectrum ' ' int2str(sp)]);
end

xlabel(legend);
ylabel(spectrum);

