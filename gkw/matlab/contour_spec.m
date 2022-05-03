% This file defines a matlab function which plots the contour plot of spectrum
% against time on a log plot, useful for checking stability of nonlinear runs.
% 
% usage: contour_spec('proj,'file','spectrum',sp)
%
% Where spectrum is one of
%
% kyspec, eflux, pflux, vflux, kyspec_em, eflux_em, pflux_em, vflux_em
% kxspec, eflux_kx, pflux_kx, vflux_kx, kxspec_em, eflux_em_kx, pflux_em_kx, vflux_em_kx
%
% sp is the species number
%
% The spectrum in file is plotted against file.krho values
% Hence this script assumes the naming conventions of the script gkwnlin and gkwpath.m
%
% TO DO: Make the plot without lines between contours (they can be removed manually afterwards).


function []=contour_spec(proj,file,spec,sp)

if ~exist('spec')
   spectrum='kyspec'
else
   spectrum=spec
end

switch spectrum
    case {'kxspec', 'eflux_kx', 'pflux_kx', 'vflux_kx','kxspec_em','eflux_em_kx','pflux_em_kx','vflux_em_kx'}
    scale=['spectrum/kxspec/' file '.kxrh'];
    legend='k_\psi \rho_i';
    sc=load([gkwpath('root',proj) scale]);
    sc=sc';
    
    otherwise
    scale=['spectrum/kyspec/' file '.krho'];
    legend='k_\theta \rho_i';    
    sc=load([gkwpath('root',proj) scale]);   
end

data=load([gkwpath('spectrum',proj) spectrum '/' file]);

time=load([gkwpath('time',proj) file]);

modes=size(sc);
modes=modes(1);

switch spectrum
    case {'kyspec','kxspec'}
        sp=1;
end

start=modes*(sp-1)+1;
fin=modes*sp;

data=data(:,start:fin);

figure
%contourf(time(:,1),sc(:,1),log(abs(data')+1e-20),15);
contourf(time(:,1),sc(:,1),log(abs(data')+1e-20),[-8:0.5:1]*modes/100);


xlabel('t (v_{th}/R)');
ylabel(legend);
title([proj ' ' file])

