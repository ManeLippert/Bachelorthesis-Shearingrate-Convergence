% This file defines a matlab function which plots the time averaged fluxes
% kx spectrum, useful for nonlinear runs.  This script assumes the naming
% conventions of the script gkwnlin, requires gkw_path.m to be setup.
% 
% fluxes_xspec(proj,runname,flux,species,n_av,normalise)
%
% WHERE:
%
% proj and runname are strings for gkwpath.m
%
% species (integer) corresponds to the species no. for which you want the flux,
% in the order of the input file
%
% n_av (integer) is the number of points to average over working backwards from
% the last point in the file. It is useful to check the saturation by eye
% in fluxes first.
%
% flux (integer) selects the flux you want. The script assumes there are exactly 3 fluxes for each species:
%   1) particle flux
%   2) heat flux
%   3) parallel momenutum flux
%
% normalize (integer) = 1 will normalize the output to total flux.
%
% EXAMPLE: fluxes_xspec('proj','runname',2,2,300,1)
%
% The spectrum in file is plotted against  runname.krho values in /spectrum/kyspec


function []=fluxes_xspec(proj,filename,flux,species,n_av,normalize)
%function []=fluxes(proj,flux,variable,n_av)


if ~exist('proj')
	proj='default';
    disp('You must provide the project name')
    return;
    %No default value
end
if ~exist('filename')
    disp('You must provide the run name')
    filename='default';
    return;
end
if ~exist('flux')
	flux='eflux_spec';
    disp('Plotting energy flux spectra')
end
if ~exist('species')
	species=1;
    disp('Plotting first species')
end
if ~exist('n_av')
    n_av=100;
    disp('Average over last 100 data points')
end
if ~exist('normalize')
    normalize=1;
    disp('Will normalize data to total flux')
end

%Related column of fluxes.dat containing totals
%This assumes three fluxes for each species
flux_no=flux+3*(species-1);
n_av=n_av-1;

switch(flux)
    case(2)
    flux='eflux_xspec' 
    case(3)
    flux='vflux_xspec'
    case(1)
    flux='pflux_xspec'
    otherwise
    disp('Invalid flux chosen')
    return;
end

input=read_gkwinput(filename,proj);
nmod=input.GRIDSIZE.nx

%scale=load([gkwpath('kxspec',proj) filename '.kxrh']);
scale=load([gkwpath('vflux_xspec',proj) 'kxrh']);
disp(['Loaded ' gkwpath('kxspec',proj) filename '.kxrh'])
flux_file=load([gkwpath(flux,proj) filename]);
disp(['Loaded ' gkwpath(flux,proj) filename])

av=mean(flux_file(end-n_av:end,:));
error=std(flux_file(end-n_av:end,:))/sqrt(n_av);
total=sum(av)

fluxes_tot=load([gkwpath('fluxes',proj) filename]);
disp(['Loaded ' gkwpath('fluxes',proj) filename]);
%Should be equal to total
total2=mean(fluxes_tot(end-n_av:end,flux_no))

%Always normalize by number of modes
%In this case total area under curve = total flux
av=av*nmod/scale(1,end);
error=error*nmod/scale(1,end);

%Optionally normalize to total flux
%In this case total area under curve = 1
if(normalize==1)
    av=av/total2;
    error=error/total2;
end

%Select columns for chosen species
start=1+(nmod*(species-1))
finish=nmod*species

label=[proj ' ' filename ' species ', sprintf('%i',species)]
errorbar(scale(1,1:nmod),av(start:finish),error(start:finish),'+-','DisplayName',label)
% Create xlabel
xlabel('kx.rho');
legend('hide')
legend('show')

% Create ylabel
if(normalize==1)
   ylabel(['Frational spectral ' flux]);
   title(['Frational spectral ' flux]);
else
   ylabel(['Spectral ' flux]); 
   title(['Spectral ' flux]);
end    
    
end