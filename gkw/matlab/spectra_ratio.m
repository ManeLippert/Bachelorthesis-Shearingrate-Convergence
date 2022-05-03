% This file defines a matlab function which plots the RATIO of time averaged
% spectra, useful for nonlinear runs.  This script assumes the naming
% conventions of the script gkwnlin.
% 
% spectra_ratio(proj,runname,flux,flux2,species,species2,first,last)
%
% WHERE:
%
% proj and runname are strings for gkwpath.m
%
% species (integer) corresponds to the species no. for which you want the flux,
% in the order of the input file
%
% first and last select the timesstep range for averaging
%
% flux selects the flux you want (same names as spectra diagnostic). 
% Species selects the species you want
%
% normalize (integer) = 1 will normalize the output to total flux.
%
% EXAMPLE: spectra_ratio('proj','runname','pflux','eflux',2,2,100,200)
%
% The spectrum in file is plotted against  runname.krho values in /spectrum/kyspec


function spectra_ratio(proj,filename,flux,flux2,species,species2,first,last)

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
    disp('You must choose flux')
    return;
end
if ~exist('flux2')
	disp('You must choose flux2')
    return;
end

input=read_gkwinput(filename,proj);
nmod=input.GRIDSIZE.nmod

%Select columns for chosen species
start=1+(nmod*(species-1));
finish=nmod*species-1;

start2=1+(nmod*(species2-1));
finish2=nmod*species2-1;

scale=load([gkwpath('kyspec',proj) filename '.krho']);
%scale=2*scale;
disp(['Loaded ' gkwpath('kyspec',proj) filename '.krho'])
flux_file=importdata([gkwpath('spectrum',proj) flux '/' filename]);
flux_file2=importdata([gkwpath('spectrum',proj) flux2 '/' filename]);

a=mean(flux_file(first:last,start:finish));
err_a=std(flux_file(first:last,start:finish))/sqrt(last-first);

b=mean(flux_file2(first:last,start2:finish2));
err_b=std(flux_file2(first:last,start2:finish2))/sqrt(last-first);

av=a./b;
%av=av./scale(2:nmod,1)';

error=(a./b).*sqrt((err_a./a).^2+(err_b./b).^2);

label=[proj ' ' filename ' ' flux ' ' num2str(species) ' / ' flux2 ' ' num2str(species2)];
%errorbar(scale(2:nmod,1),av,error,'+-','DisplayName',label);
plot(scale(2:nmod,1),av,'+-','DisplayName',label);
% Create xlabel
xlabel('ky.rho');
  
    
end