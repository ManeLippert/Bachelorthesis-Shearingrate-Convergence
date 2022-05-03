%This script makes a matlab function that extracts the AVERAGE fluxes information 
%from the fluxes folder for multiple runs in a project, averages the last data points, and plots against a scanned variable
%
%Useage: fluxes(proj,flux,variable,filename,n_av,etype)
%
%proj is a string with your project name
%Scanned variable in the form of the structure created by read_gkwinput
%Filename (optional) to select which runs of project to plot (accepts
%wildcards).  Default = All.
%n_av is the (optional) number of final data points in time.dat to average
%Default is 100.  You are advised to check convergence of your runs first!
%Using the time_traces function
%
%Flux chooses the column number to look at (optional -> default = 2, Energy flux for first species).
%At present the column layout is the following: 
% 1st: particle flux, 1st species
% 2nd: heat flux, 1st species
% 3rd: momentum flux, 1st species
% 4th: particle flux 2nd species
% ...
%
%etype (optional) chooses from electro-static ('es'), electromagnetic ('em') or compressional ('bpar')
%or neoclassical (neoc, not output yet) fluxes. The default is electro-static.
%
%Example: fluxes('latest_project',1,'ROTATION.shear_rate','nl*',100,'es')
%
%This script interfaces with Yann's scripts 
%It requires gkw_path to be set up correctly for your projects directory.

function [flux_rate]=fluxes(proj,flux,variable,filename,n_av,etype)

if ~exist('flux')
	flux=2;
    disp('Energy Flux for species 1')
end

if ~exist('proj')
	proj='default';
    disp('You must provide the project name')
    %No default value
end
if ~exist('n_av')
    n_av=100;
    disp('Average over last 100 data points')
end
if ~exist('variable')
    variable='ROTATION.shear_rate';
end
if ~exist('filename')
    filename='';
end
if ~exist('etype')
    etype='es';
end

count=0;            %Counter for output files
no_files=0;         %Counter for input files

files=dir([gkwpath('input',proj) filename]);
%files=dir(gkwpath('input',proj));
total = size(files);

switch etype
    case 'es'
        fluxes_path='fluxes'
    case 'em'
        fluxes_path='fluxes_em'
    case 'bpar'
        fluxes_path='fluxes_bpar'
    case 'neoc'
        fluxes_path='fluxes_neoc'
        disp('Not applicable yet.')
        return
end

for i = 1:total(1)
    
   if(files(i).isdir==0)
        no_files=no_files+1;
        input=read_gkwinput(files(i).name,proj);
       
        if (exist([gkwpath(fluxes_path,proj) files(i).name],'file')==2)
            files(i).name;
            count = count+1;
            str=sprintf('flux_rat(count,1)=input.%s;',variable);
            eval(str);
            time=importdata([gkwpath(fluxes_path,proj) files(i).name]);
            disp(['Loaded ' gkwpath(fluxes_path,proj) files(i).name])  
            dim=size(time);
            flux_rat(count,2)=mean(time(end-n_av:end,flux))
        else
            disp(['Skipped missing', gkwpath(fluxes_path,proj) files(i).name])
        end
        %cd ..
   end
end

if (no_files==0)
    disp('No input files found')
end    
if (count==0)
    disp('No output files found')
else    
disp(sprintf('Loaded %i output files for %i input files', count, no_files))
flux_rat=sortrows(flux_rat,1);

plot(flux_rat(:,1),flux_rat(:,2),'+-');
xlabel(variable);
ylabel(sprintf('<fluxes>'));
%legend(sprintf('%s, n av = %i',proj,n_av));
end

end
