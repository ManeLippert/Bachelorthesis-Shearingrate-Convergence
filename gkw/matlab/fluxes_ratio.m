%This script makes a matlab function that extracts the RATIO of AVERAGE fluxes information 
%from the fluxes folder for multiple runs in a project, averages the last data points, and plots against a scanned variable
%
%Usage: out=fluxes_ratio(proj,flux1,flux2,variable,filename,n_av,etype,etype2,t_av)
%
%proj is a string with your project name
%Scanned variable in the form of the structure created by read_gkwinput
%Filename (optional) to select which runs of project to plot (accepts
%wildcards).  Default = All.
%n_av is the (optional) number of final data points in time.dat to average (default 100 if not provide)
%t_av is optional time to start averaging from (overrides n_av)
%You are advised to check convergence of your runs first using the time_traces function
%
%Flux chooses the column number to look at
%Flux2 chooses the column number to normalize with
%
%etype and etype2 chooses from electro-static ('es'), electromagnetic ('em')
%or compressional ('bpar') or neoclassical (neoc, not output yet) fluxes.
%The default is electro-static.
%
%Example: fluxes_ratio('latest_project',1,2,'ROTATION.shear_rate','files*',100,'em','es')
%
%This script requires gkw_path to be set up correctly for your projects directory.

function [flux_rat]=fluxes(proj,flux,flux2,variable,filename,n_av,etype,etype2,t_av)

if ~exist('flux')
	flux=2;
    disp('Energy Flux for species 1')
end

if ~exist('flux2')
	flux=5;
    disp('Over Energy flux for species 2')
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
if ~exist('etype2')
    etype2='es';
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

switch etype2
    case 'es'
        fluxes_path2='fluxes'
    case 'em'
        fluxes_path2='fluxes_em'
    case 'bpar'
        fluxes_path2='fluxes_bpar'
    case 'neoc'
        fluxes_path2='fluxes_neoc'
        disp('Not applicable yet.')       
        return
end

for i = 1:total(1)
    
   if(files(i).isdir==0)
        no_files=no_files+1;
        input=read_gkwinput(files(i).name,proj);
       
        if ( exist([gkwpath(fluxes_path,proj) files(i).name],'file')==2 && ...
             exist([gkwpath(fluxes_path2,proj) files(i).name],'file')==2 )
        
            files(i).name;
            count = count+1;
            str=sprintf('flux_rat(count,1)=input.%s;',variable);
            eval(str);
            file=importdata([gkwpath(fluxes_path,proj) files(i).name]);
            
            if exist('t_av')
              time=load([gkwpath('time',proj) files(i).name]);
              ind=find(time(:,1) > t_av);
            else  
              ind=[length(file)-n_av:1:length(file)];
            end
            
            disp(['Loaded ' gkwpath(fluxes_path,proj) files(i).name])
            if (strcmp(etype,etype2))
                file2=file;
            else
                file2=importdata([gkwpath(fluxes_path2,proj) files(i).name]);
                disp(['Loaded ' gkwpath(fluxes_path2,proj) files(i).name])                
            end
            dim=size(file);
            dim2=size(file2);            
            flux_rat(count,2)=0;
            flux_rat(count,3)=0;
            a=mean(file(ind,flux));
            b=mean(file2(ind,flux2));
            err_a=std(file(ind,flux))/sqrt(n_av);
            err_b=std(file2(ind,flux2))/sqrt(n_av);
            flux_rat(count,2) = a/b;
            %Standard error of ratios
            flux_rat(count,3) = (a/b)*sqrt((err_a/a)^2+(err_b/b)^2);
        elseif ( exist([gkwpath(fluxes_path,proj) files(i).name],'file')~=2 && ...
                 exist([gkwpath(fluxes_path2,proj) files(i).name],'file')==2 )
            disp(['Skipped missing', gkwpath(fluxes_path,proj) files(i).name])
        elseif ( exist([gkwpath(fluxes_path,proj) files(i).name],'file')==2 && ...
                 exist([gkwpath(fluxes_path2,proj) files(i).name],'file')~=2 )
            disp(['Skipped missing', gkwpath(fluxes_path2,proj) files(i).name])            
        else
            disp(['Skipped missing', gkwpath(fluxes_path,proj) files(i).name, ...
                                     gkwpath(fluxes_path2,proj) files(i).name])            
        end
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

%renormalise if desired

hold all;
errorbar(flux_rat(:,1),flux_rat(:,2),flux_rat(:,3),'+-','DisplayName',[proj filename]);
xlabel(variable);
ylabel(['flux' sprintf('%i',flux) '/ flux' sprintf('%i',flux2)]);
end

end
