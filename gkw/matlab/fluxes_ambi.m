% This function plots flux of charge for each species and checks ambipolarity
% In, GKW ambipolarity must holds in all cases, including with poloidal asymmetries
% (and rotation gradients, which is the most complex case).  With strong rotation,
% this requires the same rotation gradient for all non trace species.
%
% Reads data from the fluxes folder for multiple runs in a project,
% averages the last data points, and plots against a scanned variable.
% Plots first 4 species if possible
%
% Usage: fluxes_ambi(proj,variable,filename,n_av,etype)
%
% proj is a string with your project name
% Scanned variable in the form of the structure created by read_gkwinput
% Filename (optional) to select which runs of project to plot (accepts wildcards).  Default = All.
% n_av is the (optional) number of final data points in fluxes.dat to average
% Default is 100.  You are advised to check convergence of your runs first!
% Using the time_traces function
% etype (optional) chooses from electro-static ('es'), electromagnetic ('em') or compressional ('bpar')
% or neoclassical (neoc, not output yet) fluxes. The default is electro-static.
%
% Example: fluxes_ambi('latest_project','SPECIES(4).dens','nl*',3)
%
% This script requires gkw_path to be set up correctly for your projects directory.

function [flux_rate]=fluxes_ambi(proj,variable,filename,n_av,etype)

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
        input=read_gkwinput(files(i).name,proj,0);
       
        if (exist([gkwpath(fluxes_path,proj) files(i).name],'file')==2)
            files(i).name;
            count = count+1;
            
            % Normalising density for fluxes is always LFS.  In some cases, if input density is FSA
            % we need to reconstruct the normalising desntiy using 1/e_0
            if (input.ROTATION.cf_trap == 1)
                if (strcmp(input.GEOM.r0_loc,'FSA'))
                    disp(input.GEOM.r0_loc)
                    try
                        fsadat=read_gkwcffsa(files(i).name, gkwpath('root',proj));
                    catch
                        error('Need the cffsa.dat file with CF effects if R0_loc = FSA')
                    end
                elseif(strcmp(input.GEOM.r0_loc,'LFS'))
                    fsadat.e_0=[1.0 1.0 1.0 1.0]
                else
                    error('Unknown R0_loc - may not be explicit in input file, try input.out?')
                end
            else
                fsadat.e_0=[1.0 1.0 1.0 1.0]
            end
            
            str=sprintf('flux_rat(count,1)=input.%s;',variable);
            eval(str);
            %flux_rat(count,1)
            
            
            time=importdata([gkwpath(fluxes_path,proj) files(i).name]);
            disp(['Loaded ' gkwpath(fluxes_path,proj) files(i).name])  
            dim=size(time);        
           
            ambi(count,1) = 0.0;
            
            for l = 1:input.GRIDSIZE.number_of_species
               % Normalising density for fluxes is always LFS.  In some cases, if input density is FSA
               % we need to reconstruct the normalising desntiy using 1/e_0
               flux_rat(count,1+l)=mean(time(end-n_av:end,(l-1)*3+1))*input.SPECIES(l).dens*input.SPECIES(l).z/fsadat.e_0(l);
               ambi(count,1) = ambi(count,1) + flux_rat(count,1+l);
            end
            
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
%flux_rat=sortrows(flux_rat,1);
% 
% plot(flux_rat(:,1)/fsadat.e_0(4),ambi(:,1),'+-','DisplayName','total');
% hold all
% 
% plot(flux_rat(:,1)/fsadat.e_0(4),flux_rat(:,2),'+-','DisplayName','e');
% plot(flux_rat(:,1)/fsadat.e_0(4),flux_rat(:,3),'+-','DisplayName','D');
% plot(flux_rat(:,1)/fsadat.e_0(4),flux_rat(:,4),'+-','DisplayName','Be');
% plot(flux_rat(:,1)/fsadat.e_0(4),flux_rat(:,5),'+-','DisplayName','W');
% 
% xlabel('W conc (LFS)');



plot(flux_rat(:,1),ambi(:,1),'+-','DisplayName','total');
hold all

plot(flux_rat(:,1),flux_rat(:,2),'+-','DisplayName',['Z =  ' num2str(input.SPECIES(1).z)]);
plot(flux_rat(:,1),flux_rat(:,3),'+-','DisplayName',['Z =  ' num2str(input.SPECIES(2).z)]);
try
  plot(flux_rat(:,1),flux_rat(:,4),'+-','DisplayName',['Z =  ' num2str(input.SPECIES(3).z)]);
end
try
  plot(flux_rat(:,1),flux_rat(:,5),'+-','DisplayName',['Z =  ' num2str(input.SPECIES(4).z)]);
end

xlabel('W conc (FSA)');


ylabel(sprintf('Flux of charge: Gamma_s*Z_s (R_{ref}/Vth_{ref} rho_*^2)'));
%legend(sprintf('%s, n av = %i',proj,n_av));
end

end
