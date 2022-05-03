% This file defines a matlab function which plots the dispersion relation
% from a batch of linear runs with nmod=1. This script assumes the naming
% conventions of the script gkwnlin, and requires gkw_path.m to be setup.
% 
% freq_spec(proj,runname,n_av)
%
% WHERE:
%
% proj and runname are strings for gkwpath.m; use wildcard in runname
% n_av is optional number of final points to average over (default 1)
%
% EXAMPLE: freq_spec('proj','runname*',5)

function [growth_rate]=freq_spec(proj,filename,n_av)

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
if ~exist('n_av')
   n_av=1;
end

%rova=1.36;
rova=1;

count=0;            %Counter for output files
no_files=0;         %Counter for input files

files=dir([gkwpath('time',proj) filename]);
total = size(files);

for i = 1:total(1)
    
   if(files(i).isdir==0)
        no_files=no_files+1;
        input=read_gkwinput(files(i).name,proj);
               
        if (exist([gkwpath('time',proj) files(i).name],'file')==2)
            count = count+1;
            time=load([gkwpath('time',proj) files(i).name]);
            disp(['Loaded ' gkwpath('time',proj) files(i).name])
            
            %disp([input.MODE.mode_box])
            growth_rate(count,2)=mean(time(end-n_av+1:end,2));
            growth_rate(count,3)=mean(time(end-n_av+1:end,3));
            if (input.MODE.mode_box==1)            
                growth_rate(count,1) = input.MODE.krhomax;
            else  %use line below if mode_box=false
                growth_rate(count,1) = input.MODE.kthrho;
            end
            
         else
            disp(['Skipped missing', gkwpath('time',proj) files(i).name])
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
growth_rate=sortrows(growth_rate,1);

subplot(2,1,1)
set(0,'defaulttextinterpreter','tex')
set(gca,'box','on','fontsize',12,'Xminortick','on','Yminortick','on')
set(gca,'TickLength',[0.015,0.07])


plot(growth_rate(:,1),growth_rate(:,2)/rova,'+-','DisplayName',[proj '.' filename]);
xlabel('k_\theta \rho_i');
ylabel('\gamma (R/v_{th})');
grid on
hold all

subplot(2,1,2)
set(0,'defaulttextinterpreter','tex')
set(gca,'box','on','fontsize',12,'Xminortick','on','Yminortick','on')
set(gca,'TickLength',[0.015,0.07])


plot(growth_rate(:,1),growth_rate(:,3)/rova,'+-','DisplayName',[proj '.' filename]);
xlabel('k_\theta \rho_i');
ylabel('\omega (R/v_{th})');
grid on
hold all

end
