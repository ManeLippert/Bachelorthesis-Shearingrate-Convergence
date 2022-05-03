% matlab function that extracts the growth rate information from the time folder
% for multiple (linear) runs in a project and plots against a scanned variable
%
% Usage: growth_rates(proj,,file,variable,n_av,optgs2)
%
% proj is a string with your project name
% Scanned variable in the form of the structure created by read_gkwinput
% n_av is the (optional) number of final data points in time.dat to average
% Default is 1.  For most linear runs time averaging should not be necssary. You are adivsed to check convergence of your runs first!
%
% Example: growth_rates('latest_project','runs*','ROTATION.shear_rate',100,optgs2)
%
% setting optgs2 = 1 reads the growth rate from the equivalent gs2 project
% (but the GKW project must still be present)
% This script interfaces with Yann's scripts and the gs2tools
% It requires gkw_path to be set up correctly for your projects directory.
% gs2 projects path is hard coded

function [growth_rat]=growth_rates(proj,filename,variable,n_av,optgs2)


if ~exist('proj')
    proj='default';
    disp('You must provide the project name')
    %No default value
end
if ~exist('variable')
    variable='ROTATION.shear_rate';
end
if ~exist('n_av')
    n_av=1;
end

if ~exist('optgs2')
    optgs2=0;
end

plotsymbol='o';
gs2runs='~/runs/gs2/';
count=0;            %Counter for output files
no_files=0;         %Counter for input files

%Currently you have to give a prefix for the filenames.
%files = dir('input/*');
files=dir([gkwpath('input',proj) filename]);
%Could be made neater to only list all files (not directories)
total = size(files);

for i = 1:total(1)
    
    if(files(i).isdir==0)
        no_files=no_files+1;
        input=read_gkwinput(files(i).name,proj);
        
        if (optgs2==0)
            
            if (exist([gkwpath('time',proj) files(i).name],'file')==2)
                count = count+1;
                str=sprintf('growth_rat(count,1)=input.%s;',variable);
                eval(str);
                time=load([gkwpath('time',proj) files(i).name]);
                disp(['Loaded ' gkwpath('time',proj) files(i).name])
                dim=size(time);
                growth_rat(count,2)=mean(time(end-n_av+1:end,2));
            else
                disp(['Skipped missing', gkwpath('time',proj) files(i).name])
            end
            
        else
            
            if (exist([gs2runs proj '/output/' files(i).name],'file'))
                count = count+1;
                str=sprintf('growth_rat(count,1)=input.%s;',variable);
                eval(str);                
                [aky, omavg, qheat, pflux, vflux, qmheat, pmflux] = read_gs2output_adia(files(i).name,[gs2runs proj '/output/'],0)
                growth_rat(count,2)=omavg(1,2)/sqrt(2);
                plotsymbol='-';
                %out(i).freq=omavg(1,1)/sqrt(2);
            else
                disp([gs2runs proj '/output/' files(i).name ' not found'])
                continue
            end
            
        end
        %cd ..
    end
end

disp(sprintf('Loaded %i output files for %i input files', count, no_files))
if (~exist('growth_rat','var'))
  growth_rat=[NaN NaN];  
else    
  growth_rat=sortrows(growth_rat,1);
end

plot(growth_rat(:,1),growth_rat(:,2),plotsymbol,'DisplayName',filename);
xlabel(variable);
ylabel('\gamma (R/v_{th})');
%legend(sprintf('%s, n av = %i',proj,n_av));

end
