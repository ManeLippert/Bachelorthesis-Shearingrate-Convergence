% This file defines a matlab function which plots the [growth rate / freq / parity] spectrum
% from a batch of linear runs with nmod=1. This script assumes the naming
% conventions of the script gkwnlin, and requires gkw_path.m to be setup.
% 
% parity_spec(proj,runname,n_av)
%
% WHERE:
%
% proj and runname are strings for gkwpath.m; use wildcard in runname
% n_av is optional number of final points to average over (default 1)
%
% If n_av <= 0, use values from the eigenvalue solver
%
% EXAMPLE: parity_spec('proj','runname*',5)

function [growth_rate]=parity_spec(proj,filename,n_av)

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
if ~exist('n_av','var')
   n_av=1;
end

%rova=1.36;
rova=1;

count=0;            %Counter for output files
no_files=0;         %Counter for input files

files=dir([gkwpath('time',proj) filename]);
total = size(files);
growth_rate(1,1)=0.0;

for i = 1:total(1)
    
   if(files(i).isdir==0)
        no_files=no_files+1;
        input=read_gkwinput(files(i).name,proj);
               
        if (exist([gkwpath('time',proj) files(i).name],'file')==2)
            time=load([gkwpath('time',proj) files(i).name]);
            disp(['Loaded ' gkwpath('time',proj) files(i).name])
            
            %disp([input.MODE.mode_box])
            %disp(std(time(end-n_av:end,2))/abs(mean(time(end-n_av:end,2))))
            
            if (size(time,1)<n_av) continue; end
                
            count = count+1;  
                     
            if  (n_av > 0)
                if (std(time(end-n_av:end,2))/abs(mean(time(end-n_av:end,2))) < 0.01 && mean(time(end-n_av:end,2)) < 10)
                    growth_rate(count,2)=mean(time(end-n_av+1:end,2));
                    growth_rate(count,3)=mean(time(end-n_av+1:end,3));
                    [p,e]=parity(proj,files(i).name);
                    growth_rate(count,4)=p;
                    growth_rate(count,5)=e;
                else
                    growth_rate(count,2)=NaN;
                    growth_rate(count,3)=NaN;
                    growth_rate(count,4)=NaN;
                    growth_rate(count,5)=NaN;
                end
            else                
                time=load([gkwpath('time',proj) files(i).name]);
                for k = 1:min(max(size(growth_rate,2),size(time,1)),2)
                    if (size(time,1) >= k)
                        if (time(k,2) > 0)
                            try
                              [p,e]=parity(proj,files(i).name,-k);                              
                            catch
                               p=NaN; e=NaN; 
                            end
                              growth_rate(count,2,k)=time(k,2);
                              growth_rate(count,3,k)=time(k,3);
                              growth_rate(count,4,k)=p;
                              growth_rate(count,5,k)=e;
                        else
                            growth_rate(count,2,k)=-NaN;
                            growth_rate(count,3,k)=-NaN;
                            growth_rate(count,4,k)=-NaN;
                            growth_rate(count,5,k)=-NaN;
                        end
                    else
                        growth_rate(count,2,k)=-Inf;
                        growth_rate(count,3,k)=-Inf;
                        growth_rate(count,4,k)=-Inf;
                        growth_rate(count,5,k)=-Inf;
                    end
                end
            end
            
            if (input.MODE.mode_box==1)            
                growth_rate(count,1,:) = input.MODE.krhomax;
            else  %use line below if mode_box=false
                growth_rate(count,1,:) = input.MODE.kthrho;
            end
            %growth_rate(count,1,:) = input.SPCGENERAL.beta_ref;
            
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

if (n_av > 0)
    
    growth_rate=sortrows(growth_rate,1);
    
    subplot(4,1,1)
    plot(growth_rate(:,1),growth_rate(:,2)/rova,'+-','DisplayName',[proj '.' filename]);
    xlabel('k_\theta \rho_i');
    ylabel('\gamma (a/v_{th})');
    grid on
    hold all
    
    subplot(4,1,2)
    plot(growth_rate(:,1),growth_rate(:,3)/rova,'+-','DisplayName',[proj '.' filename]);
    xlabel('k_\theta \rho_i');
    ylabel('\omega (a/v_{th})');
    grid on
    hold all
    
    subplot(4,1,3)
    plot(growth_rate(:,1),growth_rate(:,4)/rova,'+-','DisplayName',[proj '.' filename]);
    xlabel('k_\theta \rho_i');
    ylabel('parity p');
    grid on
    hold all
    
    subplot(4,1,4)
    plot(growth_rate(:,1),growth_rate(:,5)/rova,'+-','DisplayName',[proj '.' filename]);
    xlabel('k_\theta \rho_i');
    ylabel('EM fraction');
    grid on
    hold all
    
else
    dat=growth_rate;
    good=find(abs(dat(:,2,1))>0.01 & abs(dat(:,2,1))<1.0);  
    good2=find(abs(dat(:,2,2))>0.01 & abs(dat(:,2,2))<1.0);  
    itg=find(abs(dat(good,4,1))<0.5);
    itg2=find(abs(dat(good2,4,2))<0.5);    
    mtm=find(abs(dat(good,4,1))>0.5);    
    mtm2=find(abs(dat(good2,4,2))>0.5);
    
    datall=[dat(:,:,1); dat(:,:,2)];
    datall=sortrows(datall,1);
    goodall=find(abs(datall(:,2))>0.0001 & abs(datall(:,2))<1.0); 
    datall=datall(goodall,:);
    mtmall=find(abs(datall(:,4))>0.5);
    itgall=find(abs(datall(:,4))<0.5);
        
    %figure
    %plot(dat(itg,1,1),dat(itg,2,1),'o-','DisplayName',[proj '.' filename])
    hold all
    %plot(dat(mtm,1,1),dat(mtm,2,1),'v-','MarkerSize',8,'DisplayName',[proj '.' filename])       
    %plot(dat(good2,1,2),dat(good2,2,2),'-','DisplayName',[proj '.' filename])
    %plot(dat(itg2,1,2),dat(itg2,2,2),'o','DisplayName',[proj '.' filename])
    %plot(dat(mtm2,1,2),dat(mtm2,2,2),'v','MarkerSize',8,'DisplayName',[proj '.' filename])
    
    plot(datall(itgall,1),datall(itgall,2),'o-','DisplayName',[proj '.' filename])
    plot(datall(mtmall,1),datall(mtmall,2),'v-','MarkerSize',8,'DisplayName',[proj '.' filename])
    
    xlabel('k_\theta \rho_i');
    ylabel('\gamma (R/v_{th})');
    %ylim([0 1.1*max(dat(:,2,1))])
    %xlim([0 max(dat(:,1,1))])   
       
end


end
