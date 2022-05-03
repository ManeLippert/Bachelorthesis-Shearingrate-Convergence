% This file defines a matlab function which plots the growth rate spectrum
% from [a batch of linear runs with nmod=1] or [a single linear run with nmod > 1].
% This script assumes the naming  conventions of the script gkwnlin,
% and requires gkwpath.m to be setup.
%
% growth_rate_spec(proj,runname,n_av,optgr)
%
% WHERE:
%
% proj and runname are strings for gkwpath.m; use wildcard in runname
% n_av is optional number of final points to average over (default 1)
% optgr = 0 (optional, default) takes set of single growth rates from time.dat
% optgr = 1 takes set of single growth rates from growth.dat
% optgr = 2 takes spectrum of growth rates from single growth.dat
% optgr = 3 takes spectrum of growth rates for eigenvalue runs
%
% EXAMPLE: grdat = growth_rate_spec('proj','runname*',5)
%
% NOTE: script parity_spec does a similar job, and gives more info 
% about the mode, and has more processing for the eigenvalue solver.


function growth_rate=growth_rate_spec(proj,filename,n_av,optgr)

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
if ~exist('optgr')
    optgr=0;
end

count=0;            %Counter for output files
no_files=0;         %Counter for input files

files=dir([gkwpath('time',proj) filename])
total = size(files);

if (optgr==2)
    input=read_gkwinput(filename,proj);
    time=load([gkwpath('other',proj) filename '/growth.dat']);
    krho=load([gkwpath('kyspec',proj)  filename '.krho']);
    for imod = 1:input.GRIDSIZE.nmod
        growth_rate(imod,1) = krho(imod,1);
        growth_rate(imod,2)=mean(time(end-n_av+1:end,imod));
    end
else
    
    for i = 1:total(1)
        
        if(files(i).isdir==0)
            no_files=no_files+1;
            input=read_gkwinput(files(i).name,proj);
            
            if (exist([gkwpath('time',proj) files(i).name],'file')==2)
                count = count+1;
                
                if (input.MODE.mode_box==1)
                    growth_rate(count,1) = input.MODE.krhomax;
                else  %use line below if mode_box=false
                    growth_rate(count,1) = input.MODE.kthrho;
                end
                %freq(count,1) = growth_rate(count,1);
                
                if (optgr == 1)
                    time=load([gkwpath('other',proj) files(i).name '/growth.dat']);
                    growth_rate(count,2)=mean(time(end-n_av+1:end,1));
                elseif(optgr == 0)
                    time=load([gkwpath('time',proj) files(i).name]);
                    growth_rate(count,2)=mean(time(end-n_av+1:end,2));
                elseif(optgr == 3)  % eigenvalue solver
                    time=load([gkwpath('time',proj) files(i).name]);
                    for k = 1:min(max(size(growth_rate,2),size(time,1)),2)
                        if (size(time,1) >= k)
                            if (time(k,2) > 0)
                                growth_rate(count,k+1)=time(k,2);
                            else
                                growth_rate(count,k+1)=-NaN;      
                            end
                        else
                            growth_rate(count,k+1)=-Inf; 
                        end
                    end 
                end
                
                disp(['Loaded ' gkwpath('time',proj) files(i).name])
                
            else
                disp(['Skipped missing', gkwpath('time',proj) files(i).name])
            end
            %cd ..
        end
        
        if (no_files==0)
            disp('No input files found')
        end
        if (count==0)
            disp('No output files found')
        else
            disp(sprintf('Loaded %i output files for %i input files', count, no_files))
            growth_rate=sortrows(growth_rate,1);
            %freq=sortrows(freq,1);
        end
    end    
    
end

    plot(growth_rate(:,1),growth_rate(:,2)*sqrt(2),'+-','DisplayName',[proj '.' filename]);
    if (optgr == 3) 
    hold all    
      plot(growth_rate(:,1),growth_rate(:,3),'o-','DisplayName',[proj '.' filename]);
    end 
    xlabel('k_\theta \rho_i');
    ylabel('\gamma (R/v_{th})');
    ylabel('\gamma (R/c_s)');
end
