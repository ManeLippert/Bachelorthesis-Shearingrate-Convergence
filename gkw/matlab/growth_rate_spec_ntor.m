% This file defines a matlab function which outputs the linear spectra
% against toroidal mode number.  Uses read_gkw_all to get rhostar from inputs
% Can read the data from either
% [a batch of linear runs with nmod=1] or [a single linear run with nmod > 1].
% This script assumes the naming  conventions of the script gkwnlin,
% and requires gkwpath.m to be setup.
%
% out = growth_rate_spec_ntor(proj,runname,n_av,optgr)
%
% WHERE:
%
% proj and runname are strings for gkwpath.m; use wildcard in runname
% n_av is optional number of final points to average over (default 1)
% optgr = 0 (optional, default) takes set of single growth rates from time.dat
% optgr = 1 takes set of single growth rates from growth.dat
% optgr = 2 takes spectrum of growth rates from single growth.dat
%
% OUTPUT:
%
% Column 1 is toroidal mode number
% Column 2 is growth rate in c_s/a
% Column 3 is mode frequency in c_s/a (negative is electron direction)
% Column 4 is the parity parameter (absolute value > 0.5 for microtearing)
% Column 5 is the QL ion flux, normalised to (|phi|^2 + |Apar|^2)
% Column 6 is the QL electron flux, normalised to (|phi|^2  + |Apar|^2)
%
%
% EXAMPLE: grdat = growth_rate_spec_ntor('proj','runname*',5,0)

function [growth_rate]=growth_rate_spec_ntor(proj,filename,n_av,optgr)

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

disp(files(1).name)
dat=read_gkw_all(files(1).name,proj,1);

grfac=1/(dat.rova*sqrt(dat.SPECIES(2).temp/2));

% ntor 
kyfac=1/(2*pi*dat.rhostar*dat.geom.kthnorm);   

% to gs2 ky . rho_ref, where rho_ref = rho_i
% kyfac2 =dat.geom.e_eps_zeta*2/(dat.geom.kthnorm));

if (optgr==2)
    input=read_gkwinput(filename,proj);
    time=load([gkwpath('other',proj) filename '/growth.dat']);
    freq=load([gkwpath('other',proj) filename '/frequencies.dat']);
    
    krho=load([gkwpath('kyspec',proj)  filename '.krho']);
    for imod = 1:input.GRIDSIZE.nmod
        growth_rate(imod,1) = krho(imod,1);
        growth_rate(imod,2)=mean(time(end-n_av+1:end,imod));
        growth_rate(imod,3)=mean(freq(end-n_av+1:end,imod));        
        %growth_rate(imod,4) = parity(proj,filename,imod);
    end
else
    
    for i = 1:total(1)     
                
        if(files(i).isdir==0)
            no_files=no_files+1;
            disp(files(i).name)
            input=read_gkwinput(files(i).name,proj);
            
            if (exist([gkwpath('time',proj) files(i).name],'file')==2)
                count = count+1;
                
                if (optgr == 1)
                    time=load([gkwpath('other',proj) files(i).name '/growth.dat']);
                    growth_rate(count,2)=mean(time(end-n_av+1:end,1));
                    
                    freq=load([gkwpath('other',proj) files(i).name '/frequencies.dat']);
                    growth_rate(count,3)=mean(freq(end-n_av+1:end,1));                    
                    
                else
                    time=load([gkwpath('time',proj) files(i).name]);
                    fluxes=load([gkwpath('fluxes',proj) files(i).name]);
                    
                    nav = n_av;
                    if (size(fluxes,1) > 3000)
                       nav = 100
                    end
                    
                    growth_rate(count,2)=mean(time(end-nav+1:end,2));
                    growth_rate(count,3)=mean(time(end-nav+1:end,3));
                    growth_rate(count,5)=mean(fluxes(end-nav+1:end,2));
                    growth_rate(count,6)=mean(fluxes(end-nav+1:end,5));
                end
                
                disp(['Loaded ' gkwpath('time',proj) files(i).name])
                
                growth_rate(count,4) = parity(proj,files(i).name); 
                
                %disp([input.MODE.mode_box])
                
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
        
        if (no_files==0)
            disp('No input files found')
        end
        if (count==0)
            disp('No output files found')
        else
            disp(sprintf('Loaded %i output files for %i input files', count, no_files))
            growth_rate=sortrows(growth_rate,1);
            
        end
    end    
    
end

   %growth_rate(:,5)=growth_rate(:,1)*kyfac2;
   growth_rate(:,1)=growth_rate(:,1)*kyfac;
   growth_rate(:,2)=growth_rate(:,2)*grfac;
   growth_rate(:,3)=growth_rate(:,3)*grfac;    

   loglog(growth_rate(:,1),growth_rate(:,2),'+-','DisplayName',[proj '.' filename]);
   xlabel('n');
   ylabel('\gamma (a/c_s)');
end
