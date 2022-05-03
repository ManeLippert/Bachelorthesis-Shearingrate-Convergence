% Read the ky spectra from a GKW run and time average
% No conversions are performed, the data is read as it appears in the file
%
% Usage [out] = read_gkwspectra(flnm,proj,in); 
%
% INPUT
% flnm : the file name 
% proj : (optional) the project name 
% in   : (optional) the name to which the structure is added 


function [out]=read_gkwspectra(file,proj,in)

  % If an input structure is given copy it to the output 
  if (exist('in'))
    out  = in;
  end

  % Set the project if not given 
  if ~exist('proj') 
	proj=[];
  end

  % Always read the krho grid in as well 
  krho = load([gkwpath('krho',proj) file]);
  out.krho=krho(:,1)';
  
  try
    kxrh = load([gkwpath('kxrh',proj) file]);
    out.kxrh=kxrh(1,:);
  end
  
  
  if isfield(out,'rhostar')    
    kyfac=1/(2*pi*out.rhostar*out.geom.kthnorm);  
    out.ntor=out.krho*kyfac;    
    %out.kygene = out.krho/(out.geom.kthnorm*sqrt(out.SPECIES(2).temp/2));      
  end  
  
  out.kyspec = load([gkwpath('kyspec',proj) file]);
  out.kyspec_em=load([gkwpath('kyspec_em',proj) file]);
  out.kyeflux=load([gkwpath('eflux_spec',proj) file]);
  
  if isfield(out,'nstart')
           
     out.kyspec_av=squeeze(mean(out.kyspec(out.nstart:out.nend,:),1));
     out.kyspec_em_av=squeeze(mean(out.kyspec_em(out.nstart:out.nend,:),1));
     
     modes=out.GRIDSIZE.nmod;
  
     for sp=1:out.GRIDSIZE.number_of_species
  
       start=modes*(sp-1)+1;
       fin=modes*sp;
    
       out.kyeflux_av(sp,1:modes) = mean(out.kyeflux(out.nstart:out.nend,start:fin),1)*out.SPECIES(sp).temp       
        
     end    
     
  end  
  
   if isfield(out,'fluxMW')
      out.specout=[out.ntor' out.kyspec_av' out.kyspec_em_av' out.kyeflux_av(1,:)'*out.fluxMW(1) out.kyeflux_av(2,:)'*out.fluxMW(2)];
   end
   
end



