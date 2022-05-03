% Read the fluxes from a GKW run 
% No conversions are performed, the data is read as it appears in the file
%
% If boundary damping is used, the fluxes are also read from the 
% radial profile and averaged over the further than 2*bwidth from the boundary
% 
% Usage [out] = read_gkwflux(flnm,proj,in); 
%
% INPUT
% flnm : the file name 
% proj : (optional) the project name 
% in   : (optional) the name to which the structure is added 
% 
% OUTPUT
% out.pfluxes(n) : Electro static particle flux (n species)
% out.qfluxes(n) : Electro static energy flux (n species)
% out.vfluxes(n) : Electro static momentum flux (n species)

function [out]=read_gkwflux(flnm,proj,in)

  % If an input structure is given copy it to the output 
  if (exist('in'))
    out  = in;
  end

  % Set the project if not given 
  if ~exist('proj') 
	proj=[];
  end

  % Always read the time in as well 
  flpth=gkwpath('time',proj);
  name = [flpth flnm]; 
  if (exist(name,'file')) 
    out.time = load(name); 
  else
     out.time = 0.
  end
  
  % Set the path 
  flpth=gkwpath('fluxes',proj);

  % Read the electro-static output file 
  name = [flpth flnm]; 

  if (exist(name,'file')) 

    fluxes = load(name); 

    nspec = size(fluxes,2)/3; 
    for i = 1: nspec 
      out.pfluxes(:,i) = fluxes(:,3*(i-1) + 1);
      out.qfluxes(:,i) = fluxes(:,3*(i-1) + 2); 
      out.vfluxes(:,i) = fluxes(:,3*(i-1) + 3); 
    end; 

  else 
    disp(['File does not exist: ' name])  
    out.pfluxes = 0.;
    out.qfluxes = 0.;
    out.vfluxes = 0.;
  end; 

  % Read the electro-magnetic bit 
  name = [flpth 'em/' flnm]; 

  if (exist(name,'file')) 

    fluxes = load(name); 

    for i = 1: nspec 
      out.pfluxem(:,i) = fluxes(:,3*(i-1) + 1);
      out.qfluxem(:,i) = fluxes(:,3*(i-1) + 2); 
      out.vfluxem(:,i) = fluxes(:,3*(i-1) + 3); 
    end; 

  else 

    disp(['File does not exist: ' name])  
    out.pfluxem = 0.;
    out.qfluxem = 0.;
    out.vfluxem = 0.;

  end; 

  % Read the part connected with the field compression 
  name = [flpth 'bpar/' flnm]; 

  if (exist(name,'file')) 

    fluxes = load(name); 

    for i = 1: nspec 
      out.pfluxbpar(:,i) = fluxes(:,3*(i-1) + 1);
      out.qfluxbpar(:,i) = fluxes(:,3*(i-1) + 2); 
      out.vfluxbpar(:,i) = fluxes(:,3*(i-1) + 3); 
    end; 

  else 

    out.pfluxbpar = 0.;
    out.qfluxbpar = 0.;
    out.vfluxbpar = 0.;

  end;

  % read radial flux profiles if boundary damping is present
  if (exist('in') )
    if isfield(in.CONTROL,'spectral_radius');
      if (in.CONTROL.spectral_radius(1)=='f')
          flpth=gkwpath('radial',proj);
          dat = load([flpth '/pflux_es/' flnm]);
          dat2 = load([flpth '/eflux_es/' flnm]);
          dat3 = load([flpth '/vflux_es/' flnm]);

          out.xphi=load([gkwpath('xphi',proj) flnm]);
          out.yphi=load([gkwpath('yphi',proj) flnm]);

          %average over region away from boundary layer
          for i = 1: nspec
              is=1+ (i-1)*in.GRIDSIZE.nx; 
              ie=i*in.GRIDSIZE.nx; 

              out.pfluxes_rad(:,:,i) = dat(:,is:ie);
              out.qfluxes_rad(:,:,i) = dat2(:,is:ie);
              out.vfluxes_rad(:,:,i) = dat3(:,is:ie);

              is= is + 3*abs(in.KROOK.bwidth);
              ie= ie - 3*abs(in.KROOK.bwidth);

              %remove bounday layer and make a radial average
              is=1+ (i-1)*in.GRIDSIZE.nx + abs(2*in.KROOK.bwidth);
              ie=i*in.GRIDSIZE.nx - abs(2*in.KROOK.bwidth);

              out.pfluxes_rav(:,i) = mean(dat(:,is:ie),2);
              out.qfluxes_rav(:,i) = mean(dat2(:,is:ie),2);
              out.vfluxes_rav(:,i) = mean(dat3(:,is:ie),2);              
          end

      end
    end
  end
end



