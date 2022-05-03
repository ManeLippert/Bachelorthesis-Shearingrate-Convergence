function out=read_gkw_all(file,proj,optnstart,optnend,optamin,optBref)
% Attempts to read everything possible from a project folder
% into a the ULTIMATE single structure for a GKW NL run
%
% Usage: out=read_gkw_all(file,proj,optnstart,optnend,optamin,optBref)
%
% optnstart = timestep to begin averaging (default: mid point)
% optnend   = timestep to end averaging (default: end of run)
%
% Calculates the transport coefficients in GKW units, GB units, and SI units
% Uses values of Tref, nref and Rref from COLLISIONS for conversions
% Calculation of SI units requires Bref (from hamada file, or optBref)
% Calculation of GB units requires rminor a (from hamada file, or optamin)
%
% WARNING: The [q/v]fluxes in the structure are not identical to those 
% in the gkw fluxes file, here they are corrected by T_R/v_R factors
%
% WARNING: The dimensional values should be used with caution, 
% and the user should understand how they are calculated: 
% It is advised to perform independent checks on the numbers returned.

%Read input.out file for all namelists
out=read_gkwinput(file,proj,1);
out.file=file;
out.proj=proj;

% Alternate way to provide the dVdr information
% Read input file for comments in GEOM_PLUS namelist
% in=read_gkwinput(file,proj,0);
% if isfield(in,'GEOM_PLUS')
%   out.GEOM_PLUS=in.GEOM_PLUS;
% else
%   out.GEOM_PLUS.dvdr=NaN; 
% end

if exist('opta','var'); out.amin=optamin; end
if exist('optBref','var'); out.Bref=optBref; end  
if exist('optnstart','var');  out.nstart=optnstart; end
if exist('optnend','var');  out.nend=optnend; end

out=read_gkwflux(file,proj,out);

% Use radial fluxes subset instead of averaged flux
%out.qfluxes=out.qfluxes_rav;

try
  %out.hamada=read_hamada(out.GEOM.eqfile(14:end),gkwpath('hamada',proj));
  out.hamada=read_hamada(out.GEOM.eqfile(strfind(out.GEOM.eqfile,'hamada')+7:end),gkwpath('hamada',proj)); 
  out.Bref=out.hamada.b0exp;
  out.amin=out.hamada.amin(end);
catch    
   disp('failed to find hamada')    
end

try
  %out.geom=read_geom_plus(file,proj);
  out.geom=read_geom(file,proj);
  
 
  % If hamada is present, find rgeo and surface index
  if isfield(out,'hamada')
    ind=find(out.hamada.q<out.geom.q+0.001);
    out.hamada.ind=ind(end);
    out.geom.rgeo=out.hamada.rgeom(out.hamada.ind);
  end
  
catch;
  disp('failed to find geom')  
end

% Calculate the transport coefficents: 
% Modifies the values in [q/v]fluxes by T_R/v_R factors.
% Use these outputs with care, check how they are calculated
%try
  out=trans_coef(out);
%end 

% Always read the krho grid in as well 
try
  krho = load([gkwpath('krho',proj) file]);
  out.krho=krho(:,1)';
end

if (out.MODE.mode_box(1)=='t')
  try
    out=read_gkwspectra(file,proj,out);
  end
  
  if(out.CONTROL.spectral_radius(1)=='f')
      out.kxrh_min=2*pi/out.GRIDSIZE.lx;
  else
      out.kxrh_min=out.kxrh((out.GRIDSIZE.nx+3)/2);
  end
  
  out.geom.kxnorm=sqrt(out.geom.g_eps_eps(out.GRIDSIZE.n_s_grid/2));
  out.geom.ikxspace=out.krho(2)*out.geom.shat*out.geom.q/out.geom.eps/out.kxrh_min/out.geom.kthnorm;

end

end