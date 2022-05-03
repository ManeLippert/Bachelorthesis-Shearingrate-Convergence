% Function to calculate derived paramaters for GKW mode_box inputs
% Returns structure with derived mode quantities
%
% Usage:
%   out=mode_box_calculator(nmod,nx,krhomax,ikxspace,shat,eps,q)
%
% Alternatively takes a gkw input file structure as a single input
%

function [out in]=mode_box_calculator(varargin)

  if (nargin==1)
    strin=cell2mat(varargin(1));
    
    in.nmod     = strin.GRIDSIZE.nmod;
    in.nx       = strin.GRIDSIZE.nx;
    in.krhomax  =  strin.MODE.krhomax;
    in.ikxspace =  strin.MODE.ikxspace;
    in.shat     =  strin.GEOM.shat;
    in.eps      =  strin.GEOM.eps;
    in.q        =  strin.GEOM.q;
    
  elseif (nargin==7)
    in.nmod     = cell2mat(varargin(1));
    in.nx       = cell2mat(varargin(2));
    in.krhomax  = cell2mat(varargin(3));
    in.ikxspace = cell2mat(varargin(4));
    in.shat     = cell2mat(varargin(5));
    in.eps      = cell2mat(varargin(6));
    in.q        = cell2mat(varargin(7));
  else
    in=NaN;
    out=NaN;
    disp('Incorrect number of input arguments')
    return  
  end

  disp('WARNING: kthnorm and derived quantities only exact for s-alpha')

  out.krho_min = in.krhomax/max(in.nmod-1,1);
  out.kthnorm  = in.q/(2*pi*in.eps);
  out.ky_max   = in.krhomax/out.kthnorm;
  out.ky_min   = out.krho_min/out.kthnorm;
  out.kx_min   = 2*pi*in.shat*out.ky_min/in.ikxspace;
  out.kx_max   = out.kx_min*in.nx/2;
  out.lyn      = 2*pi/out.ky_min/out.kthnorm;
  out.lx       = 2*pi/out.kx_min;
  % Obsolete and deprecated
  % out.fft_y    = 2^floor(log2(1.5*(2*in.nmod-2))+1);
  % out.fft_x    = 2^floor(log2(1.5*(in.nx+1))+1);

end
