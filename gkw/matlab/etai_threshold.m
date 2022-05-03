function [thres limit]=etai_threshold(varargin)
% Returns the etai stability threshold in various limits
%
% [thres limit]=etai_threshold(shat,q,teti,rln,rlt,krhoi,eps)
%
% [or] [thres limit]=etai_threshold(gkwinput_structure)
%
% [or] [thres limit]=etai_threshold(shat,q,gkwinput_structure)
%                    (overrides shat and q in input_structure, e.g. for chease)
%
% gkwinput_structure with sp. 1 as ions with temp=1.0 and sp. 2 as electrons
%
% Refs: Guo and Romanelli, POF B5 520 (1993), Wesson book page 437
%       Biglari, Diamond, Rosenbluth, POF B1 109 (1989)
%       Romanelli and Briguglio. POF B2 754 (1990)
%
% Would be nice to have a similar function for TEM threshold, but simple formulae 
% including the collisionality and Te/Ti dependance do not seem to exist. 
% (e.g. Peeters POP 12 022505 (2005) and Casati POP 15 042310 (2008)) 
% Weiland fluid model could be used as an approximation.
%
% Does not include any geometry effects (e.g particular elongation is stabilising
% but Shafranov shift is destabilising
%
% FJC: 17.10.2012


if (nargin==1)     % use GKW input structure
  in=cell2mat(varargin(1));
  
  shat    = in.GEOM.shat
  q       = in.GEOM.q
  rln     = in.SPECIES(1).rln;
  rlt     = in.SPECIES(1).rlt;
  teti    = in.SPECIES(2).temp;
  krhoi   = in.MODE.kthrho;   % normalisation in Guo is same as ours
  eps     = in.GEOM.eps;
  
elseif(nargin== 3)  %use GKW input structure but override q and shat
  shat    = cell2mat(varargin(1))
  q       = cell2mat(varargin(2))
  in=cell2mat(varargin(3));    
  rln     = in.SPECIES(1).rln;
  rlt     = in.SPECIES(1).rlt;
  teti    = in.SPECIES(2).temp;
  krhoi   = in.MODE.kthrho; % normalisation in Guo is same as ours
  eps     = in.GEOM.eps; 

elseif (nargin==7)  % manually input all paramters
  shat     = cell2mat(varargin(1));
  q        = cell2mat(varargin(2));
  teti     = cell2mat(varargin(3));
  rln      = cell2mat(varargin(4));
  rlt      = cell2mat(varargin(5));
  krhoi    = cell2mat(varargin(6));
else
  thres=NaN;
  disp('etai_threshold: incorrect number of input arguments')
  return  
end

% notation used in the papers
tau=teti;
en = 1/rln;
et = 1/rlt;
bs = krhoi^2/2;   % Could not find the definition in Guo or Romanelli; but it is in Biglari

% short wavelength, moderate shear limit
if (krhoi > et^(0.25))
    
    enc = 0.9 / ((1.0 + 1.0/tau)*(1.0+2.0*shat/q));  %Wesson
    %enc = 0.5 / (1 + 1/tau);                %Guo
    
    if (en < enc)
        thres = 1.2;
    elseif (en >=enc)
        thres = (4.0/3.0)*(1.0 + 1.0/tau)*(1.0+2.0*shat/q)*en;
    else
        thres=NaN;
    end
    
    limit = 'short wavelength toroidal';
    
% long wavelength limit, toroidal branch
elseif (krhoi < et^(0.25))
    
    % long wavelength limit
    thres = 1.0 + sqrt(1.0 + (en^2 / (q^2*bs))*(1.0+1.0/tau));
    
    limit = 'long wavelength toroidal';
    
    % trapped ion limit at flat density, low collisionality
    thres2 = (1.0 + 1.0/tau)*sqrt(2*eps)*en;
        
    %if (thres2 < thres)
    %    thres = thres2;
    %    limit = 'long wavelength trapped ion';
    %end
    
end




end