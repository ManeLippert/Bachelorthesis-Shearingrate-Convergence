function [gs2mil gkwin gyroin]=gs2_gkw_gyro_miller(gs2mil,gkwin,from)
% Converts gkw or gs2 miller inputs and produces miller inputs for gyro, gs2 and gkw
%
% usage
% [gs2in gkwin gyroin] = gs2_gkw_gyro_miller(gs2mil,gkwin,from)
%
% where gs2mil is the structure containing gs2 miller params
% gkwin is the structure containing the full gkw input file
% from is either 'gs2' or 'gkw'
% when going from gs2 a gkw input should be provided with GEOM.eps already set
%
% Notes 
%  assumes gs2 with irho=2 
%  which implies gs2.rhoc = gkw.GEOM.eps*gs2.Rmaj
% 
% definitions From https://web.gat.com/theory/Gyroinput and GS2 reference input
% http://sourceforge.net/apps/mediawiki/gyrokinetics/index.php?title=GS2_Reference_Input_File
% GS2 needs to be set with bishop = 4 to use beta_prime input (or bishop = 5 to use alpha input ?)
%
% In GKW gradp_type=alpha and alpha_mhd differ when kappa/=1.  It is not yet clear which (if either)
% corresponds to GS2 alpha_input with bishop = 5.
%
% Gyro and GKW are the same, with the exception of the sdelta parameter 
% and rescalings of the radial coordinate.  Note that gyro does not have a beta_prime_input
% so for comparison in this case gkw should run with beta_prime_type='sp'. 
%
% FJC 16.07.2012

if (from == 'gs2')
    gkwin.GEOM.kappa = gs2mil.akappa;
    gkwin.GEOM.skappa = gs2mil.rmaj*gs2mil.akappri/gkwin.GEOM.kappa*gkwin.GEOM.eps;
    gkwin.GEOM.delta = sin(gs2mil.tri);
    %gkwin.GEOM.sdelta = gs2mil.rmaj*gs2mil.tripri*cos(gs2mil.tri)*gkwin.GEOM.eps/sqrt(1-gkwin.GEOM.delta^2);
    % should be equivalent to above, but a lot  simpler
    gkwin.GEOM.sdelta = gs2mil.rhoc * gs2mil.tripri;
    gkwin.GEOM.square = 0.;
    gkwin.GEOM.ssquare = 0.;
    gkwin.GEOM.dRmil = gs2mil.shift;   % GKW definition should be clarified for this one - going to assume dimensionless
    gkwin.GEOM.Zmil = 0.;
    gkwin.GEOM.dZmil = 0.;
    gkwin.GEOM.geom_type='miller';
    gkwin.GEOM.gradp_type='beta_prime';
    if isfield(gs2mil,'beta_prime_input')
      gkwin.SPCGENERAL.betaprime_type='ref';
      gkwin.SPCGENERAL.betaprime_ref=gs2mil.rmaj*gs2mil.beta_prime_input;
    end
elseif (from == 'gkw') % inversion
    gs2mil.akappa = gkwin.GEOM.kappa;
    gs2mil.akappri= gkwin.GEOM.skappa*gkwin.GEOM.kappa/gkwin.GEOM.eps;
    gs2mil.tri = asin(gkwin.GEOM.delta);
    gs2mil.tripri=gkwin.GEOM.sdelta/cos(gs2mil.tri)/gkwin.GEOM.eps*sqrt(1-gkwin.GEOM.delta^2);
    gs2mil.shift = gkwin.GEOM.dRmil;   
    gs2mil.beta_prime_input=gkwin.SPCGENERAL.betaprime_ref;  
    if (isequal(gkwin.GEOM.gradp_type,'alpha_mhd'))
      gs2mil.alpha_input=gkwin.GEOM.gradp; 
    end
else
    error('I can only go from gs2 or gkw at present')
end


% GKW Miller Not yet benchmarked against GYRO !
gyroin.aspect_ratio = 3;
gyroin.kappa = gkwin.GEOM.kappa;
gyroin.s_kappa = gkwin.GEOM.skappa;
gyroin.delta = gkwin.GEOM.delta;
gyroin.s_delta = gkwin.GEOM.sdelta*sqrt(1-gkwin.GEOM.delta^2);
gyroin.shift = gkwin.GEOM.dRmil;
gyroin.zeta = gkwin.GEOM.square;
gyroin.s_zeta = gkwin.GEOM.ssquare;
gyroin.z_mag = gkwin.GEOM.Zmil*gyroin.aspect_ratio;
gyroin.dzmag = gkwin.GEOM.dZmil;     % The gkw definition needs to be clarified - going to assume dimensionless
gyroin.geo_beta_prime_scale=1.0; % uses the beta prime in geom calculated from profiles 
% which should be equivalent to gkw.SPCGENERAL.beta_prime_type='sp'

end 
