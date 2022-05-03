% This function returns the default path for GKW runs
% assuming file locations of script gkwnlin
% 
% Usage: gkwpath('data','project')
%
% If the environment variable GKW_RUNS_PATH is set, this will be used
% as the base run directory.  If not, it will be selected by username.
%
% The environment variable can be set from within matlab as follows
%    setenv('GKW_RUNS_PATH','/some/special/gkw/runs')
%
% TO DO:
% Make this function able to return normal gkw filenames if no project 
% is provided, for interoperability with runs in which the gkwnlin folder
% structure is not used and all gkw data is in one folder.

function [flpth]=gkwpath(str,proj)

global input_data

if ~exist('proj')||isempty(proj)
	proj='default';
end

%Get usr environment variables
usr=getenv('USER');
gkw_runs_path=getenv('GKW_RUNS_PATH');

%Setup path to runs_home specific to the user.
if isempty(gkw_runs_path)
 
    switch usr
        
        case('dzar')
            runs_home = eval('pwd');
            runs_home = [runs_home '/'];
        
        case('btpp01') % Runs in Bayreuth            
            runs_home=['/home/btpp/btpp01/runs/'];
            
        case{'bt301007'}
            runs_home=['/localdisk/bt301007/runs/gkw/'];
            
        case{'user2'} %set runs_home specific to user2
            runs_home=['/path/to/gkw/runs/'];

        otherwise % generic default,
            runs_home=['~/runs/gkw/'];            
            
      end %switch usr
    
else
    runs_home = gkw_runs_path;
    runs_home = [runs_home '/'];
end
      

%Return path to project folder
switch lower(str)
    case{'time'}
        flpth=[runs_home proj '/time/'];
    case{'fluxes'}
        flpth=[runs_home proj '/fluxes/'];
    case{'fluxes_em'}
        flpth=[runs_home proj '/fluxes/em/'];
    case{'fluxes_bpar'}
        flpth=[runs_home proj '/fluxes/bpar/'];
    case{'fluxes_neoc'}
        flpth=[runs_home proj '/fluxes/neoc/'];
    case{'fluxes_det'}
        flpth=[runs_home proj '/fluxes_det/'];
    case{'parallel'}
        flpth=[runs_home proj '/parallel/'];
    case{'distr1'}
        flpth=[runs_home proj '/distr1/'];
    case{'distr2'}
        flpth=[runs_home proj '/distr2/'];
    case{'distr3'}
        flpth=[runs_home proj '/distr3/'];
    case{'distr4'}
        flpth=[runs_home proj '/distr4/'];
    case{'kxrh'}
        flpth=[runs_home proj '/grids/kxrh/'];
    case{'krho'}
        flpth=[runs_home proj '/grids/krho/'];
    case{'input'}
        flpth=[runs_home proj '/input/'];
    case{'input_out'}
        flpth=[runs_home proj '/input_out/'];
    case{'perform'}
        flpth=[runs_home proj '/perform/'];
    case{'geom'}
        flpth=[runs_home proj '/geom/'];
    case{'hamada'}
        flpth=[runs_home proj '/hamada/'];
    case{'scan'}
        flpth=[runs_home proj '/input/scan/'];
    case{'root'}
        flpth=[runs_home proj '/'];
    case{'radial'}
        flpth=[runs_home proj '/radial/'];
    case{'prof_back'}
        flpth=[runs_home proj '/radial/prof_back/'];
    case{'top'}
        flpth=[runs_home];
    case{'spectrum'}
        flpth=[runs_home proj '/spectrum/'];
    case{'phi'}
        flpth=[runs_home proj '/slices_2D/phi/'];
    case{'parallel_phi'}
        flpth=[runs_home proj '/parallel/phi_time/'];
    case{'other'}
        flpth=[runs_home proj '/other/'];
    case{'out'}
        flpth=[runs_home proj '/out/'];
    case{'den01'}
        flpth=[runs_home proj '/slices_2D/den01/'];
    case{'den02'}
        flpth=[runs_home proj '/slices_2D/den02/'];
    case{'den03'}
        flpth=[runs_home proj '/slices_2D/den03/'];
    case{'ene01'}
        flpth=[runs_home proj '/slices_2D/ene01/'];
    case{'ene02'}
        flpth=[runs_home proj '/slices_2D/ene02/'];
    case{'ene03'}
        flpth=[runs_home proj '/slices_2D/ene03/'];
    case{'apa'}
        flpth=[runs_home proj '/slices_2D/apa/'];
    case{'parallel_phi'}
        flpth=[runs_home proj '/parallel/phi_time/'];
    case{'xphi'}
        switch(usr)
            case('dzar')
                flpth=[runs_home proj '/other/'];
            otherwise
                flpth=[runs_home proj '/grids/xphi/'];
        end
    case{'yphi'}
          switch(usr)
            case('dzar')
                flpth=[runs_home proj '/other/'];
            otherwise
                flpth=[runs_home proj '/grids/yphi/'];
          end
    %All the following are probably not needed, if the method in spectra.m is used.
    case{'kyspec'}
        flpth=[runs_home proj '/spectrum/kyspec/'];
    case{'kxspec'}
        flpth=[runs_home proj '/spectrum/kxspec/'];
    case{'kyspec_em'}
        flpth=[runs_home proj '/spectrum/kyspec_em/'];
    case{'kxspec_em'}
        flpth=[runs_home proj '/spectrum/kxspec_em/'];
    case{'pflux_spec'}
        flpth=[runs_home proj '/spectrum/pflux/'];
    case{'eflux_spec'}
        flpth=[runs_home proj '/spectrum/eflux/'];
    case{'eflux_sup'}
        flpth=[runs_home proj '/spectrum/eflux_sup/'];
    case{'pflux_sup'}
        flpth=[runs_home proj '/spectrum/pflux_sup/'];
    case{'pflux_xspec'}
        flpth=[runs_home proj '/spectrum/pflux_kx/'];
    case{'eflux_xspec'}
        flpth=[runs_home proj '/spectrum/eflux_kx/'];
    case{'vflux_xspec'}
        flpth=[runs_home proj '/spectrum/vflux_kx/'];
    case{'cfdens'}
        flpth=[runs_home proj '/cfdens/'];
    otherwise
        error('gkwpath does not know this file type - please add it if required')
end
