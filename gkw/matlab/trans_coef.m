% Determine the transport coefficients
% Usage in = trans_coef(in)
% where 'in' is the structure that contains the input variables as well as the fluxes
%
% The coefficients are added to the structure
%
% Dimensional transport coefficents are also calculated if the input 
% structure contains Bref, amin, COLLSIONS.[tref,rref,nref], geom.Jacobian
% and geom.g_eps_eps:  use with gkw_read_all
%
% WARNING:  The [q/v]fluxes are rescaled from the input by the T_R/v_R factors
% Therefore this routine should not be called more than once on the same structure
%
% WARNING FOR MILLER / CHEASE: The g_eps_eps factors are only added in the dimensional quantities
%

function in = trans_coef(in)

% Check if at least one of the electro-static fluxes is present
if ~isfield(in,'pfluxes')
    disp('ERROR: no fluxes in the input structure');
end;

% Check if begin and end point are specified
if ~isfield(in,'nstart')
    disp('You have not specified a start value of the average');
    disp('To do so set structure_name.nstart = value.');
    disp('For the time being the routine will use as starting');
    disp('half of the total interval');
    ntot = size(in.pfluxes,1);
    nstart = ntot / 2;
else
    nstart = in.nstart;
end;
if ~isfield(in,'nend')
    disp('You have not specified an end point for the average');
    disp('To do so set structure_name.nend = value');
    disp('For the time being the routine will use the end point');
    disp('of the run');
    nend = size(in.pfluxes,1);
else
    nend = in.nend;
end

in.nstart=nstart;
in.nend=nend;

nspec = size(in.pfluxes,2);


% Effective particle diffusivity
for i = 1: nspec
    [in.coef.Des(i) in.coef.Deser(i) dum] = average_error(in.pfluxes(:,i),nstart,nend);
    in.coef.Des(i) = in.coef.Des(i) / in.SPECIES(i).rln;
    in.coef.Deser(i) = in.coef.Deser(i) / in.SPECIES(i).rln;
    if isfield(in.CONTROL,'nlapar');
        if (in.CONTROL.nlapar(1)=='t')
            [in.coef.Dem(i) in.coef.Demer(i) dum] = average_error(in.pfluxem(:,i),nstart,nend);
            in.coef.Dem(i) = in.coef.Dem(i) / in.SPECIES(i).rln;
            in.coef.Demer(i) = in.coef.Demer(i) / in.SPECIES(i).rln;
        end
    end
    if isfield(in.CONTROL,'nlbpar');
        if (in.CONTROL.nlbpar(1)=='t')
            [in.coef.Dbp(i) in.coef.Dbper(i) dum] = average_error(in.pfluxbpar(:,i),nstart,nend);
            in.coef.Dbp(i) = in.coef.Dbp(i) / in.SPECIES(i).rln;
            in.coef.Dbper(i) = in.coef.Dbper(i) / in.SPECIES(i).rln;
        end
    end
    
     in.gamT_ov_q(i) = mean(in.pfluxes(nstart:nend,i),1)/mean(in.qfluxes(nstart:nend,i),1)
end;

 

% Effective heat conduction coefficient
for i = 1: nspec
    % First rescale heat fluxes by TR factors
    in.qfluxes(:,i)=in.SPECIES(i).temp*in.qfluxes(:,i);    
    [in.coef.Ces(i) in.coef.Ceser(i) in.coef.Cesstd(i)] = average_error(in.qfluxes(:,i),nstart,nend);
    in.coef.Ces(i) = in.coef.Ces(i) / in.SPECIES(i).rlt;
    in.coef.Ceser(i) = in.coef.Ceser(i) / in.SPECIES(i).rlt;
    in.coef.Cesstd(i) = in.coef.Cesstd(i) / in.SPECIES(i).rlt;
    
    if isfield(in.CONTROL,'nlapar');
        if (in.CONTROL.nlapar(1)=='t')
            in.qfluxem(:,i)=in.SPECIES(i).temp*in.qfluxem(:,i);
            [in.coef.Cem(i) in.coef.Cemer(i) in.coef.Cemstd(i)] = average_error(in.qfluxem(:,i),nstart,nend);
            in.coef.Cem(i) = in.coef.Cem(i) / in.SPECIES(i).rlt;
            in.coef.Cemer(i) = in.coef.Cemer(i) / in.SPECIES(i).rlt;
            in.coef.Cemstd(i) = in.coef.Cemstd(i) / in.SPECIES(i).rlt;
        end
    end
    if isfield(in.CONTROL,'nlbpar');
        if (in.CONTROL.nlbpar(1)=='t')
            in.qfluxbpar(:,i)=in.SPECIES(i).temp*in.qfluxbpar(:,i);
            [in.coef.Cbp(i) in.coef.Cbper(i) in.coef.Cbpstd(i)] = average_error(in.qfluxbpar(:,i),nstart,nend);
            in.coef.Cbp(i) = in.coef.Cbp(i) / in.SPECIES(i).rlt;
            in.coef.Cbper(i) = in.coef.Cbper(i) / in.SPECIES(i).rlt;
            in.coef.Cbpstd(i) = in.coef.Cbpstd(i) / in.SPECIES(i).rlt;
        end
    end
end;

% Effective momentum transport coefficient
for i = 1: nspec
    % First rescale heat fluxes by vR factors
    in.vfluxes(:,i)=sqrt(in.SPECIES(i).temp/in.SPECIES(i).mass)*in.vfluxes(:,i); 
    [in.coef.Ves(i) in.coef.Veser(i) ] = average_error(in.vfluxes(:,i),nstart,nend);
    in.coef.Ves(i) = in.coef.Ves(i) / in.SPECIES(i).uprim;
    in.coef.Veser(i) = in.coef.Veser(i) / in.SPECIES(i).uprim;
    if isfield(in.CONTROL,'nlapar');
        if (in.CONTROL.nlapar(1)=='t')
            in.vfluxem(:,i)=sqrt(in.SPECIES(i).temp/in.SPECIES(i).mass)*in.vfluxem(:,i); 
            [in.coef.Vem(i) in.coef.Vemer(i) dum] = average_error(in.vfluxem(:,i),nstart,nend);
            in.coef.Vem(i) = in.coef.Vem(i) / in.SPECIES(i).uprim;
            in.coef.Vemer(i) = in.coef.Vemer(i) / in.SPECIES(i).uprim;
        end
    end
    if isfield(in.CONTROL,'nlbpar');
        if (in.CONTROL.nlbpar(1)=='t')
            in.vfluxvp(:,i)=sqrt(in.SPECIES(i).temp/in.SPECIES(i).mass)*in.vfluxbp(:,i); 
            [in.coef.Vbp(i) in.coef.Vbper(i) dum] = average_error(in.vfluxbpar(:,i),nstart,nend);
            in.coef.Vbp(i) = in.coef.Vbp(i) / in.SPECIES(i).uprim;
            in.coef.Vbper(i) = in.coef.Vbper(i) / in.SPECIES(i).uprim;
            
        end
    end
end;

% If the rotation is nonzero also calculate the effective pinch
if isfield(in.ROTATION,'vcor')
    if (in.ROTATION.vcor ~= 0)
        for i = 1: nspec
            [in.coef.Vpes(i) in.coef.Vpeser(i) dum] = average_error(in.vfluxes(:,i),nstart,nend);
            in.coef.Vpes(i) = in.coef.Vpes(i) / in.ROTATION.vcor;
            in.coef.Vpeser(i) = in.coef.Vpeser(i) / in.ROTATION.vcor;
            if isfield(in.CONTROL,'nlapar');
                if (in.CONTROL.nlapar(1)=='t')
                    [in.coef.Vpem(i) in.coef.Vpemer(i) dum] = average_error(in.vfluxem(:,i),nstart,nend);
                    in.coef.Vpem(i) = in.coef.Vpem(i) / in.ROTATION.vcor;
                    in.coef.Vpemer(i) = in.coef.Vpemer(i) / in.ROTATION.vcor;
                end
            end
            if isfield(in.CONTROL,'nlbpar');
                if (in.CONTROL.nlbpar(1)=='t')
                    [in.coef.Vpbp(i) in.coef.Vpbper(i) dum] = average_error(in.vfluxbpar(:,i),nstart,nend);
                    in.coef.Vpbp(i) = in.coef.Vpbp(i) / in.ROTATION.vcor;
                    in.coef.Vpbper(i) = in.coef.Vpbper(i) / in.ROTATION.vcor;
                end
            end
        end;
    end
end

% Convert diffusivities to GB units if amin is present
% WARNING: Gyro/TGLF uses Bunit/=Bref so there is ambiguity in the GB units for elongated plasmas
% In that case one also must convert the B in rho_ref (this is NOT done here)
if isfield(in,'amin');
    
    % Assumes Species 1 is ions with T=Tref,
    % Assumes Species 2 is electrons
    
    % TO DO convert fluxes with T_R factors
    rova = in.COLLISIONS.rref/in.amin;
    in.rova=rova;
    
    in.time_gb=in.time*rova*sqrt(in.SPECIES(2).temp/2);       
    
    in.coef_gb.Ces   = in.coef.Ces(:)*((2/in.SPECIES(2).temp)^(1.5))/rova;
    in.coef_gb.Ceser = in.coef.Ceser(:)*((2/in.SPECIES(2).temp)^(1.5))/rova;
    in.coef_gb.Cesstd= in.coef.Cesstd(:)*((2/in.SPECIES(2).temp)^(1.5))/rova;
    
    in.coef_gb.Cem   = in.coef.Cem(:)*((2/in.SPECIES(2).temp)^(1.5))/rova;
    in.coef_gb.Cemer = in.coef.Cemer(:)*((2/in.SPECIES(2).temp)^(1.5))/rova;
    in.coef_gb.Cemstd= in.coef.Cemstd(:)*((2/in.SPECIES(2).temp)^(1.5))/rova;
    
%    in.coef_gb.Cbp   = in.coef.Cbp(:)*((2/in.SPECIES(2).temp)^(1.5))/rova;
%    in.coef_gb.Cbper = in.coef.Cbper(:)*((2/in.SPECIES(2).temp)^(1.5))/rova;
%    in.coef_gb.Cbpstd= in.coef.Cbpstd(:)*((2/in.SPECIES(2).temp)^(1.5))/rova;
    
%     [in.coef_gb.Cem(:)  in.coef_gb.Cemer(:) in.coef_gb.Cemstd(:)]=...
%         [in.coef.Cem(:) in.coef.Cemer(:) in.coef.Cemstd(:)]...
%         *((2/in.SPECIES(2).temp)^(1.5))/rova;
%     
%     [in.coef_gb.Cbp(:)  in.coef_gb.Cbper(:) in.coef_gb.Cbpstd(:)]=...
%         [in.coef.Cbp(:) in.coef.Cbper(:) in.coef.Cbpstd(:)]...
%         *((2/in.SPECIES(2).temp)^(1.5))/rova;
    
end


% Convert diffusivity to SI units if Bref is present
% For total power, also need the geometry Jacobian to obtain dVdr
if isfield(in,'Bref');
    
    % assumes deuterium as mref species
    in.mref=2*1.67e-27;
    in.vthref=sqrt(in.COLLISIONS.tref)*sqrt(2*1e3*1.6e-19/in.mref);
    in.cs=in.vthref*sqrt(in.SPECIES(2).temp/2)    
    
    in.rhoref=in.mref*in.vthref/(1.6e-19*in.Bref);
    in.rhostar=in.rhoref/in.COLLISIONS.rref;
    %strictly this should appear in the dimensionless chis and chi_GB also
    in.chi_gkw = in.rhoref^2*in.vthref/in.COLLISIONS.rref/mean(in.geom.g_eps_eps);    
    in.chi_gkw2 = in.rhostar^2*in.vthref*in.COLLISIONS.rref/mean(in.geom.g_eps_eps);
    
    if isfield(in,'amin') in.chi_gb = in.chi_gkw*((2/in.SPECIES(2).temp)^(1.5))/rova; end
    
    in.coef_si.Ces   = in.coef.Ces(:)*in.chi_gkw;
    in.coef_si.Ceser = in.coef.Ceser(:)*in.chi_gkw;
    in.coef_si.Cesstd= in.coef.Cesstd(:)*in.chi_gkw;
    
    in.coef_si.Cem   = in.coef.Cem(:)*in.chi_gkw;
    in.coef_si.Cemer = in.coef.Cemer(:)*in.chi_gkw;
    in.coef_si.Cemstd= in.coef.Cemstd(:)*in.chi_gkw;
    
    %in.coef_si.Cbp   = in.coef.Cbp(:)*in.chi_gkw;
    %in.coef_si.Cbper = in.coef.Cbper(:)*in.chi_gkw;
    %in.coef_si.Cbpstd= in.coef.Cbpstd(:)*in.chi_gkw;
    
    if isfield(in.geom,'jacobian');
       for isp = 1: nspec 
        
          % convert from diffusivities using ASTRA dVdr
%           in.coef_si.Qes_MW3(isp) = in.coef_si.Ces(isp)*mean(in.geom.g_eps_eps)/in.COLLISIONS.rref^2 ...                                  
%                                   *in.GEOM_PLUS.dvdr*in.COLLISIONS.rref ...
%                                   *in.COLLISIONS.tref ... %TR is done above in fluxes
%                                   *in.COLLISIONS.nref*in.SPECIES(isp).dens...
%                                   *1.60217657e-19*1e3*1e13...
%                                   *in.SPECIES(isp).rlt;
                                  
          % convert from diffusivities: P = dVdeps n T <g^eps eps> chi R/LT                    
          in.coef_si.Qes_MW(isp) = in.coef_si.Ces(isp)*mean(in.geom.g_eps_eps)/in.COLLISIONS.rref^2 ...                                  
                                  *in.geom.jacobian*in.COLLISIONS.rref^3 ...
                                  *in.COLLISIONS.tref...  %TR is done above in fluxes
                                  *in.COLLISIONS.nref.*in.SPECIES(isp).dens...
                                  *1.60217657e-19*1e3*1e13...
                                  *in.SPECIES(isp).rlt;
                              
          % convert from diffusivities: P = dVdeps n T <g^eps eps> chi R/LT                    
          in.coef_si.Qem_MW(isp) = in.coef_si.Cem(isp)*mean(in.geom.g_eps_eps)/in.COLLISIONS.rref^2 ...                                  
                                  *in.geom.jacobian*in.COLLISIONS.rref^3 ...
                                  *in.COLLISIONS.tref...  %TR is done above in fluxes
                                  *in.COLLISIONS.nref.*in.SPECIES(isp).dens...
                                  *1.60217657e-19*1e3*1e13...
                                  *in.SPECIES(isp).rlt;
                              
          % convert from diffusivities: P = dVdeps n T <g^eps eps> chi R/LT                    
          in.coef_si.Qeser_MW(isp) = in.coef_si.Ceser(isp)*mean(in.geom.g_eps_eps)/in.COLLISIONS.rref^2 ...                                  
                                  *in.geom.jacobian*in.COLLISIONS.rref^3 ...
                                  *in.COLLISIONS.tref...  %TR is done above in fluxes
                                  *in.COLLISIONS.nref.*in.SPECIES(isp).dens...
                                  *1.60217657e-19*1e3*1e13...
                                  *in.SPECIES(isp).rlt;
          
          % convert from fluxes directly  P = <Q.nabla eps> dVdeps                  
          in.coef_si.Qes_MW2(isp) = mean(in.qfluxes(nstart:nend,isp))/in.COLLISIONS.rref...
                                  *in.geom.jacobian*in.COLLISIONS.rref^3 ...
                                  *in.COLLISIONS.tref...  %TR is done above in fluxes
                                  *in.rhostar^2*in.vthref...
                                  *in.COLLISIONS.nref*in.SPECIES(isp).dens...
                                  *1.60217657e-16... % convert Kev to joules
                                  *1e13 % convert density and MW
                              
          in.fluxMW(isp) = in.coef_si.Qes_MW(isp) / mean(in.qfluxes(nstart:nend,isp));                     
                              
          % convert from fluxes directly  P = <Q.nabla eps> dVdeps                  
%           in.coef_si.Qes_MW4(isp) = mean(in.qfluxes(nstart:nend,isp))/in.COLLISIONS.rref...
%                                   *in.GEOM_PLUS.dvdr*in.COLLISIONS.rref ...
%                                   *in.COLLISIONS.tref.*in.SPECIES(isp).temp...         
%                                   *in.rhostar^2*in.vthref...
%                                   *in.COLLISIONS.nref.*in.SPECIES(isp).dens...
%                                   *1.60217657e-16... % convert Kev to joules
%                                   *1e13 % convert density and MW
                                  
                                  
       end
    end    
        
%     [in.coef_si.Ces  in.coef_si.Ceser in.coef_si.Cesstd]=...
%         [in.coef.Ces(:) in.coef.Ceser(:) in.coef.Cesstd(:)] ...
%         *in.chi_gkw;
%     
%     [in.coef_si.Cem  in.coef_si.Cemer in.coef_si.Cemstd]=...
%         [in.coef.Cem(:) in.coef.Cemer(:) in.coef.Cemstd(:)]...
%         *in.chi_gkw;
%     
%     [in.coef_si.Cbp  in.coef_si.Cbper in.coef_si.Cbpstd]=...
%         [in.coef.Cbp(:) in.coef.Cbper(:) in.coef.Cbpstd(:)]...
%         *in.chi_gkw;
    
    in.time_si=in.time*(in.COLLISIONS.rref/in.vthref);
    
end


end

