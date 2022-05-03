function [gkwin gs2in] = gs2gkw_input(gs2_file,gkw_template,gkw_file,gkw_proj);
%
% function [gkwin gs2in] = gs2gkw_input(gs2_file,gkw_template,gkw_file,gkw_proj);
%
% Convert a GS2 input file from GS2GET to a GKW input file using a template for GKW input
%
% e.g. give in input gs2_file = 'xAUG25832_6.775_k_3_CxUB_tr_0_up1';
%
% Only for miller geometry and irho=2
% Exactly two species only 
%  Species 1 is ions (assumed deuterium)
%  Species 2 is electrons
%
% Use with care - check the output manually !
% 
% FJC 17.12.12

gs2_proj=gkw_proj;

% Read the GS2 file from GKW input folder using the generic GKW namelist read
gs2in=read_gkwinput(gs2_file,gs2_proj);
if (gs2in.theta_grid_eik_knobs.irho ~= 2) 
   error('I only do irho=2 for gs2 radial coordanate')    
end

if (gs2in.theta_grid_eik_knobs.bishop ~= 4) 
   error('I only work with bishop = 4 in gs2')    
end

%convert eqm to cla format
eql=gs2in.theta_grid_parameters;
eql.beta_prime_input=gs2in.theta_grid_eik_knobs.beta_prime_input;
eql.s_hat_input=gs2in.theta_grid_eik_knobs.s_hat_input;

gkwin = read_gkwinput(gkw_template,gkw_proj);

% if gs2in.knobs.fapar == 1;
%     gkwin.CONTROL.nlapar = '.true.'; 
% end;

gkwin.MODE.kthrho = gs2in.kt_grids_range_parameters.aky_min;
gkwin.MODE.krhomax = gs2in.kt_grids_range_parameters.aky_max;
gkwin.COLLISIONS.zeff=gs2in.parameters.zeff;
gkwin.SPCGENERAL.beta_ref = gs2in.parameters.beta;

gkwin.GEOM.shat= eql.s_hat_input;
gkwin.GEOM.q =  eql.qinp;
gkwin.GEOM.eps = eql.rhoc/eql.rmaj;
gkwin.GEOM.geom_type= 'miller';

[eqlout gkwin dum] = gs2_gkw_gyro_miller(eql,gkwin,'gs2')

%get Nref, Lref, Tref, Omega, and grad-Omega from GS2GET params
[idum nref] = unix(['grep \#nref ' [gkwpath('input',gs2_proj) gs2_file] ' | cut -c 10-21 ']);
[idum tref] = unix(['grep \#tref ' [gkwpath('input',gs2_proj) gs2_file] ' | cut -c 10-21 ']);
[idum lref] = unix(['grep \#lref ' [gkwpath('input',gs2_proj) gs2_file] ' | cut -c 10-21 ']);

nref=str2num(nref);
tref=str2num(tref);
lref=str2num(lref);

gkwin.COLLISIONS.rref=lref*eql.rmaj;
gkwin.COLLISIONS.tref=tref/1e3;
gkwin.COLLISIONS.nref=nref/1e19;

[idum omega] = unix(['grep \#omega ' [gkwpath('input',gs2_proj) gs2_file] ' | cut -c 10-21 ']);
[idum domegadr] = unix(['grep \#domegadr ' [gkwpath('input',gs2_proj) gs2_file] ' | cut -c 10-21 ']);
omega=str2num(omega);
domegadr=str2num(domegadr);
vthref=sqrt(2*tref*1e-19/(2*1.66e-27));
gkwin.ROTATION.vcor=omega*gkwin.COLLISIONS.rref/vthref;

%[idum lref] = unix(['grep l ' [gkwpath('input',proj) gs2file] ' | cut -c 10-21 '])

gkwin.SPECIES(1).mass= gs2in.species_parameters_1.mass;
gkwin.SPECIES(1).z= gs2in.species_parameters_1.z;
gkwin.SPECIES(1).temp= gs2in.species_parameters_1.temp;
gkwin.SPECIES(1).dens= 1.0;
gkwin.SPECIES(1).rlt= gs2in.species_parameters_1.tprim*eql.rmaj;
gkwin.SPECIES(1).rln= gs2in.species_parameters_2.fprim*eql.rmaj;
%This assumes Ti=Tref in GS2
gkwin.SPECIES(1).uprim= -domegadr*eql.rmaj/vthref;

gkwin.SPECIES(2).mass= gs2in.species_parameters_2.mass;
gkwin.SPECIES(2).z= gs2in.species_parameters_2.z;
gkwin.SPECIES(2).temp= gs2in.species_parameters_2.temp;
gkwin.SPECIES(2).dens= 1.0;
gkwin.SPECIES(2).rlt= gs2in.species_parameters_2.tprim*eql.rmaj;
gkwin.SPECIES(2).rln= gs2in.species_parameters_2.fprim*eql.rmaj;
%All species should have same uprim in GKW, not in CLA GS2
gkwin.SPECIES(2).uprim= -domegadr*eql.rmaj/vthref;

write_gkwinput(gkwin, gkw_file, gkw_proj,1, ['Converted from gs2 file ' gs2_file ' in matlab using gs2gkw_miller'])
 
end


