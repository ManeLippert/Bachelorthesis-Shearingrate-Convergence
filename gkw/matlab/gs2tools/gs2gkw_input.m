function [gkwin gs2in] = gs2gkw_input(gs2_file,gkw_template,gkw_proj);
%
% function [gkwin gs2in] = gs2gkw_input(gs2_file,gkw_template,gkw_proj);
%
% Convert a GS2 input file to a GKW input file using a template for GKW input
%
% e.g. give in input gs2_file = 'xAUG25832_6.775_k_3_CxUB_tr_0_up1';
%
% Only for simple geometries - other not included (they are in gkw2gs2_input)
% Exactly two species only (more are allowed in gkw2gs2_input)
%  Species 1 is ions (assumed deuterium)
%  Species 2 is electrons
%
% Use with care - check the output manually !
%
% GS2 file directory is './' (hard coded for now)
%
% FJC 01.10.11

gs2_pthin = './';
[gs2in.eql gs2in.spc gs2in.prkn gs2in.aky] = read_cgs2input(gs2_file, gs2_pthin);

gkwin = read_gkwinput(gkw_template,gkw_proj);

if gs2in.prkn.fapar == 1; gkwin.CONTROL.nlapar = '.true.'; end;

gkwin.MODE.kthrho = gs2in.aky.aky_min.*sqrt(2);

gkwin.GEOM.SHAT = gs2in.eql.s_hat_input;
gkwin.GEOM.Q = gs2in.eql.qinp;
gkwin.GEOM.EPS = gs2in.eql.rhoc;
gkwin.GEOM.GEOM_TYPE= 'circ';

gkwin.SPCGENERAL.beta = gs2in.prkn.beta;

gkwin.SPECIES(1).MASS= gs2in.spc(1).mass;
gkwin.SPECIES(1).Z= gs2in.spc(1).z;
gkwin.SPECIES(1).TEMP= gs2in.spc(1).temp;
gkwin.SPECIES(1).dens= 1.0;
gkwin.SPECIES(1).rlt= gs2in.spc(1).tprim;
gkwin.SPECIES(1).rln= gs2in.spc(2).fprim;
%This assumes Ti=Tref in GS2
gkwin.SPECIES(1).uprim= gs2in.spc(1).uprim;

gkwin.SPECIES(2).MASS= gs2in.spc(2).mass;
gkwin.SPECIES(2).Z= gs2in.spc(2).z;
gkwin.SPECIES(2).TEMP= gs2in.spc(2).temp;
gkwin.SPECIES(2).dens= 1.0;
gkwin.SPECIES(2).rlt= gs2in.spc(2).tprim;
gkwin.SPECIES(2).rln= gs2in.spc(2).fprim;
%All species should have same uprim in GKW, not in CLA GS2
%gkwin.SPECIES(2).uprim= gs2in.spc(2).uprim*gs2in.spc(2).temp;
gkwin.SPECIES(2).uprim= gs2in.spc(1).uprim;

% gkwin.SPECIES(3).MASS= gs2in.spc(3).mass;
% gkwin.SPECIES(3).Z= gs2in.spc(3).z;
% gkwin.SPECIES(3).TEMP= gs2in.spc(3).temp;
% gkwin.SPECIES(3).dens= 1.0;
% gkwin.SPECIES(3).rlt= gs2in.spc(3).tprim;
% gkwin.SPECIES(3).rln= gs2in.spc(3).fprim;
% gkwin.SPECIES(3).uprim= gs2in.spc(3).uprim;

gkwin.ROTATION.VCOR = gs2in.spc(1).upara;

gkwin.COLLISIONS.freq_override='.true.'
%ion-ion collision freq in GKW, otherwise codes use same def for coll_freq

%Assumes Tref=Te in both GKW and GS2
%gkwin.COLLISIONS.coll_freq = gs2in.spc(2).vnewk/sqrt(gs2in.spc(2).mass);

%Assumes Tref=Ti in both GKW and GS2
gkwin.COLLISIONS.coll_freq = gs2in.spc(2).vnewk*(gs2in.spc(2).temp)^(1.5)*sqrt(gs2in.spc(2).mass);

gkwin.COLLISIONS.zeff=gs2in.prkn.zeff;

end


