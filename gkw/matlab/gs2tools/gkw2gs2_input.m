function [gkwin gs2in] = gkw2gs2_input(infile,proj,gs2_template);
%
% function [gkwin gs2in] = gkw2gs2_input(gkwinfile,gkwproj,gs2template);
%
% convert a gkw input file into a gs2 input file (using a gs2 template)
% USE WITH CARE: check the outputs produced with a diff against the template
% WARNING: The GS2 template file used used the (now removed) sqrt(T/m) option
%
% Requires some knowledge of the input files and NORMALISATIONS of both codes.
%
% use gkwpath.m for gkw 
% gs2path is hard coded (all files in same folder).
%
% Assume Ti=Tref and mi=mref in both codes
% Collisions set for two species only
% Assume species 1 is ions
% Assume species 2 is electroms
%
% TO DO: make less case sensitive in read and write (always lowercase)
% FJC 17.07.2012
%

gs2_pthin = '~/gs2run/input/';
[gs2in.eql gs2in.spc gs2in.prkn gs2in.aky] = read_gs2input(gs2_template, gs2_pthin);

% Alternative under development: using gkw read / write - more general
% but GS2 structure names don't match up with namelist names
% gs2in = read_gkwinput(gs2_template,'template',2);

gkwin = read_gkwinput(infile,proj);

if gkwin.CONTROL.nlapar == 1 ;
    gs2in.prkn.fapar = 1; 
else
    gs2in.prkn.fapar = 0;
end;

if ~isfield(gkwin.GEOM,'geom_type')
  gkwin.GEOM.geom_type='s-alpha';
end

if isequal(gkwin.GEOM.geom_type,'s-alpha')
  kyfac=1;  
else
  gkwgeom = read_geom(infile,proj);    
  kthnorm_gs2=mean(gkwgeom.e_eps_zeta)*2;
  % renormalise kthnorm between GKW / GS2 (Miller /circ only ?)
  kyfac=kthnorm_gs2/gkwgeom.kthnorm;
  %disp(kyfac)
  %disp(gkwgeom.kthnorm)
end

nbspecies=gkwin.GRIDSIZE.number_of_species; 
%nbspecies=2;

%WARNING This is only correct if your GS2 setup uses the sqrt(T/m) option (no longer supported)
gs2in.aky.aky_min = gkwin.MODE.kthrho*kyfac/sqrt(2);

% First attempt to rescale ky
% kthnorm_sa=gkwin.GEOM.q/ (2*pi*gkwin.GEOM.eps);
% gs2in.aky.aky_min=gs2in.aky.aky_min*kthnorm_sa/gkwgeom.kthnorm/gkwin.GEOM.kappa;

gs2in.eql.s_hat_input = gkwin.GEOM.shat;
gs2in.eql.qinp = gkwin.GEOM.q;
gs2in.eql.pk = 2/gkwin.GEOM.q;
gs2in.eql.shat = gkwin.GEOM.shat;
gs2in.eql.epsl = 2.0; % makes a=R

%Which of these is defined changes the GS2 output to salpha / miller
gs2in.eql.rhoc = gkwin.GEOM.eps;
gs2in.eql.eps = gkwin.GEOM.eps;

if (isequal(gkwin.GEOM.geom_type,'circ')||isequal(gkwin.GEOM.geom_type,'miller'))
    gs2in.eql=rmfield(gs2in.eql,'eps');
elseif (isequal(gkwin.GEOM.geom_type,'s-alpha'))
    gs2in.eql=rmfield(gs2in.eql,'rhoc');
elseif(isequal(gkwin.GEOM.geom_type,'chease'))
    gs2in.eql=rmfield(gs2in.eql,'eps');
else    
    error('Cannot do this geom')
end

if (isfield(gkwin.GEOM,'beta_ref'))
  gs2in.prkn.beta = gkwin.SPCGENERAL.beta_ref;
elseif(isfield(gkwin.GEOM,'beta'))
  gs2in.prkn.beta = gkwin.SPCGENERAL.beta;
else
  gs2in.prkn.beta = 0.;
end

if isfield(gkwin.GEOM,'betaprime_ref')
  gs2in.xpareq.beta_prime_input = gkwin.SPCGENERAL.betaprime_ref;
  % note also BETAPRIME_TYPE='ref' in GKW input (default option)
else  
  gs2in.xpareq.beta_prime_input = 0.0;
end 
 
if ~isfield(gkwin,'ROTATION')
  gkwin.ROTATION.vcor=0;
end

if ~isfield(gkwin,'COLLISIONS')
  gkwin.COLLISIONS.coll_freq=0;
  gkwin.COLLISIONS.zeff=0;
end

for is=1:nbspecies
%gkwin.SPECIES(is).mass
gs2in.spc(is).mass = gkwin.SPECIES(is).mass;
gs2in.spc(is).z = gkwin.SPECIES(is).z;
gs2in.spc(is).temp = gkwin.SPECIES(is).temp;
gs2in.spc(is).tprim=gkwin.SPECIES(is).rlt;
gs2in.spc(is).fprim=gkwin.SPECIES(is).rln;
gs2in.spc(is).uprim=gkwin.SPECIES(is).uprim*sqrt(gs2in.spc(is).mass/gs2in.spc(is).temp);
gs2in.spc(is).upara=gkwin.ROTATION.vcor*sqrt(gs2in.spc(is).mass/gs2in.spc(is).temp);

%All species should have same uprim in GKW (normalised to vthref), not in CLA GS2 (normalised to vsp)

end

%This may also assume Ti=Tref
%gkwin.ROTATION.VCOR = gs2in.spc(1).upara;

%Assume Tref=Ti in both GKW and GS2

if (gkwin.CONTROL.collisions == 1)
    
    if (gkwin.COLLISIONS.freq_override == 1)
        %ion-ion collision freq in GKW, equivalent to GS2 vnewki
        %Assuming Te=Tref in both codes
        gs2in.spc(2).vnewk=gkwin.COLLISIONS.coll_freq/sqrt(gkwin.SPECIES(2).mass);
        
        %Assuming Ti=Tref in both codes
        gs2in.spc(2).vnewk=gkwin.COLLISIONS.coll_freq/(sqrt(gkwin.SPECIES(2).mass)*(gs2in.spc(2).temp)^(1.5));
        
    else
        
        Te=gs2in.spc(2).temp*gkwin.COLLISIONS.tref;
        ne=gkwin.COLLISIONS.nref;
        coulog = 15.94 - 0.5*log(ne/Te^2);
        %spc(2).vnewk = 0.00279 *coulog*ne/Te^1.5*Rtor*(2/Tref)^0.5;
        gs2in.spc(2).vnewk=0.00279*coulog*ne*gkwin.COLLISIONS.rref*sqrt((2/gkwin.COLLISIONS.tref))/(Te^1.5);
        
        %disp('WARNING: not checked the calculation of GS2 collisionality from nref etc')
        
        %return;
        
    end    
    
    %see insert collisionality
    % assume mref=spc(1).mass, 2 species    
    gs2in.spc(1).vnewk = gs2in.spc(2).vnewk * ...
     sqrt(gs2in.spc(2).mass/gs2in.spc(1).mass) *(gs2in.spc(2).temp/gs2in.spc(1).temp)^1.5;
else
    % Can't turn off collision model in GS2 without rewriting make_cgs2input
    gs2in.spc(1).vnewk=0.0;    
    gs2in.spc(2).vnewk=0.0;
end

gs2in.prkn.zeff=gkwin.COLLISIONS.zeff;

if (isequal(gkwin.GEOM.geom_type,'miller'))
  [gs2in.eql gkwin] = gs2_gkw_gyro_miller(gs2in.eql,gkwin,'gkw');
end

tmp=make_cgs2input(infile,gs2in.eql,gs2in.aky,nbspecies,gs2in.spc,gs2in.prkn,gs2_pthin, gs2_template);
 
% write_gkwinput(gs2in,infile,'gkw2gs2',0,'made with gkw2gs2input.m');

% Write GS2 geometry input filename for EQDSK
% ind=find(infile(:)=='r');
% unix(['perl -p -i -n -e ''s$eqdsk_filename$../../eqdsk/' infile(1:ind-2) '$g'' ' [gs2_pthin infile]]);

end


