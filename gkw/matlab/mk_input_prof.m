%
% Small matlab example on how to generate an inputfile containing the profiles for GKW 
%
% The script is an example only of how the generated file input.prof should look. The 
% script in its present form does not generate are reasonable set of input parameters 
% for GKW! 
%
% This makes sure the scripts in the GKW matlab directory can be used. This line should 
% in fact be removed and one should make sure all GKW matlab scripts can be accessed 
% elsewhere (change this line for your own machine) 
% addpath /home/btpp/bt161087/GKW/matlab

% -------------------------------------------------------------------------------------
% The part below reads the GKW input file and calculates the grid points on which 
% the data needs to be given 
%--------------------------------------------------------------------------------------

% This line reads the input.dat file from the current directory 
[GKW_in,sss_in]=read_gkwinput('input.dat');

% Set the boundary values of psi and the number of grid points 
psil     = GKW_in.GRIDSIZE.psil; 
psih     = GKW_in.GRIDSIZE.psih; 
n_x_grid = GKW_in.GRIDSIZE.nx; 
radial_boundary_conditions = GKW_in.CONTROL.radial_boundary_conditions;

% Then calculate the values of the psi grid 
if ((psil <= 0) & (radial_boundary_conditions == 'Dirichlet')) 
  psil = psih / (2.E0*n_x_grid + 1.E0);
end 
for ix = 1: n_x_grid 
  xgr(ix) = psil + (psih-psil)*(ix-0.5E0)/n_x_grid;
end;


%------------------------------------------------------------------------------------
% Define some experimental data 
%------------------------------------------------------------------------------------

r = [0. 0.5 1.0]; 
% q and shat are separately interpolated. One could change this and calculate 
% shat from the derivative of q. Then only q would need to be specified. 
% similarly one could do the same for kappa, delta, square etc. 
% Note that the data is freely invented here. There is no consistency between 
% values and derivatives. 
q = [1. 1.3 3.0]; 
shat = [0. 0.7 1.5];  
Te = [1.5 1.2 0.2]; 
rlte = [6 7 6];
kappa = [1.0 1.0 1.0]; 
delta = [0.0 0.0 0.0]; 
square = [0. 0. 0.];
skappa = [0. 0. 0.]; 
sdelta = [0. 0. 0.]; 
ssquare = [0. 0. 0.]; 
zmil = [0. 0. 0.]; 
drmil = [0. 0. 0.]; 
dzmil = [0. 0. 0.]; 
gradp = [0. 0. 0.]; 
nions = 2; 
Z(1) = 1; Z(2) = 48;
me       = 9.1e-31; 
mions(1) = 1.6e-27; 
mions(2) = 1.6e-25; 
Ti(1,:) = [1 0.7 0.2];
Ti(2,:) = [1 0.7 0.2]; 
rlti(1,:) = [6 6 6];
rlti(2,:) = [3 3 3]; 
ne  = [1 0.9 0.7]*1e19;
ni(1,:) = [1 0.9 0.7]*1e19; 
ni(2,:) = [0 0 0];
rlne = [2 2 2]; 
rlni(1,:) = [2 2 2]; 
rlni(2,:) = [2 2 2]; 

%------------------------------------------------------------------------------------
% The radial coordinates at which this data is defined must be transformed into the 
% GKW radial coordinate (r / R_0) 
%------------------------------------------------------------------------------------
r = r / 3;  

%------------------------------------------------------------------------------------
% The values of all the data needs to be interpolated onto the GKW grid 
%------------------------------------------------------------------------------------
te_gkw = spline(r,Te,xgr); 
rlte_gkw = spline(r,rlte,xgr); 
for i = 1: nions 
  ti_gkw(i,:) = spline(r,Ti(i,:),xgr); 
  rlti_gkw(i,:) = spline(r,rlti(i,:),xgr);
end; 
ne_gkw = spline(r,ne,xgr); 
rlne_gkw = spline(r,rlne,xgr); 
for i = 1: nions 
  ni_gkw(i,:) = spline(r,ni(i,:),xgr);
  rlni_gkw(i,:) = spline(r,rlni(i,:),xgr);  
end; 
q_gkw = spline(r,q,xgr); 
shat_gkw    = spline(r,shat,xgr); 
kappa_gkw   = spline(r,kappa,xgr); 
delta_gkw   = spline(r,delta,xgr); 
square_gkw  = spline(r,square,xgr); 
skappa_gkw  = spline(r,skappa,xgr); 
sdelta_gkw  = spline(r,sdelta,xgr); 
ssquare_gkw = spline(r,ssquare,xgr); 
zmil_gkw    = spline(r,zmil,xgr); 
drmil_gkw   = spline(r,drmil,xgr); 
dzmil_gkw   = spline(r,dzmil,xgr); 
gradp_gkw   = spline(r,gradp,xgr); 

%------------------------------------------------------------------------------------
% Then the quantities need to be normalized. Here simply a point is chosen in the 
% grid of GKW and the values at that location are used to normalize the profiles. 
% Note that the velocity grids of the different species can be normalized with a 
% different temperature (tgrid) but that one temperature must be chosen for the 
% timestep and field normalization (tref). Below Tref is chosen to be equal to the 
% ion temperature. (but this must not be the case). The density is normalized to the 
% electron density at the same position (nref). Again a different normalization 
% density can be applied to all species. (Note that one must always satisfy quasi-
% neutrality)  
% The mass is normalized to a chosen value (the ion mass) 
% The code does not need to know the values of mref Tref and nref. But the values 
% of tegrid tigrid negrid and nigrid need to be specified since these are relative 
% normalizations. 
%------------------------------------------------------------------------------------
ix_norm_pos = round(n_x_grid / 2); 
Tref = ti_gkw(1,ix_norm_pos); 
disp 'The temperature is normalized to: ',Tref
tigrid(1) = 1.0; 
for i = 2: nions 
  tigrid(i) = ti_gkw(i,ix_norm_pos)/Tref; 
end; 
tegrid = te_gkw(ix_norm_pos)/Tref; 

% The actual normalization 
te_gkw = te_gkw / tegrid / Tref ; 
for i = 1: nions 
  ti_gkw(i,:) = ti_gkw(i,:) / tigrid(i) / Tref; 
end; 

nref = ne_gkw(ix_norm_pos); 
disp 'The density is normalized to: ',nref
negrid = 1.0; 
for i = 1: nions 
  nigrid(i) = ni_gkw(i,ix_norm_pos)/ nref; 
  if (nigrid(i) == 0.) nigrid(i) = 1.; end;  
end; 

% actual normalization 
ne_gkw = ne_gkw / nref; 
for i = 1: nions 
  ni_gkw(i,:) = ni_gkw(i,:) / nigrid(i) / nref; 
end; 

% Set the value of mref equal to the ion mass 
mref = 1.6e-27; 
me = me / mref; 
mions = mions / mref; 

% Some graphic output 
figure(1); 
hold off; 
plot(xgr,te_gkw,'r'); 
hold on 
plot(xgr,ti_gkw,'b'); 
figure(2); 
hold off; 
plot(xgr,ne_gkw,'r'); 
hold on; 
plot(xgr,ni_gkw,'b'); 


%----------------------------------------------------------
% Generate the output file 
%----------------------------------------------------------
fileID = fopen('input.prof','w');
fprintf(fileID,'%s\n','#Geometry: psi, q, shat'); 
for i = 1: n_x_grid 
  fprintf(fileID,'%12.8f %12.8f %12.8f %12.8f %12.8f %12.8f %12.8f %12.8f %12.8f %12.8f %12.8f %12.8f %12.8f\n',[xgr(i) q_gkw(i) shat_gkw(i) kappa_gkw(i) delta_gkw(i) square_gkw(i) skappa_gkw(i) sdelta_gkw(i) ssquare_gkw(i) zmil_gkw(i) drmil_gkw(i) dzmil_gkw(i) gradp_gkw(i)]); 
end; 
fprintf(fileID,'%s\n','#Species: m , Z, ngrid, Tgrid');
fprintf(fileID,'%12.8f %i %12.8f %12.8f\n',[me -1 negrid tegrid]); 
fprintf(fileID,'%s\n','# psi density Temperature RLn RLT'); 
for i = 1: n_x_grid 
  fprintf(fileID,'%12.8f %12.8f %12.8f %12.8f %12.8f\n',[xgr(i) ne_gkw(i) te_gkw(i) rlne_gkw(i) rlte_gkw(i)]); 
end; 
for i = 1: nions 
  fprintf(fileID,'%s\n','#Species: m , Z, Tgrid, ngrid');
  fprintf(fileID,'%12.8f %i %12.8f %12.8f\n',[mions(i) Z(i) tigrid(i) nigrid(i)]);
  fprintf(fileID,'%s\n','# psi density Temperature RLn RLT');
  for j = 1: n_x_grid
    fprintf(fileID,'%12.8f %12.8f %12.8f %12.8f %12.8f\n',[xgr(j) ni_gkw(i,j) ti_gkw(i,j) rlni_gkw(i,j) rlti_gkw(i,j)]);
  end;
end; 

fclose(fileID);

 
