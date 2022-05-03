function outdens = read_gkwcfphi(sim_name, pth);
%
% function outdens = read_gkwcfdens(sim_name, pth);
%

if ~exist('pth');
[a pth] = unix('echo $GKWMAT_DIR');
end
if isempty(pth);
clear pth
[a pth] = unix('echo $GKWMAT_DIR');
end
pthsim = [pth(1:end-1) '/cfphi/'];
pthsimgeo = [pth(1:end-1) '/geom/'];

clear flnm
flnm = [pthsim sim_name];
fptr = fopen(flnm,'r');
if fptr == -1
disp(['No such simulation in ' pth]);
disp(['These are the simulations available']);
[a avsims] = unix(['ls ' pthsim]);
disp(avsims)
ystr = -1;
error(sim_name)
%return;
end;
fclose(fptr);

eval(['load ' flnm ';']);
eval(['fff = ' sim_name ';']);
%
size(fff);
%
%outdens.s = fff(:,1);
%
[a b] = size(fff);
%nb_of_species = (b-2)./2;
outdens.phi = fff(:,1);
outdens.dphi_dpsi= fff(:,2);
%
outdens.sim_name = sim_name;
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%





