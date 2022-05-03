function outdens = read_gkwcfdens(sim_name, pth);
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
pthsim = [pth(1:end-1) '/cfdens/'];
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
outdens.s = fff(:,1);
%
[a b] = size(fff);
nb_of_species = (b-2)./2;
outdens.rln(:,1:nb_of_species) = fff(:,2:nb_of_species+1);
outdens.dens(:,1:nb_of_species)= fff(:,nb_of_species+2:2*nb_of_species+1);
%
outdens.poloidal_angle = fff(:,end);
outdens.sim_name = sim_name;
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%





