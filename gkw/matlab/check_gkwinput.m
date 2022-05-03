% check the GKW input file (only quasi-neutrality, so far) and correct it if needed
%	function [out, msg]=check_gkwinput(flnm,proj)
% Inputs:
%		flnm	file name
%	proj:	project name (optional)
%		path for the input files is obtained from the gkwpath function with "proj" in argument
%


function [out]=check_gkwinput(flnm,proj)


% reads the input file
[GKWin,sss_ref]=read_gkwinput(flnm,proj);



%1) check the number of species and species fields

%2) compare the name of the fields of the ref. file to input.sample.dat

%3) check quasi-neutrality (including gradients)
n=GKWin.GRIDSIZE.number_of_species;

check=1e-8*sum([GKWin.SPECIES.dens]);
if abs(sum([GKWin.SPECIES.dens].*[GKWin.SPECIES.z]))>check,
	disp(['Quasineutrality not satisfied for file ' flnm])
end
if abs(sum([GKWin.SPECIES.rln].*[GKWin.SPECIES.z].*[GKWin.SPECIES.dens]))>check,
	disp(['Quasineutrality for the gradients not satisfied for file ' flnm])
end
