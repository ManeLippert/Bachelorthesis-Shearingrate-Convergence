% reads the content of an EXPEQ file (CHEASE units)
%	function [out]=read_expeq(flnm,flpth)
% Inputs:
%		flnm	file name for the EXPEQ file
%		flpth file path, if no path is entered, the default path is used



function [out]=read_expeq(flnm,flpth)

% default path
if ~exist('flpth')||isempty(flpth)
	%flpth='/home/space/phsgbe/runs/chease/input/';
        flpth='./'
end

% reads file
if unix(['test -e ' flpth flnm])==0
	fid = fopen([flpth flnm], 'r');
else	
	error(['The file ' flpth flnm ' does not exist' ])
end
frewind(fid);



% aspect ratio (Rmax-Rmin)/Rref/2
out.aspct=fscanf(fid,'%e',1);
% Zref/Rref
out.ZRrat=fscanf(fid,'%e',1);
% edge pressure
out.pedge=fscanf(fid,'%e',1);

% nb of points for LCFS
out.N_LCFS=fscanf(fid,'%5i',1);
% LCFS description
RZ=fscanf(fid,'%e',[2,out.N_LCFS])';
out.R=RZ(:,1)';
out.Z=RZ(:,2)';

% nb of points for profiles
out.N_rho=fscanf(fid,'%5i',1);
% NSTTP
out.NSTTP=fscanf(fid,'%5i',1);
% rho (rho_psi)
out.rho=fscanf(fid,'%e',out.N_rho);
% P'
out.pprime=fscanf(fid,'%e',out.N_rho);
% current 
out.j=fscanf(fid,'%e',out.N_rho);


fclose(fid);
