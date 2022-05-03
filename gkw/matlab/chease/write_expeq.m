% writes an EXPEQ file (CHEASE units)
%	function []=write_expeq(in,flnm,flpth)
% Inputs:
%		in	structure containing the info to be stored in the EXPEQ 
%		flnm	file name for the EXPEQ file
%		flpth file path, if no path is entered, the default path is used
%
% Warning: no check of the structure given as an input. It has to have all the required fields

	function []=write_expeq(in,flnm,flpth)

% default path
if ~exist('flpth')||isempty(flpth)
	%flpth='/home/space/phsgbe/runs/chease/input/';
        flpath='./'
end

% opens file
if unix(['test -e ' flpth flnm])==0
	disp(['Warning, file ' flnm ' already exist.'])
	ok_write=input('Overwrite (0: no, 1: yes) ?');
	if ok_write~=1, 
		disp('Exit - nothing done')
		return
	end
end
fid = fopen([flpth flnm], 'w');

fprintf(fid,'%14.7f\n',in.aspct);
fprintf(fid,'%14.7g\n',in.ZRrat);
fprintf(fid,'%14.7g\n',in.pedge);
fprintf(fid,' %d\n',in.N_LCFS);
s=size(in.R);
if s(1)>1, in.R=in.R'; end
s=size(in.Z);
if s(1)>1, in.Z=in.Z'; end
fprintf(fid,'%14.7f  %14.7f\n',[ in.R ; in.Z]);
fprintf(fid,' %d\n',in.N_rho);
fprintf(fid,' %d\n',in.NSTTP);
fprintf(fid,'%14.7f\n',in.rho);
fprintf(fid,'%14.7f\n',in.pprime);
fprintf(fid,'%14.7f\n',in.j);
fclose(fid);
