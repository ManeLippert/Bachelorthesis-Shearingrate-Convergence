% reads the GKW output file geom.dat and builds a matlab structure with the information it contains
%	function [H]=read_geom(flnm,proj)
% Inputs:
%	flnm:	file name
%	proj:	project name (optional)
%		path for the input files is obtained from the gkwpath function with "proj" argument
% This routine is not written in a general way: it only works for a specific output file structure

function [H]=read_geom(flnm,proj)


% default path
if ~exist('proj')
	proj=[];
    flpth='./'
else
flpth=gkwpath('geom',proj);
end

% reads file
if unix(['test -e ' flpth flnm])==0
	fid = fopen([flpth flnm], 'r');
else	
	error(['The file ' flpth flnm ' does not exist' ])
end
frewind(fid);
sss='';


% scalars
file=dir([flpth flnm]);
if str2num(file.date(8:11))>2011
  nscal=15;  
else
  nscal=13;
end

for ii=1:nscal
	sss=deblank(fgets(fid));
	eval(['H.' lower(sss) '=fscanf(fid,''%f'',1);'])
	fgets(fid);
end

% 1_D quantities
for ii=1:27
	sss=deblank(fgets(fid));
	eval(['H.' lower(sss) '=fscanf(fid,''%f'',H.ns);'])
	fgets(fid);
end

fclose(fid);
