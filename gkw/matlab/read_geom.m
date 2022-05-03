% reads the GKW output file geom.dat and builds a matlab structure with the information it contains
%	function [H]=read_geom(flnm,proj)
% Inputs:
%	flnm:	file name
%	proj:	project name (optional - if not use full path)
%		path for the input files is obtained from the gkwpath function with "proj" argument

function [H]=read_geom(flnm,proj)

% default path
if ~exist('proj')
	proj=[];
    flpth='';
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

aaa = fscanf(fid, '%s',1);

while ~isempty(aaa);
  fpos=ftell(fid);    
  [zzz, nel] = fscanf(fid, '%f');
  fseek(fid,fpos,'bof'); 
  H.(lower(aaa)) = fscanf(fid, '%f',nel); 
  aaa = fscanf(fid, '%s',1);
end;

H.filename = flnm;
H.proj = proj;

fclose(fid);
