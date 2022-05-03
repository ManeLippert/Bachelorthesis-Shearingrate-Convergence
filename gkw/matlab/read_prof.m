% reads the GKW output file prof_back.dat and builds a matlab structure with the information it contains
%	function [H]=read_prof(flnm,proj)
% Inputs:
%	flnm:	file name
%	proj:	project name (optional)
%		path for the input files is obtained from the gkwpath function with "proj" argument

%function [H]=read_prof(flnm,proj)

function [H] = read_prof(flnm)

% default path
%if ~exist('proj')
%	proj=[];
%    flpth='./';
%else
%  flpth=gkwpath('geom',proj);
%end

% reads file
%if unix(['test -e ' flpth flnm])==0
%	fid = fopen([flpth flnm], 'r');
%else	
%	error(['The file ' flpth flnm ' does not exist' ])
%end
fid = fopen(flnm,'r'); 
frewind(fid);

aaa = fscanf(fid, '%s',1);
countstr = [aaa];
while ~isempty(aaa);
  fpos=ftell(fid);    
  [zzz, nel] = fscanf(fid, '%f');
  fseek(fid,fpos,'bof'); 
  inum = size(findstr(countstr,aaa),2);
  H.(lower(aaa))(inum,:) = fscanf(fid, '%f',nel); 
  aaa = fscanf(fid, '%s',1);
  countstr = [countstr aaa]; 
end;

fclose(fid);
