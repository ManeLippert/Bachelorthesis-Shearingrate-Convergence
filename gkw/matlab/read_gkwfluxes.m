% read the detailed fluxes output for a GKW run
% function [grids,pflux,eflux,vflux]=read_gkwfluxes(flnm,proj)
% Inputs:
%	flnm:	file name
%	proj:	project name (optional)
%		path for the input files is obtained from the gkwpath function with "proj" in argument
%	full:	full=1 -> load time arrays, full=0 -> load only the last time step
%
% Warning: in this version, all quantities are normalised as in GKW !!!!!!!

function [grids,pflux,eflux,vflux]=read_gkwfluxes(flnm,proj)

if ~exist('proj') 
	proj=[];
end
if exist('full')~=1 
	full=0;
end

flpth=gkwpath('fluxes_det',proj);

% read grids
if unix(['test -e ' flpth flnm '.g'])==0
	fid = fopen([flpth flnm '.g'], 'r');
else	
	error(['The file ' flpth flnm '.g does not exist' ])
end
frewind(fid);

nmod=fread(fid,1,'long');
nx=fread(fid,1,'long');
nsp=fread(fid,1,'long');
nmu=fread(fid,1,'long');
nvpar=fread(fid,1,'long');
ns=fread(fid,1,'long');

grids.kthrho = fread(fid,nmod,'double');
grids.kxrho = fread(fid,nx,'double');
grids.mu = fread(fid,nmu,'double');
grids.vpar = reshape(fread(fid,ns*nmu*nvpar,'double'),ns,nmu,nvpar);
grids.s = reshape(fread(fid,nx*ns,'double'),nx,ns);

grids.intmu = fread(fid,nmu,'double');
grids.intvp = reshape(fread(fid,ns*nmu*nvpar,'double'),ns,nmu,nvpar);
grids.ints = fread(fid,ns,'double');

fclose(fid);


% read file
if unix(['test -e ' flpth flnm])==0
	fid = fopen([flpth flnm], 'r');
else	
	error(['The file ' flpth flnm ' does not exist' ])
end
frewind(fid);

pflux = reshape(fread(fid,nmod*nx*nsp*ns*nmu*nvpar,'double'),nmod,nx,nsp,ns,nmu,nvpar);
eflux = reshape(fread(fid,nmod*nx*nsp*ns*nmu*nvpar,'double'),nmod,nx,nsp,ns,nmu,nvpar);
vflux = reshape(fread(fid,nmod*nx*nsp*ns*nmu*nvpar,'double'),nmod,nx,nsp,ns,nmu,nvpar);

fclose(fid);
