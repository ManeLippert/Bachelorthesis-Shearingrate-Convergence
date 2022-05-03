% load data from the various files of a GKW scan
%	function [vars,gamma,freq,pflux,eflux,vflux,s,phi,G,flist]=read_gkwscan(flnm,proj);
% Inputs:
%	flnm	file name for the scan summary (without the '.mat' extension)
%	proj	project name (optional)
%		path for the input files is obtained from the gkwpath function with "proj" in argument
%
% Outputs:
%	G	structure with the information of the GKW input file
%	vars	structure with the values of the scanned variables
%		If the variable is a species variable, the numbers corresponding to the scanned species are appended to the variable name.
%		ex: "vars.rlt12" means that the variable "rlt" has been scanned for species "1" and "2"
%	gamma	growth rate
%	pflux	particule flux (last dimension is the species dimension)
%	eflux	heat flux (last dimension is the species dimension)
%	vflux	parallel momentum flux (last dimension is the species dimension)
%	s	parallel coordinate grid (last dimension is the grid length)
%	phi	electrostatic potential (before last dimension is the grid length, last dimension is real/imaginary part)
%	flist	cell array with file names to help making the correspondance between a point of the scan and the name of the file loaded
%Note1: the fluxes are still semi-normalised, i.e. the species dependent normalisation (factors n_s, n_s*T_s and m_s*n_s*vth_s for pflux, eflux, vflux, respectively) has been removed but the vthref*rhostar^2 factor is still there. 
%Note2: the electrostatic potential is rotated in the complex plane such as to have real(phi(s=0))=max(real(phi)) and imag(phi(s=0))=0. Same representation than GS2.

function [vars,gamma,freq,pflux,eflux,vflux,s,phi,G,flist]=read_gkwscan(flnm,proj)


% default 
if ~exist('proj')
	proj=[];
end
flpth=gkwpath('scan',proj);

% load the scan summary
if unix(['test -e ' flpth flnm '.mat'])~=0
	disp(' ')
	disp(['The file ' flpth flnm '.mat does not exist.' ])
	disp('List of available scans for that project:')
	disp(' ')
	list_gkwscan(proj)
	error(' ')
end
load([flpth flnm],'sinfo');
ntot=sinfo.nbfiles;


% get the filenames
for ii=1:ntot,
	tmp=regexp(sinfo.files{ii},'[^/]*','match');
	fl{ii}=tmp{end};
end


% build the arrays for the scanned variables
for ii=1:sinfo.nbvar
	if sinfo.var.spc{ii}==0
		sss='';
	else
		sss=num2str(sinfo.var.spc{ii}(:)');
		sss=sss(~isspace(sss));
	end
	eval(['vars.' lower(sinfo.var.name{ii}) sss '=sinfo.var.val{ii};'])
end

% add the coupled variable info
for ii=1:sinfo.nbvar
	if isfield(sinfo.var,'coupled')
	  for jj=1:sinfo.var.coupled(ii)
		eval(['vars.' lower(sinfo.var.name{ii}) sss '_cpl.name=sinfo.var.cpl(ii).name;'])
		eval(['vars.' lower(sinfo.var.name{ii}) sss '_cpl.val=sinfo.var.cpl(ii).val;'])
		eval(['vars.' lower(sinfo.var.name{ii}) sss '_cpl.spc=sinfo.var.cpl(ii).spc;'])
	  end
	end	
end


% load data
dimprod=cumprod([1 sinfo.var.nbval]);
sinfo.indx(:,1)=sinfo.indx(:,1)+1;
indx=sum((sinfo.indx-1).*repmat(dimprod(1:end-1),ntot,1),2);
G=read_gkwinput(fl{1},proj);
nbspecies=G.GRIDSIZE.number_of_species;
ns=G.GRIDSIZE.n_s_grid;

[gamma,freq]=deal(zeros([sinfo.var.nbval 1]));
[pflux,eflux,vflux]=deal(zeros([sinfo.var.nbval nbspecies]));
%[s]=deal(zeros([sinfo.var.nbval ns])); % does not work if N_s_grid is scanned
%[phi]=deal(zeros([sinfo.var.nbval ns 2]));
s=cell([sinfo.var.nbval 1]);
phi=cell([sinfo.var.nbval 1]);
flist=cell([sinfo.var.nbval 1]);

for ii=1:ntot,
  if unix(['test -e ' gkwpath('time',proj)  fl{ii} ])==0

  G=read_gkwinput(fl{ii},proj);
  nbspecies=G.GRIDSIZE.number_of_species;
  spcnumber=[1:length([G.SPECIES.z])];
  spcnumber=spcnumber([G.SPECIES.z]~=-1|~G.SPCGENERAL.adiabatic_electrons);
  ns=G.GRIDSIZE.n_s_grid;

  % growth rate
  if G.CONTROL.non_linear==0
     tmp=load([gkwpath('time',proj)  fl{ii}]);
     gamma(indx(ii))=tmp(end,2);
     if size(tmp,2)>2 
      freq(indx(ii))=tmp(end,3);
     else
      freq(indx(ii)) = NaN;
     end
  else 
     gamma(indx(ii))=NaN;
     freq(indx(ii)) = NaN;
  end

  % fluxes
  tmp=load([gkwpath('fluxes',proj)  fl{ii}]);
  pnorm=[G.SPECIES.dens];
  pflux(indx(ii)+dimprod(end)*[0:nbspecies-1])=tmp(end,1:3:end).*pnorm(spcnumber);
  enorm=[G.SPECIES.temp].*[G.SPECIES.dens];
  eflux(indx(ii)+dimprod(end)*[0:nbspecies-1])=tmp(end,2:3:end).*enorm(spcnumber);
  vnorm=[G.SPECIES.mass].*[G.SPECIES.dens].*sqrt(2*[G.SPECIES.temp]./[G.SPECIES.mass]);
  vflux(indx(ii)+dimprod(end)*[0:nbspecies-1])=tmp(end,3:3:end).*vnorm(spcnumber);

  % phi
  tmp=load([gkwpath('parallel',proj)  fl{ii}]);
  s_tmp=tmp(1:ns,1);
  phi_tmp=tmp(1:ns,2:3);
  % normalize phi to be in the same configuration than GS2
  I=find(-1e-6<s_tmp&s_tmp<1e-6);
  phi_0=phi_tmp(I,1)+i*phi_tmp(I,2);
  A=abs(phi_0);
  B=phase(phi_0);
  phi_rot=(phi_tmp(:,1)+i*phi_tmp(:,2)).*exp(-i*B)./A;
  phi_tmp(:,1)=real(phi_rot);
  phi_tmp(:,2)=-imag(phi_rot);
  % fill the arrays
%  s(indx(ii)+dimprod(end).*[0:ns-1])=s_tmp;
%  phi(indx(ii)+dimprod(end).*[0:ns-1])=phi_tmp(:,1);
%  phi(indx(ii)+dimprod(end).*[0:ns-1]+dimprod(end).*ns)=phi_tmp(:,2);
  s{indx(ii)}=s_tmp;
  phi{indx(ii)}=phi_tmp;

  end

  % files list
  flist{indx(ii)}=fl{ii};
end


% display some information
fprintf('\nScan "%s" of project "%s" loaded',flnm,proj)
fprintf('\n%i files, %i variables scanned, created on %s\n\n',ntot,sinfo.nbvar,sinfo.date)

	