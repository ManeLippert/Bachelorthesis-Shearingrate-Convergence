% modify a GKW input file 
%	function modify_gkwinput(flnm,proj,var_name,var_val,spcindx)
% Inputs:
%	flnm:		file name
%	proj:		project name (optional)
%			path for the input files is obtained from the gkwpath function with "proj" in argument
%	var_name:	name of the variables to modify (cell array)
%	var_val:	values associated to the variables to modify (cell array)
%	spcindx:	species index associated to the variables to modify (cell array). Zero it the variable is not a species variable.
%
% Warning: the existence of the variable in the input file is not checked: if the variable does not exist, nothing happens

function modify_gkwinput(flnm,proj,var_name,var_val,spcindx)

% default 
if ~exist('proj')
	proj=[];
end
flpth=gkwpath('input',proj);
if ~exist('spcindx')
	spcindx=cell(size(var_name));
	for ii=1:length(var_name)
	  spcindx{ii}=0;
	end
end

% check
if length(var_name)~=length(var_val)|length(var_name)~=length(spcindx),
	error('Inputs ''var'', ''val'' and ''spcindx'' must have the same length')
end

% loads the reference file 
[GKWin,sss_ref]=read_gkwinput(flnm,proj);
field=fieldnames(GKWin);

% store comments and replace them by !!
comments1=regexp(sss_ref,'![^\n]*','match');
sss_ref=regexprep(sss_ref,'![^\n]*','!1!');
sss_ref=[sprintf('/\n') sss_ref];
comments2=regexp(sss_ref,'/\s[^&]*&','match');
sss_ref=regexprep(sss_ref,'/\s[^&]*&','/\n!2!\n&');


% file modification
sss_new=sss_ref;
for ii=1:length(var_name),
	if iscell(var_val{ii}),
		rep=cell2mat(var_val{jj});
	else
		rep=sprintf('%e',var_val{ii});
	end
	if spcindx{ii}(1)~=0,
		for kk=1:length(spcindx{ii})
			sss_new=regexprep(sss_new,['(?<tok1>' var_name{ii} '[\s]*=[\s'']*)(?<tok2>[^\n!,;'']*)'],['$<tok1>' rep],spcindx{ii}(kk),'ignorecase');
		end
	else
		sss_new=regexprep(sss_new,['(?<tok1>' var_name{ii} '[\s]*=[\s'']*)(?<tok2>[^\n!,;'']*)'],['$<tok1>' rep],'once','ignorecase');
	end
end	

% put comments back
for jj=1:length(comments2),	%put comments2 back
	sss_new=regexprep(sss_new,'/\n!2!\n&',comments2{jj},1);
end
for jj=1:length(comments1),	%put comments1 back
	sss_new=regexprep(sss_new,'!1!',comments1{jj},1);
end
sss_new(1:2)=[];

% write file
file=[flpth flnm];
fid=fopen(file,'w');
fprintf(fid,'%s',sss_new);
fclose(fid);
