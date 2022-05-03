% Interactive script to create the input files for a GKW scan
%  
%		 function []=create_gkwscan(proj,shortfile)
%
% Inputs:
%	shortfile:	0 or 1 to use or not the name of the scanned variable to name the files (default=1)
%	proj:	project name (optional)
%		path for the input files is obtained from the gkwpath function with "proj" argument
%
% Remarks: 
%	When you are asked for the reference input file, if you type 'return', the content of the default folder is displayed
%	When you are asked fot the variable names, if you type 'return', the content of the reference file is displayed.
%	You can enter the "scan of coupled variables" mode by entering a '+' after the name of the leading variable, you will then be asked about the other variables coupled to the leading one. 
%	When you are asked for the scanned values for a given variable, you can enter any valid matlab expression.
%	If the values for a variable are strings, they have to be in  a cell. Ex: {'.true.' '.false.'}
%	The number of scanned variables is not limited, but as the variable names are used to name the files, it can generate very long file names if a large number of variables are scanned. You can switch to short file names by calling the function with shortfile=1 (default is shortfile=1)
%	The name of the variable you want to scan has to be in the reference input file, otherwise the variable will be said not to be a GKW input
 


function []=create_gkwscan(proj,shortfile)

% default 
if ~exist('proj')
	proj=[];
end
if ~exist('shortfile')||isempty(shortfile)
	shortfile=1;
end
flpth=gkwpath('input',proj);


% asks for the reference file
flnm='';
while isempty(flnm)||unix(['test -e ' flpth flnm])~=0
	flnm=input(['Enter the name of the reference input file (if not located in ' flpth ', enter the full path): '],'s');
	I=findstr(flnm,'/');
	if ~isempty(I), 
		flpth=flnm(1:I(end));
		flnm=flnm(I(end)+1:end);
	end
	if isempty(flnm)||unix(['test -e ' flpth flnm])~=0,
		fprintf(1,'\n%s\n\n%s\n',['File ' flnm ' does not exist'],['List of input files in ' flpth ':'])		
		unix(['ls ' flpth]);
		fprintf(1,'\n');
	end
end

% loads the reference file and checks the number of species
[GKWin,sss_ref]=read_gkwinput(flnm,proj);
if ~(length(GKWin.SPECIES)-GKWin.SPCGENERAL.adiabatic_electrons==GKWin.GRIDSIZE.number_of_species)
	fprintf(1,'Invalid number of species fields in the reference file:\n \tNumber of species: %1d  (adiabatic=%1d)\n\tNumber of species fields: %1d\n',GKWin.GRIDSIZE.number_of_species,GKWin.SPCGENERAL.adiabatic_electrons,length(GKWin.SPECIES));
	%return
end
field=fieldnames(GKWin);

% store comments and replace them by !!
comments1=regexp(sss_ref,'![^\n]*','match');
sss_ref=regexprep(sss_ref,'![^\n]*','!1!');
sss_ref=[sprintf('/\n') sss_ref];
comments2=regexp(sss_ref,'/\s[^&]*&','match');
sss_ref=regexprep(sss_ref,'/\s[^&]*&','/\n!2!\n&');


% asks for the scan parameters
nbvar=input('Number of variables on which you want to make a scan: ');
disp('Name of these variables: ');
varname=cell(nbvar,1);
spcindx=cell(nbvar,1);
coupled=zeros(nbvar,1);
for ii=1:nbvar,
	[var_name{ii}, spcindx{ii}, coupled(ii)]=var_to_scan(ii,sss_ref,field,GKWin);
	if coupled(ii)~=0,
		disp(['You want to couple the scan of ''' var_name{ii} ''' to other variables'])
		coupled(ii)=input(['Number of variables coupled to ''' var_name{ii} ''': ']);		
		disp('Name of these variables: ');
		for jj=1:coupled(ii)
			[varcpl_name{ii,jj}, spcindxcpl{ii,jj}, tmp]=var_to_scan(char(jj+96),sss_ref,field,GKWin);
			if tmp~=0,
				disp(['Two level coupling not implemented, sorry. No variable will be coupled to ''' varcpl_name{ii,jj} ''''])
			end
		end
	end
end

% asks for the scan parameters values
nbval=zeros(1,nbvar);
for ii=1:nbvar,
	if spcindx{ii}==0
		sss='';
	else
		sss=num2str(spcindx{ii});
		sss=sss(~isspace(sss));
	end
	var_val{ii}=input(['Enter the scan values for variable ''' var_name{ii} sss ''': ']);
	disp(' ');
	nbval(ii)=length(var_val{ii});
	for jj=1:coupled(ii)
		if spcindxcpl{ii,jj}==0
			sss1='';
		else
			sss1=num2str(spcindxcpl{ii,jj});
			sss1=sss1(~isspace(sss1));
		end
		yes=0;
		while ~yes
			varcpl_val{ii,jj}=input(['Enter the scan values for variable ''' varcpl_name{ii,jj} sss1 ''' (coupled to ''' var_name{ii} '''): ']);
			disp(' ');
			if length(varcpl_val{ii,jj})==length(var_val{ii}),
				yes=1;
			else
				disp(['The number of scan values for ''' varcpl_name{ii,jj} sss1 ''' has to be the same than for ''' var_name{ii} sss '''']);
			end
		end
	end
end	
nbval_tot=prod(nbval);
disp(['A total of ' num2str(nbval_tot) ' GKW input files will be created']);
disp(' ');


% creation of the input files 
ok_write=1;
flroot = input('Enter the filename root:  ', 's');
indx=zeros(nbval_tot,nbvar);
for ii=1:nbvar, 
	% un peu de combinatoire....
	A=repmat([1:nbval(ii)]',prod(nbval(1:ii-1)),nbval_tot/prod(nbval(1:ii)))'; 
	indx(:,ii)=A(:);
end
for ii=1:nbval_tot,
	file=[flpth flroot];
	sss_new=sss_ref;
	for jj=1:nbvar	
		if shortfile~=1,
			file=[file '_' var_name{jj}(1:min(length(var_name{jj}),10)) num2str(indx(ii,jj),'%02i')];
        else
            if max(indx(:,jj))>9
			     file=[file '_' num2str(indx(ii,jj),'%02i')];
            else 
                 file=[file '_' num2str(indx(ii,jj))];
            end
		end	

		% put main variable values
		if iscell(var_val{jj}),
			rep=cell2mat(var_val{jj}(indx(ii,jj)));
		else
			rep=sprintf('%i',var_val{jj}(indx(ii,jj)));
		end
		if spcindx{jj}(1)~=0,
			for kk=1:length(spcindx{jj})
				sss_new=regexprep(sss_new,['(?<tok1>' var_name{jj} '[\s]*=[\s'']*)(?<tok2>[^\n!,;'']*)'],['$<tok1>' rep],spcindx{jj}(kk),'ignorecase');
			end
		else
			sss_new=regexprep(sss_new,['(?<tok1>' var_name{jj} '[\s]*=[\s'']*)(?<tok2>[^\n!,;'']*)'],['$<tok1>' rep],'once','ignorecase');
		end

		% put coupled variables values
		for kk=1:coupled(jj)
			if iscell(varcpl_val{jj,kk}),
				rep=cell2mat(varcpl_val{jj,kk}(indx(ii,jj)));
			else
				rep=sprintf('%i',varcpl_val{jj,kk}(indx(ii,jj)));
			end
			if spcindxcpl{jj,kk}(1)~=0,
				for ll=1:length(spcindxcpl{jj,kk})
					sss_new=regexprep(sss_new,['(?<tok1>' varcpl_name{jj,kk} '[\s]*=[\s'']*)(?<tok2>[^\n!,;'']*)'],['$<tok1>' rep],spcindxcpl{jj,kk}(ll),'ignorecase');
				end
			else
				sss_new=regexprep(sss_new,['(?<tok1>' varcpl_name{jj,kk} '[\s]*=[\s'']*)(?<tok2>[^\n!,;'']*)'],['$<tok1>' 	rep],'once','ignorecase');
			end			
		end
	end	
	for jj=1:length(comments2),	%put comments2 back
		sss_new=regexprep(sss_new,'/\n!2!\n&',comments2{jj},1);
	end
	for jj=1:length(comments1),	%put comments1 back
		sss_new=regexprep(sss_new,'!1!',comments1{jj},1);
	end
	sss_new(1:2)=[];
	if unix(['test -e ' file])==0 & ok_write~=2
		disp(['Warning, file ' file ' already exist.'])
		ok_write=input('Overwrite (0: no, 1: yes, 2: yes to all) ?');
	end
	if ok_write==1|ok_write==2,
		fid=fopen(file,'w');
		fprintf(fid,'%s',sss_new);
		fclose(fid);
	else
		disp(['File ' file 'not created'])
		ok_write=1;
	end
	sinfo.files{ii}=file;

	% check quasi-neutrality
	tmp=regexp(sinfo.files{ii},'[^/]*','match');
	check_gkwinput(tmp{end},proj);
end	


% creation of a scan description file 
iii=1;
ext=[];
while unix(['test -e ' gkwpath('scan',proj) flroot ext '.mat'])==0
	iii=iii+1;
	ext=num2str(iii);
	if iii>200,
	  return
	end
end
sinfo.date=date;
comments=input('Enter the comments to be inserted in the scan description file: ','s');
sinfo.comments=comments;
sinfo.nbfiles=nbval_tot;
sinfo.nbvar=nbvar;
sinfo.indx=indx;
for ii=1:nbvar,
	sinfo.var.name{ii}=var_name{ii};
	sinfo.var.val{ii}=var_val{ii};
	sinfo.var.nbval(ii)=nbval(ii);
	sinfo.var.spc{ii}=spcindx{ii};
	sinfo.var.coupled(ii)=coupled(ii);
	if coupled(ii)>0
	  for jj=1:coupled(ii)
		sinfo.var.cpl(ii).name{jj}=varcpl_name{ii,jj};
		sinfo.var.cpl(ii).val{jj}=varcpl_val{ii,jj};
		sinfo.var.cpl(ii).spc{jj}=spcindxcpl{ii,jj};
	  end
	else
		sinfo.var.cpl(ii).name='no coupled variable';
		sinfo.var.cpl(ii).val=[];
		sinfo.var.cpl(ii).spc=[];
	end	
end

save([gkwpath('scan',proj) flroot ext],'sinfo')
disp(['Description of the scan saved in ' gkwpath('scan',proj) flroot ext '.mat'])




%%%%%%%%%%%%%%%%%%%%%%%%%% sub-function to read the variable names
function [var_name,spcindx,coupled]=var_to_scan(ii,sss_ref,field,GKWin)

yes=0;
while ~yes,	
	var_name=input([num2str(ii) ')   '],'s');
	% check whether the variable name is empty
	if isempty(var_name),
		regexprep(regexprep(sss_ref,'!\d!',''),'\n[\s]*','\n') %displays the content of the ref file
		continue;
	end
	% check whether the variable is coupled to another one or not 
	if var_name(end)=='+'
		coupled=1;
		var_name=var_name(1:end-1);
	else
		coupled=0;
	end
	% find whether the variable name exists and if it is a species parameter
	for jj=1:length(field),	
		eval(['I=find(strcmpi(''' var_name ''',fieldnames(GKWin.' field{jj} ')));']);
		if ~isempty(I), 
			yes=1; % var_name{ii} is a GKW variable
			if strcmp('SPECIES',field{jj}), % test if it is a species variable
				disp(' ')
				disp(['''' var_name ''' is a species parameter. The species are:'])
				disp(' ')
				for kk=1:length(GKWin.SPECIES),
					str=sprintf('Species %i:\tZ=%i\t\tmass=%0.5e',kk,GKWin.SPECIES(kk).z,GKWin.SPECIES(kk).mass);
					disp(str)
				end
				disp(' ')	
				spcindx=input(['Ref. number of the species you want to scan (enter an array to scan several species simultaneously): ']);	
			else
				spcindx=0;
			end
		end
	end
	if yes==0,			
		disp(['Variable ''' var_name ''' is not a GKW input.'])
		disp(['Enter a new variable name or type ''return'' to display the content of the ref. input file: '])
	end
end
