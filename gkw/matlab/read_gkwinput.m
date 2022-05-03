% Reads the GKW input file and builds a matlab structure with the information it contains
%	function [GKWin file_text]=read_gkwinput(flnm,proj,out)
% Inputs:
%   	flnm:	file name
%       proj:	project name (optional) using structure of gkwnlin, else read from present directory
%       out:    integer (optional) if present, read input.out from specified project instead.
%               The GKW input.out file always contains a complete listing of all variables as used by the code.
%       Path for the input files is obtained from the gkwpath function with "proj" in argument.


function [GKW_in,sss_in]=read_gkwinput(flnm,proj,sw)

in_or_out='input';

%Optionally use output files from specified project.
if ~exist('sw')
    sw=0;
elseif (sw==1)    
    disp('Reading input.out...')
    in_or_out='input_out';
    % old gkwnlin format 
    if ~exist(gkwpath(in_or_out,proj),'dir')
       in_or_out='out'; 
       flnm=[flnm '.out'];
    end            
    sw=1;
elseif (sw==2)
    disp('Reading GS2 input');
end

% default path
if ~exist('proj')
	proj=[];
    flpth='./';
else
    flpth=gkwpath(in_or_out,proj);
end

if sw==2
   flpth=gs2path(in_or_out,proj);
end

% reads file
if unix(['test -e ' flpth flnm])==0
	fid = fopen([flpth flnm], 'r');
else	
	error(['The file ' flpth flnm ' does not exist' ])
end
frewind(fid);
sstmp='';
sss_in='';
while ischar(sstmp)
	sss_in=[sss_in sstmp];
	sstmp=fgets(fid);
end
fclose(fid);

% removes comments
sss=regexprep(sss_in,'![^\n]*',''); %removes all what is after a '!'
sss=[sprintf('/\n') sss];
sss=regexprep(sss,'/\s[^&]*&','/\n&'); % also removes what is between two namelist blocks, even if it is not commented

%forbidden character " in matlab eval function -> replace with '
sss=regexprep(sss,'"','''');
sss=regexprep(sss,')','');  % complex numbers
sss=regexprep(sss,'(','');

% finds the indexes of the different input blocks delimited by &..../
[I_start I_end]=regexp(sss,'&\w[^&]*','start','end');
field=regexp(sss,'&\w[^\s]*','match');


% builds the matlab structure

% loop over namelists
for ii=1:length(I_start),
	sub=regexp(sss(I_start(ii):I_end(ii)),'\<(?<field>[\w]*)[\s]*=[\s]*(?<value>[^\s,]*)','names');
	if exist('GKW_in')&isfield(GKW_in,field{ii}(2:end)),
		eval(['n=size(GKW_in.' field{ii}(2:end) ',2);'])
		n=n+1;
	else
        	n=1;
    end
    
    % loop over variables
	for jj=1:length(sub),
		if length(regexp(sub(jj).value,'\D'))==1
			if sub(jj).value(1)=='T'&sub(jj).value(end)=='T', % T values
				sub(jj).value='true';
			elseif sub(jj).value(1)=='F'&sub(jj).value(end)=='F', % F values
				sub(jj).value='false';
			end
		end
		
		if sub(jj).value(1)=='.'&sub(jj).value(end)=='.', % .false. and .true. values
			if strcmp(sub(jj).value(2),'t'),
				sub(jj).value='true';
			else
				sub(jj).value='false';
			end
		
		% detect string inputs that are not enclosed by '' or "" in the GKW input file. Look
		% at the number of characters that are not a digit or ./+/- in the string.
		% if it is bigger than 2, it is considered to be a string and avoids numbers (like 1.234e+05)
		% But it will fail if the path+filename is 1 character.
		% Could  also fail for some filenames containing many characters ./+/-/   
		%elseif length(regexp(sub(jj).value,'[A-Za-z]'))>2&&(sub(jj).value(1)~='''')
		elseif length(regexp(sub(jj).value,'[^\d\.\+\-]'))>1&&(sub(jj).value(1)~='''')
			sub(jj).value=['''' sub(jj).value ''''];
        end
        try	    
		    eval(['GKW_in.' field{ii}(2:end) '(n).' lower(sub(jj).field) '=' sub(jj).value ';'])
        catch            
            disp(['READ FAIL *** GKW_in.' field{ii}(2:end) '(n).' lower(sub(jj).field) '=' sub(jj).value ';'])
            %error('read_gkwinput: cannot eval the above command')
        end
	end

end

end
