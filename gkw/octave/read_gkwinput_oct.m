% Reads the GKW input file and builds a matlab structure with the information it contains
%   function [GKW_in, file_text] = read_gkwinput_oct(flnm, proj, out)
%
% Inputs:
%       flnm:   file name
%       proj:   project name (optional) using structure of gkwnlin, else read from present directory
%       out:    integer (optional) if present, read input.out from specified project instead.
%               The GKW input.out file always contains a complete listing of all variables as used by the code.
%       Path for the input files is obtained from the gkwpath function with "proj" in argument.


function [GKW_in,sss_in]=read_gkwinput_oct(flnm, proj, sw)

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

  % Read the file, the complete content is stored in  sss_in
  % Open the file
  if unix(['test -e ' flpth flnm])==0
    fid = fopen([flpth flnm], 'r');
  else
    error(['The file ' flpth flnm ' does not exist' ])
  end
  % Initialize the values
  frewind(fid);
  sstmp='';
  sss_in='';
  % Read the content
  while ischar(sstmp)
    sss_in=[sss_in sstmp];
    sstmp=fgets(fid);
  end
  % Close the file.
  fclose(fid);

  % Now start proccessing of the content of the file.

  % removes comments
  sss=regexprep(sss_in,'![^\n]*',''); %removes all what is after a '!'
  sss=[sprintf('/\n') sss];
  sss=regexprep(sss,'/\s[^&]*&','/\n&'); % also removes what is between two namelist blocks, even if it is not commented

  %forbidden character " in matlab eval function -> replace with '
  sss=regexprep(sss,'"','''');
  sss=regexprep(sss,'\)','');  % complex numbers
  sss=regexprep(sss,'\(','');

  % finds the indexes of the different input blocks delimited by &..../
  [I_start I_end]=regexp(sss,'&\w[^&]*','start','end');
  field=regexp(sss,'&\w[^\s]*','match');

  % builds the matlab structure

  % loop over namelists
  for ii=1:length(I_start)
    sub=regexp(sss(I_start(ii):I_end(ii)),'[\s]*(?<field>[\w\d]*)[\s]*=[\s]*(?<value>[^\s,]*)','names');

    if exist('GKW_in') && isfield(GKW_in,field{ii}(2:end)),
      eval(['n=size(GKW_in.' field{ii}(2:end) ',2);']);
      n=n+1;
    else
      n=1;
    end

    % If only one element in the namelist, the result (field/value) is not recognized as cell array. Thus it
    % has to be manually transformed to a cell array.
    if (! iscell(sub.field))
      sub.field = mat2cell(sub.field,1);
      sub.value = mat2cell(sub.value,1);
    endif
    % loop over variables
    for jj=1:length(sub.value)
      if length(regexp(sub(1).value{jj},'\D'))==1
        if sub.value{jj}(1)=='T' && sub.value{jj}(end)=='T', % T values
          sub.value{jj}='true';
        elseif sub.value{jj}(1)=='F' && sub.value{jj}(end)=='F', % F values
          sub.value{jj}='false';
        end
      end

      if strcmp(sub.value{jj}, '.true.') || strcmp(sub.value{jj}, '.false.')
        if strcmp(sub.value{jj}(2),'t'),
          sub.value{jj}='true';
        else
          sub.value{jj}='false';
        end

      % detect string inputs that are not enclosed by '' or "" in the GKW input file. Look
      % at the number of characters that are not a digit or ./+/- in the string.
      % if it is bigger than 2, it is considered to be a string and avoids numbers (like 1.234e+05)
      % But it will fail if the path+filename is 1 character.
      % Could  also fail for some filenames containing many characters ./+/-/
      %elseif length(regexp(sub(jj).value,'[A-Za-z]'))>2&&(sub(jj).value(1)~='''')
      elseif length(regexp(sub.value{jj},'[^\d\.\+\-]', 'match')) > 1 && (sub.value{jj}(1)~='''')
        sub.value{jj}=['''' sub.value{jj} ''''];
      end

      eval(['GKW_in.' field{ii}(2:end) '(n).' lower(sub.field{jj}) '=' sub.value{jj} ';']);
    end
  end

end
