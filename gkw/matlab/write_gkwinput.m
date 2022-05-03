function write_gkwinput(instr, filename, proj, check, comment)
% Function that writes GKW input files from a matlab structure
% in the format produced by read_gkwinput
%
% usage:  write_gkwinput(struc, filename, proj, check, comment)
% e.g.    write_gkwinput(cyc_str, 'cyinput', 'cyclone', 1, 'AUG shot 27788')
%
% where:
%   structure - the matlab structure with the input data (examples in input_refs.mat)
%   filename  - the filename for the output input file
%   proj      - the project folder to write to (uses gkwpath.m and gkwnlin)
%   check     - optional integer to run the input checker executable on the file created
%   comment   - optional comment to insert into input file
%
% (the checker executable is built with "make checker"
%    executable will be $GKW_HOME/run/input_check*.x
%    script is in $GKW_HOME/scripts/gkw_check_input_new)

if (~exist('proj')==1)
    proj='writes';
end

if (~exist('check')==1)
    check=0;
end

if (~exist('comment')==1)
    comment='GKW Input file written by write_gkwinput.m';
end

pth=gkwpath('root',proj);
pthsim=gkwpath('input',proj);

%Check if the project directory exists
if (~exist(pth,'dir'))
    display([pth ' is not a directory!']);
    display(['use gkwnlin -p proj to generate']);
    return;
elseif (~exist(pthsim,'dir'))
    display([pthsim ' does not exist...']);
    display([pth ' is not a project directory?']);
    display(['use gkwnlin -p proj to generate']);
    return;
end

clear flnm
flnm = [pthsim filename];

[ya aaa] = unix(['ls ' flnm]);

if ya ~= 0;

else;
    disp(['ABORT: The simulation  ' [flnm] '  already exists.']);
    return
    yrpl=input('Do you want to replace the input file ? [yes/no]  ','s');
    if yrpl(1)~='y'; return; end;
end;

fptr = fopen(flnm,'w');

clear para paranm  namelists  equalctr;
frewind(fptr);

%Write the header comemnt
fprintf(fptr,'! %s\n',comment);

%Get namelist names
namelists=fieldnames(instr);

%loop over namelists
for jnm=1:length(namelists);

    %evaluate this namelist in the structure
    equalctr=eval(['instr.' namelists{jnm}]);

    for klm=1:length(equalctr);
        fprintf(fptr,'&%s\n',namelists{jnm});
        para=eval(['instr.' namelists{jnm} '(' num2str(klm) ')']);

        %Get variable names
        paranm=fieldnames(para);

        %loop over variables and write each one
        for ij=1:length(paranm);

            %numbers
            if isnumeric(getfield(para,paranm{ij}));
                lngtpara =  length(getfield(para,paranm{ij}));;
                if lngtpara == 1;
                    % If value is an integer, write an integer (assumes compiler promotes ok)
                    % note that all floats read by read_gkwinput have property isfloat
                    % but some longer integers also have isfloat, and not
                    % isinteger (why?)
                    % so cannot seperate this way
                    if round(getfield(para,paranm{ij})) == getfield(para,paranm{ij});
                        fprintf(fptr,' %s= %i,\n',paranm{ij},getfield(para,paranm{ij}));
                    else; %write a float
                        fprintf(fptr,' %s= %12.8f,\n',paranm{ij},getfield(para,paranm{ij}));
                    end;
                else
                    fprintf(fptr,' %s= ',paranm{ij});
                    dum=getfield(para,paranm{ij});
                    for jpn = 1:lngtpara;
                        fprintf(fptr,' %12.8f,',dum(jpn));
                    end
                    fprintf(fptr,'\n');
                end

            %logicals all have isinteger property if read_gkwinput was used
            elseif(islogical(getfield(para,paranm{ij})));
                if(getfield(para,paranm{ij}));
                    fprintf(fptr,' %s= %s,\n',paranm{ij},'.true.');
                else
                    fprintf(fptr,' %s= %s,\n',paranm{ij},'.false.');
                end

            %strings which are actually logicals    
            elseif(strcmp(getfield(para,paranm{ij}),'.false.') | strcmp(getfield(para,paranm{ij}),'.true.' ));
                fprintf(fptr,' %s= %s,\n',paranm{ij},getfield(para,paranm{ij}));

            %strings
            else;
                fprintf(fptr,' %s= ''%s'',\n',paranm{ij},getfield(para,paranm{ij}));
            end;
        end;
        fprintf(fptr,' /\n');
    end;
end;


fclose(fptr);
display([flnm ' written succesfully']);


% check the input using the executable made with "make checker"
% executable is in $GKW_HOME/run/input_check*.x
% script is in $GKW_HOME/scripts/gkw_check_input_new
if (check == 1)
    [a gkw_home] = unix('echo $GKW_HOME');
    gkw_home=gkw_home(1:end-1);
    unix([gkw_home '/scripts/gkw_check_input_new ' flnm]);
end

end