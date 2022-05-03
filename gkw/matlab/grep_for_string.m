function [present str] = grep_for_string(proj,file,string)
% function [present str] = grep_for_string(proj,file,string)
% check for given string in gkw screen output file
    
    present = false;
    
    fil = [gkwpath('out',proj) file];
    
    if ~(exist(fil,'file'))
        error(['File not found: ' fil])
    end

    [err str]=system(['grep -i "' string '" ' fil]);
    if (size(str>0) & err==0); present = true; end

end