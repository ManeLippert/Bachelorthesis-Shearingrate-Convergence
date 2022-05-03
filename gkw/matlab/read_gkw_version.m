function [integer string] = read_gkw_version(file,proj)
% [integer string] =read_gkw_version(file,proj)
% reads the GKW version number used for a run from input.out
% uses the gkwnlin folder structure
    
    [err str]=system(['grep -i "GKW version" ' gkwpath('input_out',proj) file]);
    string = str(40:end);
    integer = str2num(str(40:44));

end