% Read all the output connected with a run 
% (spectra still need to be programmed, sorry ...) 
% Usage name = read_gkwrun(flnm,proj) 
% Input flnm : filename of the run 
%       proj : project name 
% Output name : structure that contains all the data 

function name = read_gkwrun(flnm,proj) 

% Read the input (always check for input.out)
flpth = gkwpath('out',proj);

if (exist([flpth flnm],'file')) 
  out = 1; 
  name = read_gkwinput(flnm,proj,out);
else 
  name = read_gkwinput(flnm,proj);
end

% Read the fluxes  
name = read_gkwflux(flnm,proj,name); 
