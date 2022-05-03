% This file defines a matlab function which plots time averaged kyspectrum
% on a log log plot, useful for nonlinear runs.
% 
% useage: kyspec(n_av,'file')
%
% Where n_av is the number of points to average over working backwards from
% the last point in the file. It is useful to check the saturation by eye
% in fluxes first.
%
% The spectrum in file is plotted against  file.krho values
% Hence this script assumes the naming conventions of the script gkwnlin
% If not provided, will try and read kyspec and krho as default filenames.

function []=kyspec(n_av,file)

if ~exist('n_av')
    n_av=100;
end

if ~exist('file')
	file='kyspec';
    file2='krho';
else
    file2=[file '.krho']
end

kyspc=load(file);
dim=size(kyspc);

for i=1:dim(2)
   av(i)=0; 
end    

for j=0:n_av-1
    for i=1:dim(2)
        av(i)=av(i)+kyspc(dim(1)-j,i);
    end
end

av(:)=av(:)/n_av;

scale=load(file2);

loglog(scale(2:dim(2),1),av(2:dim(2)),'+-','DisplayName',file)

% Create xlabel
xlabel('ky.rho');

% Create ylabel
ylabel('Amplitude');

end