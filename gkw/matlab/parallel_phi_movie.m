% This file defines a matlab function which plots the parallel potential as
% a movie against time.  This script assumes the naming conventions of the script gkwnlin.
% 
% parallel_phi(proj,runname)
%
% WHERE:
%
% proj and runname are strings for gkwpath.m
%
% n_av (integer) is the number of points to average over working backwards from
% the last point in the file. It is useful to check the saturation by eye
% in fluxes first.

function []=parallel_phi_movie(proj,filename)
%function [M]=parallel_phi_movie(proj,filename)

if ~exist('proj')
	proj='default';
    disp('You must provide the project name')
    return;
    %No default value
end
if ~exist('filename')
    disp('You must provide the run name')
    filename='default';
    return;
end

input=read_gkwinput(filename,proj);
ns=input.GRIDSIZE.n_s_grid
scale=[1:1:ns];
scale=-0.5+(scale-0.5)/ns;

data_file=load([gkwpath('parallel_phi',proj) filename]);
%data_file=load('parallel_phi.dat')
disp(['Loaded ' gkwpath('parallel_phi',proj) filename])

%figure(99)
max1=max(data_file);
max2=max(max1);
min1=min(data_file);
min2=min(min1);

axis manual

for j=1:size(data_file,1)

start=1;
finish=ns;    
    
label=[proj,' ',filename,' Norm:\phi'];
plot(scale(start:finish),data_file(j,start:finish),'+-','DisplayName',label)
axis([-0.5,0.5,min2,max2])

%M(j) = getframe;
 
pause(0.05)

end

end