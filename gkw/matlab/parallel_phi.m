% This file defines a matlab function which plots the time averaged parallel potential
% useful for nonlinear runs.  This script assumes the naming conventions of the script gkwnlin.
% 
% parallel_phi(proj,runname,start,fin,norm)
%
% WHERE:
%
% proj and runname are strings for gkwpath.m
% Norm is an integer - 1 to normalise
%
% n_av (integer) is the number of points to average over working backwards from
% the last point in the file. It is useful to check the saturation by eye
% in fluxes first.

function [av]=parallel_phi(proj,filename,beg,fin,norm)


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
if ~exist('n_av')
    n_av=100;
    disp('Average over last 100 data points')
end
if ~exist('norm')
    norm=0
end

input=read_gkwinput(filename,proj);
ns=input.GRIDSIZE.n_s_grid
scale=[1:1:ns];
scale=-0.5+(scale-0.5)/ns;

data_file=load([gkwpath('parallel_phi',proj) filename]);
%data_file=load('parallel_phi.dat')
disp(['Loaded ' gkwpath('parallel_phi',proj) filename])

av=mean(data_file(beg:fin,:));
error=std(data_file(beg:fin,:))/sqrt(n_av);

% for i=1:ns
% norm(:,i)=sqrt(data_file(:,i).^2+data_file(:,i+ns).^2);
% end
% 
% av_norm=mean(norm(end-n_av:end,:));
% error_norm=std(norm(end-n_av:end,:))/sqrt(n_av);

start=1;
finish=ns;

if norm==1
   normal=sum(av(start:finish))
   %norm=sum(norm(:))
else
  normal=1000 
end

%input.ROTATION.shear_rate

%label=[proj,' ',filename,' Norm:\phi, \gamma_E:', num2str(input.ROTATION.shear_rate)]
label=[proj,' ',filename]
%errorbar(scale(start:finish),av(start:finish),error(start:finish),'+-','DisplayName',label)
plot(scale(start:finish),av(start:finish)/normal,'+-','DisplayName',label)

axis([-0.5,0.5,0,1.2*max(av)/normal])

hold all;


% Create xlabel
xlabel('s');
legend('hide')
legend('show')

% Create ylabel
ylabel('A. U.');
title('|\phi|^2 Parallel potential');
 
    
end