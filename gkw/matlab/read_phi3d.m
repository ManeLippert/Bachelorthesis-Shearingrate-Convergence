% Reads the GKW 3d potential spectral output SPc3dXXXXX files and time
% averages, returning as a 3D array. Assumes directory structure of gkwnlin
% 
% usage: phi3d=read_phi3d(proj,runname,start);
%
% where proj is a string for the project folder name
% runname is the string for the run name
% start is an integer - file to start averaging from
% (you might want to skip those before nonlinear saturation)
%
% This script interfaces with Yann's scripts 
% It requires gkw_path to be set up correctly for your projects directory.

function [phi3d]=read_phi3d(proj,runname,start)

if ~exist('proj')
	proj='default';
    disp('You must provide the project name')
    %No default value
end
if ~exist('start')
    start=1;
    disp('Average over all data files')
end


files=dir([gkwpath('other',proj) runname '/Spc3d*']);
total = size(files)
start
input=read_gkwinput(runname,proj);
nx=input.GRIDSIZE.nx
ns=input.GRIDSIZE.n_s_grid
nmod=input.GRIDSIZE.nmod
kxrh=load([gkwpath('kxspec',proj) runname '.kxrh']);
krho=load([gkwpath('kyspec',proj) runname '.krho']);
phi3d=zeros(nmod,nx,ns);

for i = start:total(1)
   %i
   fid = fopen([gkwpath('other',proj) runname '/' files(i).name],'r');
   AA = fread(fid,'double');
   phi3d=phi3d+reshape(AA,nmod,nx,ns);   
   fclose(fid);
end

%normalise
phi3d=phi3d/(total(1));

end