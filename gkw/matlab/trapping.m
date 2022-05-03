% Calculates and plots trapping condition from GKW grid outputs.
% Assumes gkwnlin folder structure
%
% usage: trapping(proj,filename,species)

function []=trapping(proj,filename,species)

if ~exist('proj')
	proj='default';
    disp('You must provide the project name')
    return
end

if ~exist('filename')
    disp('You must provide the project name')
    return
end

if ~exist('species')
    species=1
    disp('Plotting for species 1')
end

flpth=gkwpath('other',proj)
%messy by works for now
direc=pwd;
cd([flpth filename]);
geom=read_geom('geom.dat');
cd (direc);
cfen=load([flpth filename '/cfen.dat']);
%vpar=load([gkwpath('distr1',proj) filename]);
%vperp=load([gkwpath('distr2',proj) filename]);
%vpar=vpar(1,:)
%vperp=vperp(1,:)
vpar=[-4.5:0.01:4.5];

cfenhfs=max(cfen(:,7+species))
cfenlfs=min(cfen(:,7+species))

func=((vpar.^2+(cfenlfs-cfenhfs))/(geom.bmax/geom.bmin-1)).^(0.5);
func=real(func);

plot(vpar,func,'w--','LineWidth',2)
%plot(vpar,func,'--','LineWidth',2)

end


