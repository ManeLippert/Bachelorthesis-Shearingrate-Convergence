%This matlab function the plots the 2D velocity space diagnostic
%abs(f(vparallel,vperp) of mode(1,1) for one species at one s point.
%only meaningful for mode_box false runs
%
%Usage: velplot(proj,filename,species)
%
% Species optional (default 1)
% Using a negative species will instead plot a 1D cross cut for the first mu point
%
%This script gkwpath to be set up correctly for your projects directory.
%
% To overplot the (centrifugal) trapping condition:
% hold all
% trapping(proj,filename,species)

function []=velplot(proj,filename,spc)

if ~exist('spc')
    spc=1;
end

if spc==0
  species='';    
elseif spc<10
  species=['.sp0' num2str(abs(spc)) '.dat']      
else
    disp('Can only manage first 9 species')
    return
end 
        

if ~exist('proj')
	proj='default';
    disp('You must provide the project name')
    %No default value
end

if ~exist('filename')
    filename='';
end

input=read_gkwinput(filename,proj);

d1=load([gkwpath('distr1',proj) filename]);
d2=load([gkwpath('distr2',proj) filename]);
d3=load([gkwpath('distr3',proj) filename species]);
d4=load([gkwpath('distr4',proj) filename species]);

if (spc>0) figure; end
set(gcf,'defaulttextinterpreter','Tex')

%Locked z axis, fdisi/phi
%contourf(d1,d2,sqrt(d3.*d3+d4.*d4),[0:0.0002:0.0026],'DisplayName',filename);
%Locked z axis, fdisi
%contourf(d1,d2,sqrt(d3.*d3+d4.*d4),[0:0.0001:0.0026],'DisplayName',filename);
%contourf(d1,d2,sqrt(d3.*d3+d4.*d4),[0:0.00001:0.0016],'DisplayName',filename);
%Auto z axis
if (spc>0); 
  contourf(d1,d2,sqrt(d3.*d3+d4.*d4),20,'DisplayName',filename); 
  %contourf(d1,d2,d3.*d3-d4.*d4,20,'DisplayName',filename);
  ylabel('v_\perp/v_{th}');
else
  %plot(d1(1,:),sqrt(d3(1,:).*d3(1,:)+d4(1,:).*d4(1,:)))
  plot(d1(1,:),d3(1,:).*d3(1,:)-d4(1,:).*d4(1,:))
  ylabel('f (A.U.)');
end

xlabel('v_{||}/v_{th}');

%xlim([-1.5 1.5]);
%ylim([0 2.5]);
title(['u = ' num2str(input.ROTATION.vcor) ])

end
