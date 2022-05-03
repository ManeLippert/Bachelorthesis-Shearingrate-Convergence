% This file defines a matlab function which can extract data for a single mode from the mess
% that is parallel.dat with multiple modes and plots.
% This script assumes the naming conventions of the script gkwnlin and
% requires gkwpath.m to be configured.
% 
% out=paralleldat(proj,runname,imod,ix,is,col,optplot,opteiv)
%
% col contents
% 1: sgr         the length along the field line 
% 2-3: phi       the potential: (re and im)
% 4-5: apar      the parallel component of the vector potential (re and im)
% 6-7: dens      the perturbed density (re and im)
% 8-9: tpar      the perturbed parallel temperature (re and im)
% 10-11: tperp   the perturbed perpendicular temperature (re and im)
% 12-13: wflow   the perturbed parallel flow velocity (re and im)
% 14-15: bpar    the perturbed parallel (compressional) magnetic feild (re and im)
%
% negative 2,4,6,8,10,12,14 will give abs(complex)
%
% WHERE:
%
% proj and runname are strings for gkwpath.m
% optplot 1 or zero turns the plot on or off
%
% YC: you may want to consider the modifications commented below. 
% The idea is to reshape data_file once for all to avoid loading multiple 
% instance of parallel.dat when you want to look at different modes
%
% FC: alternatively, use read_gkwparallel for reading in and reshaping the file once.


function [out]=paralleldat(proj,filename,imod,ix,is,col,optplot,opteiv)

%YC [kxrh,krho,out]=paralleldat(proj,filename)

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

if ~exist('optplot')
   optplot=1
end

% if ~exist(imod)
%     disp('You must select both modes')
%     return;
% end
% if ~exist(ix)
%     disp('You must select both modes')
%     return;
% end
% if ~exist(col)
%     disp('potential')
%     col=2
% end

input=read_gkwinput(filename,proj);
%input=read_gkwinput('input.dat');

ns=input.GRIDSIZE.n_s_grid
nx=input.GRIDSIZE.nx
nmod=input.GRIDSIZE.nmod
%YC nsp=input.GRIDSIZE.number_of_species;
kxrh(1,1:nx)=0;

scale=[1:1:ns];
scale=-0.5+(scale-0.5)/ns;

%load kxrh 
%load krho

if exist([gkwpath('kyspec',proj) filename '.krho'], 'file')
 krho=load([gkwpath('kyspec',proj) filename '.krho']);
 disp(['Loaded ' gkwpath('kyspec',proj) filename '.krho'])
else 
  krho(1,1)=input.MODE.kthrho
end
if exist([gkwpath('kxspec',proj) filename '.kxrh'],'file')
 kxrh=load([gkwpath('kxspec',proj) filename '.kxrh']);
disp(['Loaded ' gkwpath('kxspec',proj) filename '.kxrh'])
end 

disp(['selected mode: kx: ' num2str(kxrh(1,ix)) ', kyrho: ' num2str(krho(imod,1)) ])

if exist('opteiv','var')
  data_file=load([gkwpath('other',proj) filename '/parallel' num2str(opteiv,'%03i') '.dat']);  
  disp([gkwpath('other',proj) filename 'parallel' num2str(opteiv,'%03i') '.dat'])
else
  data_file=load([gkwpath('parallel',proj) filename]);
  disp(['Loaded ' gkwpath('parallel',proj) filename])
end
%data_file=load('parallel.dat');

% YC these two lines could replace the end of the file:
%     ncol=size(data_file,2);
%     out=reshape(data_file,[ns nx nmod nsp ncol]);
% 'out' is an array of size ns*nx*nmod*nsp*ncol
% Therefore to plot the imaginary part of the potential for mode imod, ix and the first species:
%     isp=1; col=3;
%     plot(tmp(:,ix,imod,isp,1),tmp(:,ix,imod,isp,col))

%calculate right rows of file
start_row=(is-1)*ns*nx*nmod+(imod-1)*ns*nx+(ix-1)*ns+1;
out(:,1)=data_file(start_row:start_row+ns-1,1)  ;

% restore the unrotated complex number
% rotate=get_rotate(proj,filename);
% rotate=rotate(1)+i*rotate(2);
% 
% complex=data_file(start_row:start_row+ns-1,col)+i*data_file(start_row:start_row+ns-1,col+1);  
% complex=complex*rotate;

%subplot(1,2,1)

label=[proj ' ' filename ' ix: ' num2str(ix) ' imod: ' num2str(imod) ' species: ' num2str(is) ];
if (col>0)
  if optplot
     plot(data_file(start_row:start_row+ns-1,1)-(nx-1)/2-1+ix,data_file(start_row:start_row+ns-1,col),'-','DisplayName',label)
     %hold all
     %plot(data_file(start_row:start_row+ns-1,1)-(nx-1)/2-1+ix,real(complex))
  end
  
out(:,2)=data_file(start_row:start_row+ns-1,col);
end

if (col<0)
  if optplot plot(data_file(start_row:start_row+ns-1,1)-(nx-1)/2-1+ix,(data_file(start_row:start_row+ns-1,-col).^2+data_file(start_row:start_row+ns-1,-col+1).^2).^(0.5),'-','DisplayName',label)    
  end
out(:,2)=(data_file(start_row:start_row+ns-1,-col).^2+data_file(start_row:start_row+ns-1,-col+1).^2).^(0.5);
end

%subplot(1,2,2)
%plot(abs(fft(out(:,2))),'+-')

function t=get_rotate(proj,file)
    
    [err str]=system(['grep rotated -A 1 ' gkwpath('out',proj) file ' | tail -n1']);
    %t=str(7:end-2)
    t=str2num(str(7:end-2));

end


end