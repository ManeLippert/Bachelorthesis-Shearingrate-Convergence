function [p,e]=parity(proj,filename,imod)
% returns the tearing parity parameter and EM ratio calculated from parallel.dat
%
% [p,e]= parity(file,proj,imod) 
%
% where
%     p =  | int complex_apar ds | / int  | complex_apar | ds
%     could alternatively be integrated over theta if making direct comprison with other codes
%
%     e =  int | complex_apar | ds | / int |complex_phi | ds
%
% if imod is negative, use values from the eigenvalue solver
%
% FJC 7.12.12

if ~exist('proj')
	proj='default';
    disp('You must provide the project name')
    return;
    %No default value
    p=NaN;
end
if ~exist('filename')
    disp('You must provide the run name')
    filename='default';
    return;
    p=NaN;
end
if ~exist('imod')
    imod=1;
end

is=1;
eiv=0;
if (imod < 0)
    imod = abs(imod);
    eiv = 1;
end
    
if (exist([gkwpath('parallel',proj) filename],'file')==2 | eiv == 1)
    
    input=read_gkwinput(filename,proj);
    ns=input.GRIDSIZE.n_s_grid;
    nx=input.GRIDSIZE.nx;
    nmod=input.GRIDSIZE.nmod;
    ix=(nx+1)/2;
    
    %calculate right rows of file
    start_row=(is-1)*ns*nx*nmod+(imod-1)*ns*nx+(ix-1)*ns+1;
    
    if (eiv == 0)
        data=load([gkwpath('parallel',proj) filename]);
        %data=load('parallel.dat')
        disp(['Loaded ' gkwpath('parallel',proj) filename])
    elseif(eiv == 1)
        data=load([gkwpath('other',proj) filename '/parallel' num2str(imod,'%03i') '.dat']);
        disp(['Loaded ' gkwpath('other',proj) filename '/parallel' num2str(imod,'%03i') '.dat'])
    end
        
    % AVOID using i as a variable.  It is the imaginary unit.
    clear i
    apar=data(start_row:start_row+ns-1,4)+i*data(start_row:start_row+ns-1,5);
    phi=data(start_row:start_row+ns-1,2)+i*data(start_row:start_row+ns-1,3);
    sgrid=data(start_row:start_row+ns-1,1);
    %theta=get from geom.dat
    
    p=abs(trapz(apar,sgrid))/trapz(abs(apar),sgrid);
    e=trapz(abs(apar),sgrid)/trapz(abs(phi),sgrid);
    
else
    p=NaN;
    e=NaN;
end

end


