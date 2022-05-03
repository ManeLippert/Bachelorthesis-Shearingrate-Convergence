% Reads the GKW 2D files when output as binary (xy_bin = T)
% with stacked time slices in one file         (xy_estep = T)
%
% usage:
% dat=read_file2d(ny,nx,nt,filename)
%
% %example 1 (real space):
% load xphi
% load yphi
% load time.dat
% dat=read_file2d(size(xphi,1),size(xphi,2),size(time,1),'den01_bin')
% contourf(xphi,yphi,dat(:,:,1))
%
% %example 2 (velocity space):
% load distr1.dat % vpar grid
% load distr2.dat % vperp grid
% load time.dat
% dat2=read_file2d(size(distr1,2),size(distr1,1),size(time,1),'efluxes_vspace01_bin')
% contourf(distr1',distr2',dat2(:,:,1))
%
% WARNING: Endianness of the data might need correcting in some situations
% see doc fread for how to do this within matlab
% Also assumes output unformatted output written by fortran
% is double precision with 32 bit header integers
%
% To DO: Allow seek to read only correct frames
% To DO: HDF5 would make this obselete...!

function dat=read_file2d(ny,nx,nt,file)

 %i
 fid = fopen(file);
 
 % Preallocate for speed
 dat=NaN*zeros(ny,nx,nt);
 
 for it=1:nt
   n1  = fread(fid,1,'int32');
   AA  = fread(fid,ny*nx,'double');
   size(AA);
   try
      n2  = fread(fid,1,'int32');
   catch
     it  
     error('read_file2d: Cannot read timestep number above')  
   end
   dat(:,:,it)  = reshape(AA,ny,nx);    
 end   
 
 fclose(fid);
end

