% This script is for making larger 2D movies out or larger datasets.
% Requires lots of memory and CPU time
% The avi file write works for large files but needs a version of matlab > 2010b 
% Define start and fin as integers to select first and last files to read
%
% FJC. 11.11.12

% Do not CD to a directory with 1000's of files in it 
% because matlab will (slowly) search it in path !
dir2='./spectrum/phi/run_name'
dir='./other/run_name'

xphi=load([dir2 '/xphi']);
xphi=circshift(xphi,[0 size(xphi,2)/2])


yphi=load([dir2 '/yphi']);
xmin = min(min(xphi));
xmax = max(max(xphi));
ymin = min(min(yphi));
ymax = max(max(yphi));

%clear M;
% To reset the color scales
%clear min1 min2 min3 min4 min5 shift1 shift2 shift3 shift4 shift5
close all; 

% Change which files to plot here
file=[dir '/PFlesr0003_']
file2=[dir2 '/phi']
file3=[dir '/den02_']
file4=[dir '/ene02_']
file5=[dir '/pac02_']
file6=[dir '/ene01_']

% join 5 colormaps together (matlab only allows one colormap per figure)
map1=colormap('jet');
map2=colormap('gray');
map3=colormap('bone');
map4=colormap('copper');
map5=colormap('hot');
map6=colormap('pink');
cmap=[map1; map2; map3; map4; map5; map6];

% Use something like this is you don't want to renormalise the data in time
%   dat=load(name);
%   maxv=max(maxv,max(max(dat)));
%   minv=min(minv,min(min(dat)));
%   norm = sum(sum(dat));
%   figure(2);
%   plot([i],[norm],'bo');
%   end;
%   maxv = max(maxv,abs(minv));
%   vv = linspace(-1.1*maxv, 1.1*maxv, 20);

% Use to make the figure a bit smaller
sfac=1.0;
scrsz = get(0,'ScreenSize');
h=figure('Position',[1 scrsz(2)*sfac scrsz(3)*sfac scrsz(4)*sfac]);
colormap(cmap)
ha=tight_subplot(2,3,[0 0.01],[0.03 0.03],[0.03 0.01]);

for j= start:fin
  i = j
  il = j - start + 1;
  %il = j;
 
  str=sprintf('%06i',i)
  str2=sprintf('%04i',i)
  
  try
   dat=load([file str2]);
  catch
   dat=dat_prev;
  end

  try
    dat2=load([file2 str]);
  catch  
    dat2=dat_prev;
  end
  
  try 
   dat3=load([file3 str]);
  catch
   dat3=dat3_prev;
  end  
  
  try
    dat4=load([file4 str]);
  catch
    dat4=dat4_prev;  
  end
  
  try
    dat5=load([file5 str]);
  catch
    dat5=dat5_prev;
  end
  
  try
    dat6=load([file6 str]);
  catch
    dat6=dat6_prev;  
  end
  
  dat_prev=dat;
  dat2_prev=dat2;
  dat3_prev=dat3;
  dat4_prev=dat4;
  dat5_prev=dat5;
  dat6_prev=dat6;
  
  
  % shift all datasets onto same colormap 
  % non overlapping contiguous dataets each with range 1-2*shift
  shift=0.021;
  fac=1-2*shift;
  levels=48;
  
  % Use rolling average to keep the color scales smooth
  max1(mod(i,5)+1)=max(max(dat));
  max2(mod(i,5)+1)=max(max(dat2));
  max3(mod(i,5)+1)=max(max(dat3));
  max4(mod(i,5)+1)=max(max(dat4));
  max5(mod(i,5)+1)=max(max(dat5));
  max6(mod(i,5)+1)=max(max(dat6));
  min1(mod(i,5)+1)=min(min(dat));
  min2(mod(i,5)+1)=min(min(dat2));
  min3(mod(i,5)+1)=min(min(dat3));
  min4(mod(i,5)+1)=min(min(dat4));
  min5(mod(i,5)+1)=min(min(dat5));
  min6(mod(i,5)+1)=min(min(dat6));
  
  iuse=find(min1~=0);
  
  norm1=mean(max1(iuse))-mean(min1(iuse));
  norm2=mean(max2(iuse))-mean(min2(iuse));
  norm3=mean(max3(iuse))-mean(min3(iuse));
  norm4=mean(max4(iuse))-mean(min4(iuse));
  norm5=mean(max5(iuse))-mean(min5(iuse));
  norm6=mean(max6(iuse))-mean(min6(iuse));
  shift1=mean(min1(iuse));
  shift2=mean(min2(iuse));
  shift3=mean(min3(iuse));
  shift4=mean(min4(iuse));
  shift5=mean(min5(iuse));    
  shift6=mean(min6(iuse));    
    
  %normalise onto [0 1]
  dat=fac*(dat-shift1)/norm1;
  dat2=fac*(dat2-shift2)/norm2;
  dat3=fac*(dat3-shift3)/norm3;
  dat4=fac*(dat4-shift4)/norm4;
  dat5=fac*(dat5-shift5)/norm5;  
  dat6=fac*(dat6-shift6)/norm6; 

  %avoid saturation of color tables (due to time average of color table)
  dat(find(dat<0))=shift; dat2(find(dat2<0))=shift; dat3(find(dat3<0))=shift;  dat4(find(dat4<0))=shift; dat5(find(dat5<0))=shift; dat6(find(dat6<0))=shift; 
  dat(find(dat>1))=1-shift; dat2(find(dat2>1))=1-shift; dat3(find(dat3>1))=1-shift;  dat4(find(dat4<0))=1-shift; dat5(find(dat5>1))=1-shift; dat6(find(dat6>1))=1-shift;

  % stack up contiguously
  dat2=(dat2+1);
  dat3=(dat3+2);
  dat4=(dat4+3);
  dat5=(dat5+4); 
  dat6=(dat6+5);
  
  %maxv=max(max(max(dat5)));
  %minv=min(min(min(dat)));
  
  %maxv=5;
  maxv=6;
  minv=0;  
  
  figure(h);

  axes(ha(1))
  %subplot(1,5,1)  
  contourf(xphi,yphi,dat,levels);   
  set(ha(1),'CLim',[minv maxv]);
  axis equal;
  shading flat;  
  axis([xmin xmax ymin ymax]);
  %title([file ' ' num2str(i)])
  title(['\Gamma_W   ' num2str(i)])
  %colormap(gca,'gray')
  %if (i==start) freezeColors; end
  %set(gca,'YTickLabel',[])
    
  axes(ha(2))
  %subplot(1,5,2)  
  contourf(xphi,yphi,dat2,levels*2); 
  set(ha(2),'CLim',[minv maxv]);
  axis equal;
  shading flat;  
  axis([xmin xmax ymin ymax]);
  %title([file2 ' ' num2str(i)])
  title(['\phi   ' num2str(i)])
  %colormap(gca,'summer')
  %if (i==start) freezeColors; end
  set(gca,'YTickLabel',[])
  
  axes(ha(3))
  %subplot(1,5,3)  
  contourf(xphi,yphi,dat3,levels*2);  
  set(ha(3),'CLim',[minv maxv]);
  axis equal;
  shading flat;  
  axis([xmin xmax ymin ymax]);
  %title([file3 ' ' num2str(i)])
  title(['n_e   ' num2str(i)])
  %colormap(gca,'bone')
  %if (i==start) freezeColors; end
  set(gca,'YTickLabel',[])
  
  axes(ha(4))
  %subplot(1,5,4)  
  contourf(xphi,yphi,dat4,levels);  
  set(ha(4),'CLim',[minv maxv]);
  axis equal;
  shading flat;  
  axis([xmin xmax ymin ymax]);
  %title([file4 ' ' num2str(i)])  
  xlabel(['T_e   ' num2str(i)])
  %colormap(gca,'copper')
  %if (i==start) freezeColors; end
  %set(gca,'YTickLabel',[])
  
  axes(ha(5))
  %subplot(1,5,5)  
  contourf(xphi,yphi,dat5,levels);
  set(ha(5),'CLim',[minv maxv]);
  axis equal;
  shading flat;  
  axis([xmin xmax ymin ymax]);
  %title([file5 ' ' num2str(i)])
  xlabel(['J_{||}   ' num2str(i)])
  %colormap(gca,'pink')
  %if (i==start) freezeColors; end
  set(gca,'YTickLabel',[])
  
  axes(ha(6))
  %subplot(1,5,5)  
  contourf(xphi,yphi,dat6,levels);
  set(ha(6),'CLim',[minv maxv]);
  axis equal;
  shading flat;  
  axis([xmin xmax ymin ymax]);
  %title([file5 ' ' num2str(i)])
  xlabel(['T_i ' num2str(i)])
  %colormap(gca,'pink')
  %if (i==start) freezeColors; end
  set(gca,'YTickLabel',[])
        
  %caxis([-0.8*maxv 0.8*maxv]);
   
  M(il) = getframe(h);
  %clear(name)
  
  if (mod(i,100)==0) save Movie_bck.mat M -v7.3; end
end;

% The version is needed to save data > 2GB (also compresses it)
save movie_out.mat M -v7.3

% % This bit to write the movie to a compressed avi.
% % VideoWriter must be used then the M array is greater than 2GB.
% % VideoWriter only in versions of matlab > 2010b
 clear OBJ
 OBJ = VideoWriter('movie_out.avi')
 OBJ.Quality=70;
 OBJ.FrameRate=10;
 open(OBJ)
% 
% % Write frame by frame (this mehtod can also be used as the frames are made)
% % This would avoid the need to store all of M in memory
%for ik=1:length(M)
%    writeVideo(OBJ,M(ik));
% %    disp(num2str(ik))
% %end
% 
% % Or write the whole movie in one go
 writeVideo(OBJ,M)
 close(OBJ);
