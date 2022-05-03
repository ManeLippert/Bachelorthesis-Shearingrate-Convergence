load xphi;
load yphi;
xmin = min(min(xphi));
xmax = max(max(xphi));
ymin = min(min(yphi));
ymax = max(max(yphi));
clear M;
close all; 
maxv = 0; 
minv = 0;
figure(1);
%colormap('bone')
figure(2);
hold on ;
file='phi'

for j= start:fin 
  i = j
  s = int2str(i);
  if (i<10) 
    name = strcat(file,'00000',s)
  else 
    if (i<100) 
      name = strcat(file,'0000',s)
    else
      if (i<1000) 
        name = strcat(file,'000',s)
      else
        if (i<10000) 
          name = strcat(file,'00',s)
        else
          if (i<100000) 
            name = strcat(file,'0',s)
          else
            name = strcat(file,s)
          end
        end
      end
    end
  end   
  dat=load(name);
  maxv=max(maxv,max(max(dat)));
  minv=min(minv,min(min(dat)));
  norm = sum(sum(dat));
  figure(2);
  plot([i],[norm],'bo');
end;
maxv = max(maxv,abs(minv));
vv = linspace(-1.1*maxv, 1.1*maxv, 20); 
figure(1); 

for j= start: fin
  i = j
  s = int2str(i);
  if (i<10) 
    name = strcat(file,'00000',s)
  else 
    if (i<100) 
      name = strcat(file,'0000',s)
    else
      if (i<1000) 
        name = strcat(file,'000',s)
      else
        if (i<10000) 
          name = strcat(file,'00',s)
        else
          if (i<100000) 
            name = strcat(file,'0',s)
          else
            name = strcat(file,s)
          end
        end
      end
    end
  end   
  dat=load(name);
  figure(1);
  shading flat;
  axis equal; 
  axis([xmin,xmax,ymin,ymax]); 
  contourf(xphi,yphi,dat,vv); 
  axis equal;
  axis([xmin xmax ymin ymax]);
  caxis([-0.8*maxv 0.8*maxv]);
  shading flat; 
  M(j) = getframe;
  clear(name)
end;
