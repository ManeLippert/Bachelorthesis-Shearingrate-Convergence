load xphi;
load yphi;
xmin = min(min(xphi));
xmax = max(max(xphi));
ymin = min(min(yphi));
ymax = max(max(yphi));
clear N;
close all; 
maxv = 0; 
minv = 0;
figure(1);
%colormap('bone')
figure(2);
hold on ;
for j= 1: nframes 
  i = j-1;
  s = int2str(i);
  if (i<10) 
    name = strcat('fdi00000',s)
  else 
    if (i<100) 
      name = strcat('fdi0000',s)
    else
      if (i<1000) 
        name = strcat('fdi000',s)
      else
        if (i<10000) 
          name = strcat('fdi00',s)
        else
          if (i<100000) 
            name = strcat('fdi0',s)
          else
            name = strcat('fdi',s)
          end
        end
      end
    end
  end   
  ll = ['load ' name ];
  eval(ll);
  dd = ['maxv=max(maxv,max(max(' name ')))'];
  eval(dd);
  dd = ['minv=min(minv,min(min(' name ')))'];
  eval(dd);
  dd = ['norm = sum(sum(' name '))'];
  eval(dd);
  %figure(2);
  %plot([i],[norm],'bo');
end;
maxv = max(maxv,abs(minv));
vv = linspace(-1.1*maxv, 1.1*maxv, 20); 
figure(1); 
load fdi000000
contourf(xphi,yphi,fdi000000,vv)
axis equal; 
shading flat; 
for j= 1: nframes
  i = j-1
  s = int2str(i);
  if (i<10) 
    name = strcat('fdi00000',s)
  else 
    if (i<100) 
      name = strcat('fdi0000',s)
    else
      if (i<1000) 
        name = strcat('fdi000',s)
      else
        if (i<10000) 
          name = strcat('fdi00',s)
        else
          if (i<100000) 
            name = strcat('fdi0',s)
          else
            name = strcat('fdi',s)
          end
        end
      end
    end
  end   
  ll = ['load ' name ];
  eval(ll);
  figure(1);
  axis equal; 
  axis([xmin,xmax,ymin,ymax]); 
  pp = strcat('contourf(xphi,yphi,',name,',20)'); 
  eval(pp);
  axis equal;
  %axis([xmin xmax ymin ymax]);
  %caxis([-0.8*maxv 0.8*maxv]);
  shading flat; 
  N(j) = getframe;
  clear(name)
end;
