clear;
xx=load('xphi');
yy=load('yphi');
%prof = load('prof_back.dat');
xysize=size(xx);
x1 = xysize(2)/2;
x2 = x1+1;
y1 = xysize(1)/2;
y2 = y1+1;
bmin = 0.84;
radp = xysize(2);
torp = xysize(1);
rhostar = 5.0E-3;
factor=1;
factor_r=2;
nmod = 1;

%Set to 1 for flux tube
is_flux_tube = 0;
%Set to 1 for 2nd order integration
second_order = 1;
%Suppresses the plot
poinc_plot = 0;
%Filters out the high frequency modes
filter = 0;
%Change this parameter to determine how many modes are NOT filtered out
%The zero mode is left
fil_leave = 1;
%Is the run data linear
is_linear = 0;
mn=1;

if(is_flux_tube==1)
    radp_geom=1;
else
    radp_geom=radp;
end

nturns = 300;

GEOM = read_geom('geom.dat');

ns = GEOM.ns;
ns_real = 2*ns*factor;

epsizet = reshape(GEOM.e_eps_zeta,ns,radp_geom);
epsis = reshape(GEOM.e_eps_s,ns,radp_geom);
ffun = reshape(GEOM.f,ns,radp_geom);
epsizet2(2*factor*ns,radp_geom)=0; 
epsis2(2*factor*ns,radp_geom) =0; 
ffun2(2*factor*ns,radp_geom) =0;
q = GEOM.q; 
s_grid = GEOM.s_grid;
ds = (GEOM.s_grid(ns)-GEOM.s_grid(ns-1))/factor;
dzeta = 1/torp; %2*pi/torp;%The 2pi is normalised out
s0 = ns/2;

s_grid_i(ns*2*factor)=0;
s_grid_i(ns*factor)=-ds/2;
s_grid_i(ns*factor+1)=ds/2;

for i = 1:1:(ns*factor-1)
    ii=ns*factor+1+i;
    s_grid_i(ii)=s_grid_i(ii-1)+ds/2; 
    s_grid_i(ns*factor-i)=-s_grid_i(ii);
end
s_grid_i=s_grid_i';
s_grid2 = interp(s_grid,factor);

if(is_flux_tube)
    dx = xx(1,2)-xx(1,1);
else
    dx = GEOM.eps(radp)-GEOM.eps(radp-1);
end


%Equally spaced inital points
psi0(1:xysize(2)) = GEOM.eps;
zeta0 = 0.5*ones([radp 1]);

fid2 = fopen('Apara00000010','r');
fid = fopen('Apara00000010','r');
a = fread(fid2,'double');
ph = fread(fid,'double');
phi3d = reshape(ph,[torp radp ns]);
%Reshape the data into 3d slab
apa3d = reshape(a,[torp radp ns]);
buff(ns)=0;
apa3d2(torp,radp,2*factor*ns)=0;
cc(torp,radp,2*factor*ns)=0;dd(torp,radp,2*factor*ns)=0;
rat(radp,2*factor*ns)=0;
%s_grid=interp(s_grid,factor);

if(torp<128)
    apa3d3(128,radp,ns)=0.0;
    if(is_linear)
      ft2(1,radp)=0;
      for i=1:ns
          ft = fft(10*apa3d(:,:,i),[],1);
          ft2(1+mn,:)=ft(1,:);
          ft2(128-(mn-1),:)=ft(1,:);
          apa3d3(:,:,i)=ifft(ft2,[],1);
      end
    else
      for i=1:ns
        ft = fft(apa3d(:,:,i),[],1);
        ft2(1,:) = ft(1,:);
        for j=1:(nmod-1)
          ft2(1+j,:)=ft(1+j,:);
          ft2(128+1-j,:)=ft(torp+1-j,:);
        end   
        apa3d3(:,:,i)=ifft(ft2,[],1);
      end
    end
    
   torp = 128;
   dzeta = 1/torp;
   apa3d2(torp,radp,2*factor*ns)=0;
   cc(torp,radp,2*factor*ns)=0;dd(torp,radp,2*factor*ns)=0;
   clear apa3d
   apa3d = apa3d3;
   clear apa3d3
end

if(filter)
    for i=1:ns
      ft = fft(apa3d(:,:,i),[],1);
      ft((2+fil_leave):(torp-fil_leave),:)=0.0;
      ft(1,:)=0.0;
      apa3d(:,:,i)=ifft(ft,[],1);
    end
    clear ft
    %This filters out high frequency stuff along the field line
    for i=1:torp
      dum(:,:) = apa3d(i,:,:);
      ft = fft(dum,[],2);
      ft(:,(ns/2-1):(ns/2+3))=0.0;
      apa3d(i,:,:)=ifft(ft,[],2);
    end
end      


%Interpolate to make the data finer in the s direction
for j=1:radp
    for i=1:torp
        buff(:) = apa3d(i,j,:);
        buff2(1:factor*ns*2) = interp1(s_grid,buff,s_grid_i,'spline');
        apa3d2(i,j,:)=buff2;
    end
end



for j=1:radp_geom
    buff(:) = ffun(:,j);
    buff2(1:factor*ns*2) = interp1(s_grid,buff,s_grid_i,'spline');
    ffun2(:,j) = buff2;
    
    buff(:) = epsizet(:,j);
    buff2(1:factor*ns*2) = interp1(s_grid,buff,s_grid_i,'spline');
    epsizet2(:,j) = buff2;
    
    buff(:) = epsis(:,j);
    buff2(1:factor*ns*2) = interp1(s_grid,buff,s_grid_i,'spline');
    epsis2(:,j) = buff2;
end

%Pre-calculate the gradients of Aparallel and the ratio of
%tensors
for i=1:2*factor*ns
  rat(:,i) = epsizet2(i,:)./ffun2(i,:);   
  aa = apa3d2(:,:,i);     
  %Gradient in radial coordinate
  cc(:,:,i) = gradient(aa)/dx;
  %Gradient in zeta
  dd(:,:,i) = gradient(aa')'/dzeta;
end
   
%cc(:,:,:)=0;
%dd(:,:,:)=0;
%This is the initial condition of the field lines. Equally spaced in radius
%and at a constant zeta
psi0 = psi0';
psi=interp(psi0,factor_r);
zeta=interp(zeta0,factor_r);

psi(1:4) = 0.15;

zeta(1) = 0.2;
zeta(2) = 0.7;
zeta(3) = 0.8;
zeta(4) = 0.3;

num_field_lines = size(psi);
num_field_lines = num_field_lines(1)
%Arrays holding the coordinates after each turn
zeta_n(num_field_lines,nturns)=0;
psi_n(num_field_lines,nturns)=0;

coord_store{factor*ns*nturns} = 0;
coord_store_cart{factor*ns*nturns} = 0;

if(poinc_plot==1)
  figure(7);plot(psi,zeta,'r.','MarkerSize',2);
end

%Coordinates for interpolation
[xcuns ycuns]=meshgrid(psi0,[1:1:torp]/torp);

%The integral part
for kk=1:nturns
  for i=factor*s0+1:factor*ns
  %for i=1:factor*ns 
  
      ii = 2*i - 1;
      ii_mid = ii+1;
      
      %Check for NaNs
      psi(isnan(psi)) = 0;
      zeta(isnan(zeta)) = 0;
      
      %Must put radial interpolation here somewhere
      %Radial interpolation of ratio of metric elements (they are invariant
      %in zeta so only 1D is needed
      rat_local(1:radp) = rat(1:radp,ii);
      rat_interp = interp1(psi0,rat_local,psi,'spline'); 
      
      %2D interpolation in the plane
      dd_interp = interp2(xcuns,ycuns,dd(:,:,ii),psi,zeta,'spline');
      cc_interp = interp2(xcuns,ycuns,cc(:,:,ii),psi,zeta,'spline');
      %Dirty hack
      dd_interp(isnan(dd_interp))=0; 
      cc_interp(isnan(cc_interp))=0;
            
      if(second_order==1)
        k1_psi = -rhostar*rhostar*(2*rat_interp).*dd_interp*ds/2;
        k1_zeta = rhostar*rhostar*(2*rat_interp).*cc_interp*ds/2;
 
        k1_psi = psi + k1_psi;
        k1_zeta = zeta + k1_zeta;
        
        k1_zeta(gt(k1_zeta,1.0)) = k1_zeta(gt(k1_zeta,1.0)) - 1;
        k1_zeta(lt(k1_zeta,0.0)) = k1_zeta(lt(k1_zeta,0.0)) + 1;
      
        rat_local(1:radp) = rat(1:radp,ii_mid);
        rat_interp = interp1(psi0,rat_local,psi,'spline'); 

        %2D interpolation in the plane
        dd_interp = interp2(xcuns,ycuns,dd(:,:,ii_mid),k1_psi,k1_zeta,'spline');
        cc_interp = interp2(xcuns,ycuns,cc(:,:,ii_mid),k1_psi,k1_zeta,'spline');

        %Dirty hack
        dd_interp(isnan(dd_interp))=0; 
        cc_interp(isnan(cc_interp))=0; 
      end
      
      psi = psi - rhostar*rhostar*(2*rat_interp).*dd_interp*ds;
      zeta = zeta + rhostar*rhostar*(2*rat_interp).*cc_interp*ds;

      psi(gt(psi,max(psi0))) = 0;
      psi(lt(psi,min(psi0))) = 0;
      
      zeta(gt(zeta,1.0)) = zeta(gt(zeta,1.0)) - 1;
      zeta(lt(zeta,0.0)) = zeta(lt(zeta,0.0)) + 1;

      coord_store{(kk-1)*factor*ns+i} = [psi(5) zeta(5) 2*pi*s_grid2(i)]; 
      q_interp = interp1(psi0,q,psi(5),'spline');
      phi_coord = 2*pi*s_grid2(i)*q_interp - 2*pi*zeta(5);
      x = (3+psi(5)*cos(2*pi*s_grid2(i)))*cos(phi_coord);
      y = (3+psi(5)*cos(2*pi*s_grid2(i)))*sin(phi_coord);
      z = psi(5)*sin(2*pi*s_grid2(i));
      coord_store_cart{(kk-1)*factor*ns+i} = [x y z];
  end
  
  %Boundary condition shift must go here it seems. The shift in zeta
  %is q
  kk
  q_interp = interp1(psi0,q,psi,'spline');
  zeta = zeta + (q_interp-floor(q_interp));
  for cu=1:num_field_lines
    if(zeta(cu)>1)
       zeta(cu)=zeta(cu)-1;
    end
  end
  psi_n(:,kk) = psi;
  zeta_n(:,kk) = zeta;
  
  if(poinc_plot==1)
    figure(7);hold on;plot(psi,zeta,'k.','MarkerSize',2)
    %figure(2);hold on;scatter(psi,zeta,1,interp(psi0,factor_r));
  end
  
  for i=1:factor*s0
      
      ii = 2*i - 1;
      ii_mid = ii+1;
      
      %Check for NaNs
      psi(isnan(psi)) = 0;
      zeta(isnan(zeta)) = 0;
      
      %Must put radial interpolation here somewhere
      %Radial interpolation of ratio of metric elements (they are invariant
      %in zeta so only 1D is needed
      rat_local(1:radp) = rat(1:radp,i);
      rat_interp = interp1(psi0,rat_local,psi,'spline'); 
      
      %2D interpolation in the plane
      dd_interp = interp2(xcuns,ycuns,dd(:,:,ii),psi,zeta,'spline');
      cc_interp = interp2(xcuns,ycuns,cc(:,:,ii),psi,zeta,'spline');
      %Dirty hack
      dd_interp(isnan(dd_interp))=0; 
      cc_interp(isnan(cc_interp))=0;
      
      if(second_order==1)
        k1_psi = -rhostar*rhostar*(2*rat_interp).*dd_interp*ds/2;
        k1_zeta = rhostar*rhostar*(2*rat_interp).*cc_interp*ds/2;

        k1_psi = psi + k1_psi;
        k1_zeta = zeta + k1_zeta;

        k1_zeta(gt(k1_zeta,1.0)) = k1_zeta(gt(k1_zeta,1.0)) - 1;
        k1_zeta(lt(k1_zeta,0.0)) = k1_zeta(lt(k1_zeta,0.0)) + 1;
      
        rat_local(1:radp) = rat(1:radp,ii_mid);
        rat_interp = interp1(psi0,rat_local,psi,'spline'); 
        %2D interpolation in the plane
        dd_interp = interp2(xcuns,ycuns,dd(:,:,ii_mid),k1_psi,k1_zeta,'spline');
        cc_interp = interp2(xcuns,ycuns,cc(:,:,ii_mid),k1_psi,k1_zeta,'spline');
        %Dirty hack
        dd_interp(isnan(dd_interp))=0; 
        cc_interp(isnan(cc_interp))=0; 
      end
      
      psi = psi - rhostar*rhostar*(2*rat_interp).*dd_interp*ds;
      zeta = zeta + rhostar*rhostar*(2*rat_interp).*cc_interp*ds;
      
      psi(gt(psi,max(psi0))) = 0;
      psi(lt(psi,min(psi0))) = 0;
      
      zeta(gt(zeta,1.0)) = zeta(gt(zeta,1.0)) - 1;
      zeta(lt(zeta,0.0)) = zeta(lt(zeta,0.0)) + 1;
      
      coord_store{(kk-1)*factor*ns+i} = [psi(5) zeta(5) 2*pi*s_grid2(i)]; 
      q_interp = interp1(psi0,q,psi(5),'spline');
      phi_coord = 2*pi*s_grid2(i)*q_interp - 2*pi*zeta(5);
      x = (3+psi(5)*cos(2*pi*s_grid2(i)))*cos(phi_coord);
      y = (3+psi(5)*cos(2*pi*s_grid2(i)))*sin(phi_coord);
      z = psi(5)*sin(2*pi*s_grid2(i));
      coord_store_cart{(kk-1)*factor*ns+i} = [x y z];
  end
end

save('Poincare1.mat','psi_n','zeta_n','coord_store_cart');