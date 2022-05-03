%*******************************************************************************************************************
%This programm integrates and thereby traces fieldlines in a tokamak reactor.
%It creates a Poincare plot in the outboard midplane and stores the 3d cartesian coordinates of a specific fieldline.
%It was written by Benjamin Schobert
%*******************************************************************************************************************

%Inputs:        xphi
%               geom.dat
%               Apara00000***
%               kzeta                   (only for local runs)
%               lxn                     (only for local runs)

%Outputs:       Poincare1.mat           contains:
%               psi_n                   all psi coordinates for each turn
%               zeta_n                  all zeta coordinates for each turn
%               coord_store_cart        the cartesian coordinates for a specific field plot
%               figure(1)               the Poincare plot in the outboard midplane

clear;                                  %clear memory
clf(figure(1));                         %clear figure(1) Poincare plot
clf(figure(3));                         %clear psimaxdisp plot

rhostar = 5e-3;                         %normalised gyro radius (only for global runs)
factor_s=1;                             %refinement factor in s direction
factor_r=2;                             %refinement factor in psi and zeta direction
nmod = 21;                              %number of modes
mn=1;                                   %only for linear runs
nturns = 100;                           %number of turns
fl_3d = 5;                              %this is the number of the field line, whose 3d coordinates will be saved
dotsize = 2;                            %this is the dotsize in the Poincare plot. I recommend 2 for octave plots
                                        %and 0.5 for png oder pdf or eps
nfields = 1;                            %number of different A_para fields that are used (only for time dependency)
napara = 200;                           %number of the Apara file (Apara00000***)
                                        
%Flags
%Set to 1 for flux tube or 0 for global run
is_flux_tube = 1;

%Filters out the high frequency modes if set to 1
filter = 0;

%Change this parameter to determine how many modes are NOT filtered out
%The zero mode is left
fil_leave = 1;

%Set to 1 if run data is linear or 0 if not
is_linear = 0;

%Set to 1 for additional plot of the diffusion coefficient
diff_plot = 1;

%Set to 1 for additional plot of the maximal displacement in psidirection after half a turn
maxdisp_plot = 1;

%Set to 1 for saving all plots as eps files with latex code
save_plots = 0;

%Switch for initial condition
%IC1 is a line where zeta = 0.5 and the field lines are equally spaced in psi.
%IC2 is a grid where the field lines are seeded at zeta = 1/(2*ic2zetasize):1/ic2zetasize:1-1/(2*ic2zetasize) and eqally spaced in psi.

ic = 1;
ic2zetasize = 20;                       %this is the zeta ic gridsize (only relevant for ic = 2)


%**********************************************************************************************************
%The code
%**********************************************************************************************************

function [dx, dy] = pbcgrad(m)                  %m is a 2d matrix, dx is the gradient of the colums, dy the gradient of the rows
  nrows = size(m)(1);                           %it calculates the derivative with a five point stencil and pbc in both directions
  ncol = size(m)(2);
  dx(nrows,ncol) = 0;
  dy(nrows,ncol) = 0;
    for j = 1:ncol
      jm1 = j-1;
      jm2 = j-2;
      jp1 = j+1;
      jp2 = j+2;
      if(j==1)
        jm1 = ncol;
        jm2 = ncol-1;
      endif
      if(j==2)
        jm2 = ncol;
      endif
      if(j==ncol-1)
        jp2=1;
      endif
      if(j==ncol)
        jp1=1;
        jp2=2;
      endif
      dx(:,j)= (-m(:,jp2)+8*m(:,jp1)-8*m(:,jm1)+m(:,jm2))/(12);
    endfor
  
    for i = 1:nrows
      im1 = i-1;
      im2 = i-2;
      ip1 = i+1;
      ip2 = i+2;
      if(i==1)
        im1 = nrows;
        im2 = nrows-1;
      endif
      if(i==2)
        im2 = nrows;
      endif
      if(i==nrows-1)
        ip2=1;
      endif
      if(i==nrows)
        ip1=1;
        ip2=2;
      endif
      dy(i,:)= (-m(ip2,:)+8*m(ip1,:)-8*m(im1,:)+m(im2,:))/(12);
    endfor
endfunction


if(is_flux_tube==1)
    radp_geom=1;                                        %for flux tube is the radial grid size 1
    lx = load('lxn');                                   %radial boxsize in psi direction
    kzeta = load('kzeta');                              %kzeta stores the wave vectors in zeta direction multiplied with rhostar
    kzeta_min = kzeta(2,1);                             %this returns the minimal non zero wave vector
    int_factor = kzeta_min/(2*pi);                      %this is the factor needed for the numerical integration
    ly = 1/int_factor;
else
    radp_geom=radp;
    int_factor = rhostar*rhostar;
end

GEOM = read_geom('geom.dat');                           %reads geom.dat
ns = GEOM.ns;                                           %ns is the gridsize in s direction (parallel B)
eps = GEOM.eps;                                         %eps is the ratio r/R
shat = GEOM.shat;                                       %shat is the electromagnetic shear
xphi=load('xphi');                                      %loads xphi into matrix
xphi_size=size(xphi);                                   %stores size of xphi in form [rows, colums]
radp = xphi_size(2);                                    %number of colums = radial gridsize
torp = xphi_size(1);                                    %number of rows = toridal gridsize


epsizet = reshape(GEOM.e_eps_zeta,ns,radp_geom);                %generates matrix from e_eps_zeta with ns rows and radp_geom colums
ffun = reshape(GEOM.f,ns,radp_geom);                            %generates matrix from F with ns rows and radp_geom colums
epsizet2(2*factor_s*ns,radp_geom)=0;                            %generates matrix filled with zeros, double number of rows
ffun2(2*factor_s*ns,radp_geom) =0;
q = GEOM.q;                                                     %the savety factor q
s_grid = GEOM.s_grid;
ds = 1/(factor_s*ns);                                           %stepwidth in s direction
dzeta = 1/torp;                                                 %2*pi/torp;     %The 2pi is normalised out
s0 = ns/2;

s_grid_i(ns*2*factor_s)=0;                                      %generates a vector with zeros
s_grid_i(ns*factor_s)=-ds/2;                                    %writes middle entries in the vector
s_grid_i(ns*factor_s+1)=ds/2;

for i = 1:1:(ns*factor_s-1)
    ii=ns*factor_s+1+i;
    s_grid_i(ii)=s_grid_i(ii-1)+ds/2;                           %fills second half of the vector with values increasing by ds/2
    s_grid_i(ns*factor_s-i)=-s_grid_i(ii);                      %fills lower half of the vector result (-3,-2,-1,1,2,3)
end

s_grid2 = (-0.5+(ds/2)):ds:(0.5-(ds/2));

s_grid_i=s_grid_i';                                             %transposes s_grid vector
s_grid2=s_grid2';


if(is_flux_tube)                                                %calculate stepwidth dpsi
    dpsi = lx/radp;
    shift_factor = (kzeta_min*q*shat)/(2*pi*eps);               %this is the factor for the zeta displacement due to the shear and the local coordinates
else
    dpsi = eps(radp)-eps(radp-1);               
end


%This is the initial condition of the field lines. 
%It depends on the variable ic, ic = 1 is the default case.

if (is_flux_tube)
  psi0 = -lx/2+dpsi/2:dpsi:lx/2-dpsi/2;               %generates psi0 vector of length lx with radp equidistant entries
else
  psi0(1:radp) = eps;
endif

zeta0 = 0.5*ones([radp 1]);                                             %ones([n 1]) returns a vector length n filled with ones, so every entry will be 0.5

psi_refined_grid = 1/(2*factor_r):1/factor_r:radp-1/(2*factor_r);       %interpolation in psi and zeta direction
zeta_refined_grid = 1/(2*factor_r):1/factor_r:radp-1/(2*factor_r);
psi = interp1(psi0, psi_refined_grid,'spline','extrap');
zeta = interp1(zeta0, zeta_refined_grid,'spline','extrap');

switch(ic)
  case 1
    
  case 2
      psi_refined_grid = 1/(2*factor_r*ic2zetasize):1/(factor_r*ic2zetasize):radp-1/(2*factor_r*ic2zetasize);
      psi = interp1(psi0, psi_refined_grid,'spline','extrap');
      zeta0 = 1/(2*ic2zetasize):1/ic2zetasize:1-1/(2*ic2zetasize);
      zeta = repmat(zeta0,1,factor_r*radp);
  otherwise
  
endswitch

psi0 = psi0';                                           %psi0 gets transposed

if(is_flux_tube)
  num_field_lines = size(psi);
  num_field_lines = num_field_lines(2);
else
  num_field_lines = size(psi);
  num_field_lines = num_field_lines(1);
end

%Arrays holding the coordinates after each turn
zeta_n(num_field_lines,nturns*nfields)=0;
psi_n(num_field_lines,nturns*nfields)=0;
diff_psi_n(num_field_lines,nturns*nfields)=0;
diff_zeta_n(num_field_lines,nturns*nfields)=0;
diff_psi(num_field_lines)=0;
diff_zeta(num_field_lines)=0;
deltapsi(num_field_lines)=0;
deltazeta(num_field_lines)=0;
halfpsi(num_field_lines)=0;
psimaxdisp(2*nturns*nfields)=0;

coord_store{factor_s*ns*nturns*nfields} = 0;            %stores coordinates of a specific fieldline at each step in each turn (psi, zeta, f(s))
coord_store_cart{factor_s*ns*nturns*nfields} = 0;       %stores the cartesian coordinates of the same fieldline


for ntime = 0:1:nfields-1                               %this is the time loop, for changing the field

if((napara+ntime)<100)
  fid = fopen(['Apara000000' int2str(napara+ntime)],'r');  %open Apara00000*** for reading
else
  fid = fopen(['Apara00000' int2str(napara+ntime)],'r');
endif
apara = fread(fid,'double');                            %read binary with precision double
apa3d = reshape(apara,[torp radp ns]);                  %Reshape the data into 3d slab
buff(ns)=0;                                             %new vector buff
apa3d2(torp,radp,2*factor_s*ns)=0;                      %new 3d matrix
cc(torp,radp,2*factor_s*ns)=0;                          %later needed for gradient in radial direction
dd(torp,radp,2*factor_s*ns)=0;                          %later needed for gradient in toroidal direction
rat(radp_geom,2*factor_s*ns)=0;                         %later needed for ratio of efun and ffun

if(torp<128)                                            %if the toridal gridsize is below 128 it will be refined
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
        ft = fft(apa3d(:,:,i),[],1);                    %fast fourier transform the 2d matrix with third index i, [] non specified size, 1 operate along colums
        ft2(1,:) = ft(1,:);                             
        for j=1:(nmod-1)                                
          ft2(1+j,:)=ft(1+j,:);
          ft2(128+1-j,:)=ft(torp+1-j,:);
        end   
        apa3d3(:,:,i)=ifft(ft2,[],1);
      end
    end
    
   torp_new = 128;
   dzeta = 1/torp_new;
   apa3d2(torp_new,radp,2*factor_s*ns)=0;
   cc(torp_new+2,radp+2,2*factor_s*ns)=0;
   dd(torp_new+2,radp+2,2*factor_s*ns)=0;
   clear apa3d;
   apa3d = apa3d3;
   clear apa3d3;
else
  torp_new = torp;
endif


if(filter)
    for i=1:ns
      ft = fft(apa3d(:,:,i),[],1);
      ft((2+fil_leave):(torp_new-fil_leave),:)=0.0;
      ft(1,:)=0.0;
      apa3d(:,:,i)=ifft(ft,[],1);
    end
    clear ft
    
    %This filters out high frequency modes along the field line
    for i=1:torp_new
      dum(:,:) = apa3d(i,:,:);
      ft = fft(dum,[],2);
      ft(:,(ns/2-1):(ns/2+3))=0.0;
      apa3d(i,:,:)=ifft(ft,[],2);
    end
end      


%Interpolate to make the data finer in the s direction
for j=1:radp
    for i=1:torp_new
        buff(:) = apa3d(i,j,:);
        buff2(1:factor_s*ns*2) = interp1(s_grid,buff,s_grid_i,'spline','extrap');
        apa3d2(i,j,:)=buff2;
    end
end

for j=1:radp_geom

    buff(:) = ffun(:,j);
    buff2(1:factor_s*ns*2) = interp1(s_grid,buff,s_grid_i,'spline','extrap');
    ffun2(:,j) = buff2;
    
    buff(:) = epsizet(:,j);
    buff2(1:factor_s*ns*2) = interp1(s_grid,buff,s_grid_i,'spline','extrap');
    epsizet2(:,j) = buff2;
    
end


%Pre-calculate the gradients of Aparallel and the ratio of tensors
for i=1:2*factor_s*ns
  rat(:,i) = epsizet2(i,:)./ffun2(i,:);
  %Gradients in both directions
  [dx,dy] = pbcgrad(apa3d2(:,:,i));
  dx = [dx(:,radp-1),dx,dx(:,1)];                       %periodically enlarge the matrices one row in each dimension
  dy = [dy(:,radp-1),dy,dy(:,1)];
  dx = [dx(torp_new-1,:);dx;dx(1,:)];
  dy = [dy(torp_new-1,:);dy;dy(1,:)];
  cc(:,:,i) = dx/dpsi;
  dd(:,:,i) = dy/dzeta;
end


%The Poincare plot
    figure(1);
    hold on;
    xlabel('$\psi_N$');
    ylabel('$\zeta_N$');
    axis([psi0(1)-dpsi/2 psi0(radp)+dpsi/2 0 1]);                         %set axis ranges
    grid;
    plot(psi(1:length(psi)/6),zeta(1:length(zeta)/6),'k.','MarkerSize',dotsize);                                                                    %black dots
    plot(psi(1+1*floor(length(psi)/6):2*floor(length(psi)/6)),zeta(1+1*floor(length(zeta)/6):2*floor(length(zeta)/6)),'b.','MarkerSize',dotsize);   %blue dots
    plot(psi(1+2*floor(length(psi)/6):3*floor(length(psi)/6)),zeta(1+2*floor(length(zeta)/6):3*floor(length(zeta)/6)),'g.','MarkerSize',dotsize);   %green dots
    plot(psi(1+3*floor(length(psi)/6):4*floor(length(psi)/6)),zeta(1+3*floor(length(zeta)/6):4*floor(length(zeta)/6)),'r.','MarkerSize',dotsize);   %red dots
    plot(psi(1+4*floor(length(psi)/6):5*floor(length(psi)/6)),zeta(1+4*floor(length(zeta)/6):5*floor(length(zeta)/6)),'c.','MarkerSize',dotsize);   %cyan dots
    plot(psi(1+5*floor(length(psi)/6):length(psi)),zeta(1+5*floor(length(zeta)/6):length(zeta)),'m.','MarkerSize',dotsize);                         %magenta dots


%Coordinates for interpolation
[xcuns ycuns]=meshgrid([psi0(1)-dpsi:dpsi:psi0(radp)+dpsi],[-dzeta/2:dzeta:1+dzeta/2]);


%The integral part
for kk=0:1:nturns-1
  for i=factor_s*s0+1:factor_s*ns                                               %starting at the outboard midplane and going anticlockwise
  
      ii = 2*i - 1;                                                             %the odd gridpoints are for first order integration
      ii_mid = ii+1;                                                            %the even gridpionts for second order
      
      
      %Radial interpolation of ratio of metric elements (they are invariant in zeta so only 1D is needed)
      if(is_flux_tube)
        rat_local = rat(ii);                                                    %only one rat vector if local
        rat_interp = rat_local;                                                 %rat_local is scalar if local
      else
        rat_local(1:radp)= rat(1:radp,ii);
        rat_interp = interp1(psi0,rat_local,psi,'spline','extrap'); 
      end
      
      %2D interpolation in the plane
      dd_interp = interp2(xcuns,ycuns,dd(:,:,ii),psi,zeta,'linear',0);           %values that need extrapolation will be set 0
      cc_interp = interp2(xcuns,ycuns,cc(:,:,ii),psi,zeta,'linear',0);           %creates vector of length 2*radp, psi and zeta are treated as the coordinates of scattered points
            
      %second order integration
        k1_psi = -int_factor*(2*rat_interp).*dd_interp*ds/2;
        k1_zeta = int_factor*(2*rat_interp).*cc_interp*ds/2;
 
        k1_psi = psi + k1_psi;
        k1_zeta = zeta + k1_zeta;
        
        k1_zeta(gt(k1_zeta,1.0)) = k1_zeta(gt(k1_zeta,1.0)) - 1;
        k1_zeta(lt(k1_zeta,0.0)) = k1_zeta(lt(k1_zeta,0.0)) + 1;
      
        if(is_flux_tube)
          rat_local=rat(ii_mid);
          rat_interp = rat_local;
        else
        rat_local(1:radp) = rat(1:radp,ii_mid);
        rat_interp = interp1(psi0,rat_local,psi,'spline','extrap');
        end

        %2D interpolation in the plane
        dd_interp = interp2(xcuns,ycuns,dd(:,:,ii_mid),k1_psi,k1_zeta,'linear',0);
        cc_interp = interp2(xcuns,ycuns,cc(:,:,ii_mid),k1_psi,k1_zeta,'linear',0);       
      %end of second order integration
     
      deltapsi = - int_factor*(2*rat_interp).*dd_interp*ds;
      deltazeta = + int_factor*(2*rat_interp).*cc_interp*ds;
      psi += deltapsi;
      halfpsi += deltapsi;
      zeta += deltazeta;
      diff_psi += deltapsi;
      diff_zeta += deltazeta;

      if(is_flux_tube)
        psi(gt(psi,lx/2)) = psi(gt(psi,lx/2)) - lx;                             %periodic boundary condition for psi
        psi(lt(psi,-lx/2)) = psi(lt(psi,-lx/2)) + lx;
      else
        psi(gt(psi,max(psi0))) = 0;                                             %treatment for out of bound psi in global
        psi(lt(psi,min(psi0))) = 0;
      end
      
      zeta(gt(zeta,1.0)) = zeta(gt(zeta,1.0)) - 1;                              %periodic boundary condition for zeta
      zeta(lt(zeta,0.0)) = zeta(lt(zeta,0.0)) + 1;

      coord_store(kk*factor_s*ns+i+ntime*nturns) = [psi(fl_3d) zeta(fl_3d) 2*pi*s_grid2(i)]; 
      if(is_flux_tube==0)
        q_interp = interp1(psi0,q,psi,'spline','extrap');
        phi_coord = 2*pi*s_grid2(i)*q_interp - 2*pi*zeta(fl_3d);
        x = (3+psi(fl_3d)*cos(2*pi*s_grid2(i)))*cos(phi_coord);
        y = (3+psi(fl_3d)*cos(2*pi*s_grid2(i)))*sin(phi_coord);
        z = psi(fl_3d)*sin(2*pi*s_grid2(i));
        coord_store_cart((kk-1)*factor_s*ns+i+ntime*nturns) = [x y z];
      end
  end
 
  psimaxdisp((2*kk)+1+(2*ntime*nturns))= max(abs(halfpsi));                       %This saves the most displaced field line (in psi) after half a turn
  halfpsi = 0;
  
  %Boundary condition shift in the inboard midplane
  if(is_flux_tube)
    for jj=1:num_field_lines
      zeta(jj) += shift_factor*psi(jj);                         %The shift in zeta generated by the shearing of the coordinates
      if(zeta(jj)>1)
        zeta(jj)=zeta(jj)-floor(zeta(jj));
      endif
      if(zeta(jj)<0)
        zeta(jj)=zeta(jj)-floor(zeta(jj));
      endif
    endfor
  else
    q_interp = interp1(psi0,q,psi,'spline','extrap');           %The shift in zeta is q for global runs
    zeta += q_interp-floor(q_interp);
    for jj=1:num_field_lines
      if(zeta(jj)>1)
        zeta(jj)=zeta(jj)-1;
      endif
    endfor
  endif
  
  
  for i=1:factor_s*s0                                         %starting in the inboard midplane and going anticlockwise
      
      ii = 2*i - 1;
      ii_mid = ii+1;
      
      
      %Radial interpolation of ratio of metric elements (they are invariant in zeta so only 1D is needed)
      if(is_flux_tube)
        rat_local = rat(i);                                            
        rat_interp = rat_local;
      else
        rat_local(1:radp)= rat(1:radp,i);
        rat_interp = interp1(psi0,rat_local,psi,'spline','extrap'); 
      end
      
      
      %2D interpolation in the plane
      dd_interp = interp2(xcuns,ycuns,dd(:,:,ii),psi,zeta,'linear',0);
      cc_interp = interp2(xcuns,ycuns,cc(:,:,ii),psi,zeta,'linear',0);
      
      
      %second order integration
        k1_psi = -int_factor*(2*rat_interp).*dd_interp*ds/2;
        k1_zeta = int_factor*(2*rat_interp).*cc_interp*ds/2;

        k1_psi = psi + k1_psi;
        k1_zeta = zeta + k1_zeta;

        k1_zeta(gt(k1_zeta,1.0)) = k1_zeta(gt(k1_zeta,1.0)) - 1;
        k1_zeta(lt(k1_zeta,0.0)) = k1_zeta(lt(k1_zeta,0.0)) + 1;
      
        if(is_flux_tube)
          rat_local=rat(ii_mid);
          rat_interp = rat_local;
        else
        rat_local(1:radp) = rat(1:radp,ii_mid);
        rat_interp = interp1(psi0,rat_local,psi,'spline','extrap');
        end
        
        %2D interpolation in the plane
        dd_interp = interp2(xcuns,ycuns,dd(:,:,ii_mid),k1_psi,k1_zeta,'linear',0);
        cc_interp = interp2(xcuns,ycuns,cc(:,:,ii_mid),k1_psi,k1_zeta,'linear',0);
      %end of second order integration
     
      deltapsi = - int_factor*(2*rat_interp).*dd_interp*ds;
      deltazeta = + int_factor*(2*rat_interp).*cc_interp*ds;
      psi += deltapsi;
      halfpsi += deltapsi;
      zeta += deltazeta;
      diff_psi += deltapsi;
      diff_zeta += deltazeta;
            
      if(is_flux_tube)
        psi(gt(psi,lx/2)) = psi(gt(psi,lx/2)) - lx;                             %periodic boundary condition for psi
        psi(lt(psi,-lx/2)) = psi(lt(psi,-lx/2)) + lx;
      else
        psi(gt(psi,max(psi0))) = 0;                                             %treatment for out of bound psi in global
        psi(lt(psi,min(psi0))) = 0;
      end
      
      zeta(gt(zeta,1.0)) = zeta(gt(zeta,1.0)) - 1;                              %periodic boundary condition for zeta
      zeta(lt(zeta,0.0)) = zeta(lt(zeta,0.0)) + 1;

      
      coord_store{kk*factor_s*ns+i+ntime*nturns} = [psi(fl_3d) zeta(fl_3d) 2*pi*s_grid2(i)]; 
      if(is_flux_tube==0)
        q_interp = interp1(psi0,q,psi,'spline','extrap');
        phi_coord = 2*pi*s_grid2(i)*q_interp - 2*pi*zeta(fl_3d);
        x = (3+psi(fl_3d)*cos(2*pi*s_grid2(i)))*cos(phi_coord);
        y = (3+psi(fl_3d)*cos(2*pi*s_grid2(i)))*sin(phi_coord);
        z = psi(fl_3d)*sin(2*pi*s_grid2(i));
        coord_store_cart{kk*factor_s*ns+i+ntime*nturns} = [x y z];
      end
  end
  
  kk                                                            %output current turn after each turn
  psimaxdisp((2*kk)+2+(2*ntime*nturns))= max(abs(halfpsi));         %This saves the most displaced field line (in psi) after half a turn
  halfpsi = 0;
  
  %This generates the plot in the outboard midplane
    figure(1);
    hold on;
    plot(psi(1:length(psi)/6),zeta(1:length(zeta)/6),'k.','MarkerSize',dotsize);                                                                    %black dots
    plot(psi(1+1*floor(length(psi)/6):2*floor(length(psi)/6)),zeta(1+1*floor(length(zeta)/6):2*floor(length(zeta)/6)),'b.','MarkerSize',dotsize);   %blue dots
    plot(psi(1+2*floor(length(psi)/6):3*floor(length(psi)/6)),zeta(1+2*floor(length(zeta)/6):3*floor(length(zeta)/6)),'g.','MarkerSize',dotsize);   %green dots
    plot(psi(1+3*floor(length(psi)/6):4*floor(length(psi)/6)),zeta(1+3*floor(length(zeta)/6):4*floor(length(zeta)/6)),'r.','MarkerSize',dotsize);   %red dots
    plot(psi(1+4*floor(length(psi)/6):5*floor(length(psi)/6)),zeta(1+4*floor(length(zeta)/6):5*floor(length(zeta)/6)),'c.','MarkerSize',dotsize);   %cyan dots
    plot(psi(1+5*floor(length(psi)/6):length(psi)),zeta(1+5*floor(length(zeta)/6):length(zeta)),'m.','MarkerSize',dotsize);                         %magenta dots
    %plot(psi(55),zeta(55),'rx','MarkerSize',2);                %this is for plotting a single FL and number each turn
    %if(kk>120)
    %number = num2str(kk);
    %text(psi(55)+0.2,zeta(55),number);
    %endif
  
  %Save the current coordinates
    psi_n(:,kk+1+ntime*nturns) = psi;                                            %store psi and zeta after each turn
    zeta_n(:,kk+1+ntime*nturns) = zeta;
    diff_psi_n(:,kk+1+ntime*nturns) = diff_psi;                                  %stores differences for diffusion calculations
    diff_zeta_n(:,kk+1+ntime*nturns) = diff_zeta;
 
endfor
endfor

if(is_flux_tube)
  save('Poincare1.mat','psi_n','zeta_n','coord_store','diff_psi_n','diff_zeta_n','psimaxdisp','ly');
else
  save('Poincare1.mat','psi_n','zeta_n','coord_store_cart','diff_psi_n','diff_zeta_n','psimaxdisp');
end

if(maxdisp_plot)
  figure(3);
  plot(psimaxdisp);
  xlabel('number of halfturns');
  ylabel('maximal displacement $\Delta \psi_N$');
  avg_pmd = mean(psimaxdisp);
  legend({['average =' num2str(avg_pmd)]});
  if(save_plots)
    print('epslatex', 'figure3.tex');                     %create eps file with plot and latex file with labels for the max displacement plot
  endif
endif


if(save_plots)
  figure(1);
  print('epslatex', 'figure1.tex');                       %create eps file with plot and latex file with labels for the Poincare plot
endif

if(diff_plot)
  FL_Diff
endif