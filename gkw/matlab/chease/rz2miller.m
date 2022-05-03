% Computes the Miller parameters from  R,Z flux surfaces description
% 
%     [r,Rmil,Zmil,k,d,z,dRmildr,dZmildr,sk,sd,sz,Rout,Zout,chi2]=rz2miller(R,Z,code,doplots)
%
% Inputs:
%  R,Z:   flux surfaces description in [m]
%         Size of the array assumed to be (nrho,npol), take care of it: no check performed!
%         The array is also assumed to have no double points
%  code:  Switch between different code conventions.
%         Available are: 'gkw', 'gyro', 'tglf' 
%  doplots: if 1 (default) perform plots to check the quality of the parametrisation
%
% Outputs:
%  r,Rmil,Zmil,k,d,z,dRmildr,dZmildr,sk,sd,sz: Miller parameters, as described in the GKW manual
%  Rout, Zout:  parametrised flux surfaces in [m]
%  chi2: parametrisation quality (chi2>0.05 starts to be quite bad, switch to full equilibrium advised!)
%
% Triangularity and squareness obtained by minimising the distance between the input flux surfaces and the 
% parametrised solutions.
%
% Requires interpos, available from https://crppwww.epfl.ch/~sauter/interpos/
%
%
% YC 22.11.2012

function [r,Rmil,Zmil,k,d,z,dRmildr,dZmildr,sk,sd,sz,Rout,Zout,chi2]=rz2miller(R,Z,code,doplots)

if (~exist('doplots','var')||isempty(doplots))
 doplots=1;
end 

if (~exist('code','var')||isempty(code))
 code='gkw';
end 

if ~(strcmp(code,'gkw')|strcmp(code,'gyro')|strcmp(code,'tglf'))
 disp(['Miller parametrisation not available for code=' code])
 disp('Available options are: gkw, gyro and tglf')
 return
end

if ~exist('doplots')||isempty(doplots)
 doplots = true;
end 

% number of points for refined grids
Nt=1000;

%% Rmil, Zmil
Nr=size(R,1);

% elevation, as defined in Candy PPCF2009, supposed to be more robust than (Zmax+Zmin)/2
Ztmp1 = Z; 
Ztmp2 = circshift(Z,[0 1]);
Ztmp1(:,end+1) = Ztmp1(:,1);
Ztmp2(:,end+1) = Ztmp2(:,1);
dZ = 0.5*(abs(diff(Ztmp1,[],2)) + abs(diff(Ztmp2,[],2)));
Zmil = sum(R.*Z.*dZ,2)./sum(R.*dZ,2); 

% reference major radius: first rough estimate
Rmil = (max(R,[],2) + min(R,[],2))./2;

% now get the accurate value (code dependent) and the minor radius
r = NaN.*Rmil;
for ii=1:Nr
 Ztmp1 = Z(ii,R(ii,:)>Rmil(ii)); [Ztmp1 I]=sort(Ztmp1); 
 Rtmp1 = R(ii,R(ii,:)>Rmil(ii)); Rtmp1=Rtmp1(I);
 Ztmp2 = Z(ii,R(ii,:)<Rmil(ii)); [Ztmp2 I]=sort(Ztmp2); 
 Rtmp2 = R(ii,R(ii,:)<Rmil(ii)); Rtmp2=Rtmp2(I);

 switch code
  case 'gkw'  % Rmax and Rmin are extrema over the whole flux surface
   Zg=linspace(min(Ztmp1),max(Ztmp1),Nt);
   Rg=interpos(Ztmp1,Rtmp1,Zg,-0.01);
   Rmax = max(Rg);
   Zg=linspace(min(Ztmp2),max(Ztmp2),Nt);
   Rg=interpos(Ztmp2,Rtmp2,Zg,-0.01);
   Rmin = min(Rg);
   
  case {'gyro','tglf'} % Rmax and Rmin calculated at centroid elevation
   Rg=interpos(Ztmp1,Rtmp1,[Zmil(ii) Zmil(ii)],-0.01);
   Rmax = Rg(1);
   Rg=interpos(Ztmp2,Rtmp2,[Zmil(ii) Zmil(ii)],-0.01);
   Rmin = Rg(1);
 end

 Rmil(ii) = (Rmax + Rmin)/2;
 r(ii) = (Rmax - Rmin)/2;
end

%% Computes a poloidal angle (th) zero at LFS midplane, positive going upward 
% Defined as
%   R-Rmil = a cos(th)
%   Z-Zmil = a sin(th)

Rr = R - repmat(Rmil,1,size(R,2));
Zz = Z - repmat(Zmil,1,size(Z,2));

ratio = Zz./Rr;
th = atan(ratio);

% pi shifts to have th in [0, 2*pi]
th(Rr<0)=th(Rr<0)+pi;
th(Rr>0&Zz<0)=th(Rr>0&Zz<0)+2*pi;

% Sort from 0 to 2*pi
[th,Isort]=sort(th,2,'ascend');
for ii = 1:Nr, 
 R(ii,:) = R(ii,Isort(ii,:)); 
 Z(ii,:) = Z(ii,Isort(ii,:)); 
 Rr(ii,:) = Rr(ii,Isort(ii,:)); 
 Zz(ii,:) = Zz(ii,Isort(ii,:));  
end


% detect phase jumps
[I,J]=find(abs(diff(th,[],2))>3);
if ~isempty(I) % phase jump detected
  disp('Large jumps in th, something went wrong')
  disp('Check your R,Z inputs, program stopped')
  return
end

% local minor radius
a = sqrt(Rr.^2 + Zz.^2);


%% Now get the shape parameters
% average minor radius
r = (max(R,[],2) - min(R,[],2))./2;

% elongation - could be improved as for the calculation of Rmil (i.e. interpolation on a finer grid to get Zmax and Zmin)
k = (max(Z,[],2) - min(Z,[],2))./(2*r);

% get delta and zeta from a minimisation
ths=linspace(0,2*pi,Nt);
[d,z,chi2]=deal(0.*k);
[Rout,Zout]=deal(zeros(Nr,Nt));
x0=[0 0];

opt=optimset('TolFun',1e-5,'TolX',1e-5);
for ii=1:Nr
 as = interpos(th(ii,:),a(ii,:),ths,0.,-1,2*pi);  

 [xout,chi2(ii)] = fminsearch(@(x) da(x,r(ii),k(ii),ths,as),x0,opt); 
 
 d(ii) = sin(xout(1));
 z(ii) = xout(2);

 Rout(ii,:) = Rmil(ii) + r(ii).*cos(ths + xout(1).*sin(ths)); 
 Zout(ii,:) = Zmil(ii) + r(ii).*k(ii).*sin(ths+xout(2).*sin(2.*ths));

 x0 = xout;

end
chi2=sqrt(chi2).*100;
if any(chi2>0.05)
 disp('Warning, for chi2>0.05 parametrisation quality is poor (see plots), better to use the full geometry')
 doplots=1;
end

%% Computes the radial derivatives 

[tmp1 tmp2]=interpos(13,r,Rmil,r,-5);
dRmildr=tmp2;

[tmp1 tmp2]=interpos(13,r,Zmil,r,-5);
dZmildr=tmp2;

[tmp1 tmp2]=interpos(13,r,k,r,-5);
sk=r./k.*tmp2;

[tmp1 tmp2]=interpos(13,r,d,r,-10);
sd=r./sqrt(1-d.^2).*tmp2;
switch code
 case 'gkw'  
  sd=r./sqrt(1-d.^2).*tmp2;

 case {'gyro','tglf'} 
  sd=r.*tmp2;
end

[tmp1 tmp2]=interpos(13,r,z,r,-10);
sz=r.*tmp2;

%% plots
if doplots

 figure
 plot(r./max(r),chi2)
 hold on
 plot(r./max(r),repmat(0.05,1,Nr),'r')
 xlabel('r/r_{max}')
 ylabel('chi2')
 axis([0 1 0 0.5]) 

 figure
 plot(r./max(r),d,'b')
 hold on
 plot(r./max(r),z,'b--')
 xlabel('r/r_{max}') 
 ylabel('\delta  and  \zeta')

 figure
 h1=plot(R',Z','b');
 hold on
 h2=plot(Rout',Zout','r');
 xlabel('R')
 ylabel('Z')
 legend([h1(1),h2(1)],'full','Miller')
 axis equal

end 


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function chi2=da(x0,r0,k0,ths,as)

if abs(sin(x0(1)))>1 | abs(x0(2))>0.7
 chi2=NaN;
else
 Nt0=1000;
 t0 = linspace(0,2*pi,Nt0);
 t0(end)=[];
 Rr0 = r0.*cos(t0 + x0(1).*sin(t0));
 Zz0 = r0.*k0.*sin(t0 + x0(2).*sin(2*t0));

 th0 = atan(Zz0./Rr0);

 % pi shifts to have th in [0, 2*pi]
 th0(Rr0<0)=th0(Rr0<0)+pi;
 th0(Rr0>0&Zz0<0)=th0(Rr0>0&Zz0<0)+2*pi;

 a0 = sqrt(Rr0.^2 + Zz0.^2);

 a0s = interpos(th0,a0,ths,0.,-1,2*pi);  

 chi2 = sum(((a0s-as)./as).^2)./Nt0.^2; % power of 2 to punish large excursions

end


