% Calculate the miller's parameters from an hamada file and write them in a GKW input file
%	function [H,M] = hamada_miller(ham_file,ham_path,gkw_file,proj,rw,square_sw,zmil_sw,epsilon)
% Inputs :
%	ham_file : name of the hamada file
%	ham_path : path of the hamada file
%	gkw_file : name of the GKW input file
%	proj : name of the project where the GKW input file is
%	rw : 'r' if you only want to calculate the Miller's parameters, 'rw' if you also want to write them in a GKW file 
%   square_sw: (optional) integer switch, use 0 to set squareness and deriv to zero
%   zmil_sw: (optional) integer switch, use 0 to set Zmil and deriv to zero
%   epsilon (optional) select eps surface (chease value normalised by R_axis) to select
%
% Outputs :
%	H : contains the informations of hamada file
%	M : contains all the Miller's parameters
%
% WARNING : this version works fine with zeta ~ 0 (for numerical surfaces). High zeta produces surfaces that cannot be resolved using zmax etc... could be done with an fminsearch but not sure if it is worth.
% There is no known problems with surfaces perfectly described using miller's parameters.
% For up-down asymmetric equilibiria it is more common to calculate upper and lower values of
% each parameter, and then take the mean.
%
%
% YC 22.11.2012: 
% WARNING, this routine is not very robust and fails to properly calculate delta and zeta in many cases (particularly for real equilibrium)
% Spurious radial variations are obtained leading to random values of sdelta and szeta
% Use of rz2miller.m suggested


function [H,M] = hamada_miller(ham_file,ham_path,gkw_file,proj,rw,square_sw,zmil_sw,epsilon)

disp('Routine disabled (unreliable) - use of rz2miller.m suggested')
H=[];M=[];
return

if ~exist('rw')
  rw='rw';
end

if ~exist('square_sw')
  square_sw = 1;
end

if ~exist('zmil_sw')
  zmil_sw = 1;
end


H = read_hamada(ham_file,ham_path);

N1 = H.npsi;
N2 = H.ns;

H.r = H.r';
H.z = H.z';
H.theta = H.theta';
H.drdpsi = H.drdpsi';
H.dzdpsi = H.dzdpsi';

[rmax dum1] = max(H.r);
[rmin dum2] = min(H.r);

eps = (rmax-rmin)/2;

deps(2:N1-1) = abs(eps(3:N1)-eps(1:N1-2));

iinteg1 = sum(abs(H.z(3:N2,:)-H.z(1:N2-2,:)).* H.r(2:N2-1,:) .*H.z(2:N2-1,:));
iinteg2 = sum(abs(H.z(3:N2,:)-H.z(1:N2-2,:)).* H.r(2:N2-1,:));

Z0 = iinteg1 ./ iinteg2;
if (zmil_sw ==0) Z0(:) = 0; end

R0 = (rmax+rmin)/2;

zmin = min(H.z);

[zmax1 tmp2] = max(H.z);
[lim1 the1] = min(tmp2);
[lim2 the2] = max(tmp2);

% Used to have a smooth definition of delta (should be increased if delta presents discontinuities due to small discretisation)
Ninterp = 10;
dtheta = 2*pi / (N2 * Ninterp);
theta = (0:dtheta:2*pi);
N3 = length(theta);

for ith = 1:N1

  R1(:,ith) = interp1(H.theta(:,ith),H.r(:,ith),theta,'spline');
  Z1(:,ith) = interp1(H.theta(:,ith),H.z(:,ith),theta,'spline'); 
  [zmax tmp] = max(Z1(:,ith));
  [zmax3 tmp3] = max(Z1(1:tmp-20,ith));
  
  if (zmax - 1E-06 < zmax3 & zmax3 < zmax + 1E-06)
      [zmax4 tmp4] = max(Z1(1:tmp3-20,ith));
      if (zmax3 - 1E-06 < zmax4 & zmax4 < zmax3 + 1E-06)     
	delta(ith) = (R0(ith)-R1(tmp3,ith))./eps(ith); 
        kappa(ith) = (zmax-zmin(ith))./(2*eps(ith));	
      end
  else
    delta(ith) = (R0(ith)-R1(tmp,ith))./eps(ith); 
    kappa(ith) = (zmax-zmin(ith))./(2*eps(ith));  
  end  
  
end


for ith =1:N2

  R(ith,:) = R0+eps.*cos(H.theta(ith,:)+asin(delta).*sin(H.theta(ith,:)));

end

for ith = 1:N1
    
  R2(:,ith) = interp1(H.theta(:,ith),R(:,ith),theta,'spline');
  
  [a tmp] = min(abs(theta-pi/4));
  [a tmp1] = min(abs(R2(tmp,ith)-R1(1:round(N3/2),ith)));
  zeta(ith) = asin(Z1(tmp1,ith)/(kappa(ith)*eps(ith))) - pi/4;
  if (square_sw==0) zeta(:) = 0; end
     
end


for ith =1:N2 

  Z(ith,:) = Z0 + eps.*kappa.*sin(H.theta(ith,:)+zeta.*sin(2*H.theta(ith,:)));
  
end

%  for ith = 1:N1
%  
%      Z2(:,ith) = interp1(H.theta(:,ith),Z(:,ith),theta,'spline');
%  
%      thetatest(:,ith) = atan((Z2(:,ith)-Z0(1))./(R2(:,ith)-R0(1)));
%      diff(:,ith) = abs(thetatest(2:end,ith)-thetatest(1:end-1,ith));
%      [a b1] = max(diff(:,ith));
%      [a b2] = max(diff(1:b1-10,ith));
%      thetatest(b2+1:end,ith) = thetatest(b2+1:end,ith) + pi; 
%      thetatest(b1+1:end,ith) = thetatest(b1+1:end,ith) + pi;        
    
%      for ii = 1:N3
%        [a b(ii)] = min(abs(theta - thetatest(ii,ith)));
%      end
    
%      R1(:,ith) = R1(b(:),ith);
%      Z1(:,ith) = Z1(b(:),ith);          
    
%      Norme(ith) = sum(sqrt((R1(:,ith) - R2(:,ith)).^2 + (Z1(:,ith) - Z2(:,ith)).^2));
%  end

%  for ith = 1:N1
%  
%      x0 = [kappa(ith), delta(ith), zeta(ith)];
%      func = @(x)sum(sqrt((R1(:,ith) - R0(ith)-eps(ith)*cos(theta(:) + asin(x(2)) * sin(theta(:)))).^2 + (Z1(:,ith) - Z0(ith) - eps(ith)*x(1)*sin(theta(:) + x(3)*sin(2*theta(:)))).^2));
%      [x,fval,exitflag] = fminsearch(func,x0);
%      kappa(ith) = x(1);
%      delta(ith) = x(2);
%      zeta(ith) = x(3);
%      R2(:,ith) = R0(ith) + eps(ith) * cos(theta + asin(delta(ith))*sin(theta));
%      Z2(:,ith) = Z0(ith) + eps(ith) * kappa(ith) * sin(theta + zeta(ith) * sin(2* theta));
%      Norme2(ith) = sum(sqrt((R1(:,ith) - R2(:,ith)).^2 + (Z1(:,ith) - Z2(:,ith)).^2));
%  
%    
%  end

szeta(2:N1-1) = eps(2:N1-1).*(zeta(3:N1)-zeta(1:N1-2))./(eps(3:N1)-eps(1:N1-2));

skappa(2:N1-1) = eps(2:N1-1).*(kappa(3:N1)-kappa(1:N1-2))./(kappa(2:N1-1).*(eps(3:N1)-eps(1:N1-2)));

sdelta(2:N1-1) = eps(2:N1-1).*(delta(3:N1)-delta(1:N1-2))./((eps(3:N1)-eps(1:N1-2)).*sqrt(1-delta(2:N1-1).^2));

dRmil(2:N1-1) = (R0(3:N1)-R0(1:N1-2))./(eps(3:N1)-eps(1:N1-2));

dZmil(2:N1-1) = (Z0(3:N1)-Z0(1:N1-2))./(eps(3:N1)-eps(1:N1-2));

% Determination of alpha
theta = theta';
theta = repmat(theta(:),1,N1);
eps1 = repmat(eps(1,:),N3,1);

dRdth(2:N3-1,:) = (R1(3:N3,:) - R1(1:N3-2,:))./(theta(3:N3,:)-theta(1:N3-2,:));
dZdth(2:N3-1,:) = (Z1(3:N3,:) - Z1(1:N3-2,:))./(theta(3:N3,:)-theta(1:N3-2,:));

dRdr(:,2:N1-1) =  (R1(:,3:N1)-R1(:,1:N1-2))./abs(eps1(:,3:N1)-eps1(:,1:N1-2));
dZdr(:,2:N1-1) =  (Z1(:,3:N1)-Z1(:,1:N1-2))./abs(eps1(:,3:N1)-eps1(:,1:N1-2));

d2Rdrdth(2:N3-1,2:N1-1) = (dRdr(3:N3,2:N1-1)-dRdr(1:N3-2,2:N1-1))./abs(theta(3:N3,2:N1-1)-theta(1:N3-2,2:N1-1));
d2Zdrdth(2:N3-1,2:N1-1) = (dZdr(3:N3,2:N1-1)-dZdr(1:N3-2,2:N1-1))./abs(theta(3:N3,2:N1-1)-theta(1:N3-2,2:N1-1));

S(1,:)=abs(sum(dtheta*dRdth(2:N3-1,:).*Z1(2:N3-1,:)));

Rg(1,:)=sum(dtheta*R1(2:N3-1,:).*sqrt((dRdth(2:N3-1,:)).^2+(dZdth(2:N3-1,:)).^2))./sum(dtheta*sqrt((dRdth(2:N3-1,:)).^2+(dZdth(2:N3-1,:)).^2));

V(1,:)=2*pi*Rg(1,:).*S(1,:);

dSdr(1,2:N1-1)=sum(dtheta*(dZdr(2:N3-1,2:N1-1).*dRdth(2:N3-1,2:N1-1)+Z1(2:N3-1,2:N1-1).*d2Rdrdth(2:N3-1,2:N1-1)));

dRgdr(1,2:N1-1)=sum(dtheta*(dRdr(2:N3-1,2:N1-1).*sqrt((dRdth(2:N3-1,2:N1-1)).^2+(dZdth(2:N3-1,2:N1-1)).^2)+R1(2:N3-1,2:N1-1).*(d2Rdrdth(2:N3-1,2:N1-1).*dRdth(2:N3-1,2:N1-1)+d2Zdrdth(2:N3-1,2:N1-1).*dZdth(2:N3-1,2:N1-1))./sqrt((dRdth(2:N3-1,2:N1-1)).^2+(dZdth(2:N3-1,2:N1-1)).^2)))./sum(dtheta*sqrt((dRdth(2:N3-1,2:N1-1)).^2+(dZdth(2:N3-1,2:N1-1)).^2))-sum(dtheta*R1(2:N3-1,2:N1-1).*sqrt((dRdth(2:N3-1,2:N1-1)).^2+(dZdth(2:N3-1,2:N1-1)).^2)).*sum((d2Rdrdth(2:N3-1,2:N1-1).*dRdth(2:N3-1,2:N1-1)+d2Zdrdth(2:N3-1,2:N1-1).*dZdth(2:N3-1,2:N1-1))*dtheta./sqrt((dRdth(2:N3-1,2:N1-1)).^2+(dZdth(2:N3-1,2:N1-1)).^2))./sum(dtheta*sqrt((dRdth(2:N3-1,2:N1-1)).^2+(dZdth(2:N3-1,2:N1-1)).^2)).^2;

H.damindpsi = H.damindpsi';
vprime(2:N1-1)=2*pi*(dRgdr(2:N1-1).*S(2:N1-1)+Rg(2:N1-1).*dSdr(2:N1-1)).*H.damindpsi(2:N1-1);

H.dpdpsi = H.dpdpsi';

alpha(2:N1-1) = -2*vprime(2:N1-1)*4*pi*10^(-07).*H.dpdpsi(2:N1-1).*sqrt(V(2:N1-1)./(2*(pi^2)*R0(2:N1-1)))/(2*pi)^2;

% Normalised Miller's parameters for all flux surfaces
M.eps = eps./R0;
M.kappa = kappa;
M.delta = delta;
M.zeta = zeta;
M.skappa(2:N1-1) = skappa(2:N1-1);
M.sdelta(2:N1-1) = sdelta(2:N1-1);
M.szeta(2:N1-1) = szeta(2:N1-1);
M.R0 = R0;
M.Z0 = Z0./R0;
M.dRmil(2:N1-1) = dRmil(2:N1-1);
M.dZmil(2:N1-1) = dZmil(2:N1-1);
M.alpha = alpha;

if exist('epsilon') % return miller params structure for a single flux surface in format for gkw file
    
    clear M
    
    epsn = eps./R0;            % miller epsilon normalisation
    epsn_chease = eps/H.r0exp; % chease epsilon normalisation
    
    [a b] = min(abs(epsilon - epsn_chease));
    M.eps = epsn(b);
    M.kappa = kappa(b);
    M.delta = delta(b);
    M.square = zeta(b);
    M.skappa = skappa(b);
    M.sdelta = sdelta(b);
    M.ssquare = szeta(b);
    %M.R0 = R0(b);
    M.Zmil = Z0(b)/R0(b);
    M.dRmil = dRmil(b);
    M.dZmil = dZmil(b);
    M.gradp = alpha(b);
    M.gradp_type = 'alpha';
    M.shat = eps(b)*H.dqdpsi_chk(b)/(H.damindpsi(b)*H.q(b));
    M.q = H.q(b);
    
    plot(H.r(:,b),H.z(:,b),'b');
    hold on;
    plot(R(:,b),Z(:,b),'r');
    axis equal;
    title('Surfaces de flux');
    xlabel('R');
    ylabel('Z');
    
else
    
    plot(H.r,H.z,'b');
    hold on;
    plot(R,Z,'r');
    title('Surfaces de flux');
    axis equal;
    xlabel('R');
    ylabel('Z');
    
end

if (rw == 'rw')
  %Warning this epsilon is not the usual GKW normalised epsilon !
  lim_dwn_eps = num2str(eps(1),4);
  lim_up_eps = num2str(eps(end),4);
  epsilon = input(['Choose amin = R0*eps (not normalised) between ' lim_dwn_eps ' and ' lim_up_eps ' : ']);

  [a b] = min(abs(epsilon - eps));
  if (b == 1 | b == N1)
    epsilon = input(['Choose another amin : ']);
  end
  modify_gkwinput(gkw_file,proj,{'kappa','delta','square','skappa','sdelta','ssquare','dRmil','dZmil','gradp','gradptype','Q','SHAT'},{M.kappa(b),M.delta(b),M.zeta(b),M.skappa(b),M.sdelta(b),M.szeta(b),M.dRmil(b),M.dZmil(b),M.alpha(b),'alpha',H.q(b),eps(b)*H.dqdpsi(b)/(H.damindpsi(b)*H.q(b))},{0,0,0,0,0,0,0,0,0,0,0,0});

end
