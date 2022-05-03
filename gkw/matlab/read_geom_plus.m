function G = read_geom_plus(flnm,proj)
% geom = read_geom_plus(filename,proj)
% Reads gkw geom file into a matlab structure
%
% and calculates some additional derived quantities 
% such as the flux surface area (normalised by Rref**2)
%
% Warning:
% dVdeps is Jacobian*Rref^3, and this is the value that should 
% be used for making dimensional power fluxes, 
% not the flux surface area.
%
% YC, first version, March 2013
% FC, commit,        April 2013

% read the GKW geom file
G=read_geom(flnm,proj);

% read the GKW version number
G.version = read_gkw_version(flnm,proj);

% work on one poloidal turn only
Istart=find(G.s_grid>=-0.5 & G.s_grid<=0.5,1,'first');
Iend=find(G.s_grid>=-0.5 & G.s_grid<=0.5,1,'last');
Ns=Iend-Istart+1;

% contravariant metric tensor
Gup=zeros(3,3,Ns);
jac_test=zeros(1,Ns);
for ii=Istart:Iend
    Gup(:,:,ii-Istart+1) = [G.g_eps_eps(ii)   G.g_eps_zeta(ii)   G.g_eps_s(ii);
        G.g_eps_zeta(ii)  G.g_zeta_zeta(ii)  G.g_zeta_s(ii);
        G.g_eps_s(ii)     G.g_zeta_s(ii)     G.g_s_s(ii)];
    jac_test(ii-Istart+1)=1./det(Gup(:,:,ii-Istart+1)); % -> just to check that the Jacobian is a flux function as expected
end
G.Gup=Gup./G.rref.^2;


% covariant metric tensor
Gdown=zeros(3,3,Ns);
for ii=1:Ns
    Gdown(:,:,ii)=inv(Gup(:,:,ii));
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Method 1 - area from arc length and centroid calculation
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

R=G.r(Istart:Iend)';  %.*G.rref;
Z=G.z(Istart:Iend)';  %.*G.rref;

if ( G.version < 3201 )
  disp('Using Z fix for chease geometry, cf issue 190')
  Z=G.z(Istart:Iend)/G.rref;   
  G.z=G.z./G.rref;  
end

Rmid=0.5*(max(R)+min(R));
Zmid=0.5*(max(Z)+min(Z));
R=R-Rmid;
Z=Z-Zmid;

% centroid position

Ztmp1 = Z;
Ztmp2 = circshift(Z,[0 1]);
Ztmp1(:,end+1) = Ztmp1(:,1);
Ztmp2(:,end+1) = Ztmp2(:,1);
dZ = 0.5*(abs(diff(Ztmp1,[],2)) + abs(diff(Ztmp2,[],2)));
Z0 = sum(abs(R).*Z.*dZ,2)./sum(abs(R).*dZ,2);
%Z0a = trapz(abs(R).*Z,dZ,2)./trapz(abs(R),dZ,2);
Z0 = Z0+Zmid;


Rtmp1 = R;
Rtmp2 = circshift(R,[0 1]);
Rtmp1(:,end+1) = Rtmp1(:,1);
Rtmp2(:,end+1) = Rtmp2(:,1);
dR = 0.5*(abs(diff(Rtmp1,[],2)) + abs(diff(Rtmp2,[],2)));
R0 = sum(R.*abs(Z).*dR,2)./sum(abs(Z).*dR,2);
%R0a = trapz(R.*abs(Z),dR,2)./trapz(abs(Z),dR,2);
R0=R0+Rmid;

% Note: A=sum(abs(R).*dZ,2)=sum(abs(Z).*dR,2) = cross section area

% cross-section arc length
dl = sqrt(dR.^2+dZ.^2);
L1 = sum(dl,2);

% flux surface area
S1 = L1.*2*pi.*R0;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Method 2 - same but taking the integral over phi for each length element instead of calculing the centroid position
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
dS=dl.*2.*pi.*(R+Rmid);
S2=sum(dS);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Method 3 - total area from the integral of the area element dA
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
ds=G.s_grid(2)-G.s_grid(1); % assumes uniform s grid
dzeta=1; % no quantity depends on zeta, the integral over zeta is therefore a factor unity

dA = sqrt(squeeze(Gdown(2,2,:).*Gdown(3,3,:)-Gdown(2,3,:).^2)).*ds.*dzeta; % This is in fact the jacobian of the metric tensor excluding the radial coordinate, i.e. keeping only a 2x2 matrix

% Seems to be wrong
S3=sum(dA);

% Store the derived quantities wanted for export
G.S1 = S1;
G.S2 = S2;
G.surface_area=S2;
G.L1 = L1;
%G.S3 = S3;
%G.Gup=Gup;
%G.Gdown=Gdown;

end



