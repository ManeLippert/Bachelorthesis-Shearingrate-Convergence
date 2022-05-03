function [gkwmil eql yMill]= eqdsk_miller(eqdsk_file,eqdsk_path,epsilon,optplot,eqm);
%
% function [gkwinput gs2input gs2_mill_all_surfs] = (eqdsk_file,eqdsk_path,epsilon,optplot,eqm);
%
% builds the Miller equilibrium from eqdsk file (or hamada, in principle)
%
% Requires interpos, available from https://crppwww.epfl.ch/~sauter/interpos/x
%
% WARNING: this produces quite different results to hamada_miller
% due to either better smoothing / interpolation or averaging over top and bottom halves
%
% Derives GS2 miller parameters then converts to GKW at the end
% in output, GKW miller parameters
%
% optplot = 0 : no plots
% optplot = 1 : plot eqdsk surfaces
% optplot = 2 : plot eqdsk + contourc
% optplot = 3 : plot miller surfaces + eqdsk + contourc (default)
%
% eqm = 'eqdsk' read eqdsk file (default)
% eqm = 'hamada' read gkw-chease hamada file (not yet working)
%
% TO DO: also calculate q and shat from eqdsk
%
% cla, 11.01.2007
% fjc, 31.07.2012
%

if ~exist('optplot'); optplot = 3; end
if ~exist('eqm'); eqm='eqdsk'; end

%number or surfaces to remove (adjust if some give problems)
bnd=3;

%% Get the flux surfaces into same format
%% and calculate the q profile
% yout = RZsurfs(eqdg, shot, t1,t2,1, expeq, editeq); close(gcf);
%Zmag = yout.Zin(1);

if (isequal(eqm,'hamada'))
    error('hamada option still has bugs')
    ham = read_hamada(eqdsk_file,eqdsk_path);
    Zmag=mean(ham.z(:,1));
    for ij=1:size(ham.r,2)
      %first flux surface must be outermost ?
      ik = size(ham.r,2) - ij + 1;
      yout.flx_srfs(ij).r=ham.r(:,ik)';
      yout.flx_srfs(ij).z=ham.z(:,ik)';
    end 
    bnd=8;
elseif (isequal(eqm,'eqdsk'))
    yout = read_eqdsk_cla2(eqdsk_file,eqdsk_path,optplot);    
    Zmag = yout.zmaxis;
else
    error('only hamada or eqdsk supported')
end

Ar_mid = NaN.*ones(1,length(yout.flx_srfs));
ARmaj = Ar_mid; Aelong = Ar_mid; Atri = Ar_mid;
Br_mid = NaN.*ones(1,length(yout.flx_srfs));
BRmaj = Br_mid; Belong = Br_mid; Btri = Br_mid;
ABr_mid = NaN.*ones(1,length(yout.flx_srfs));
ABRmaj = ABr_mid; ABelong = ABr_mid; ABtri = ABr_mid;
%
%
for ij =  bnd:length(yout.flx_srfs)-bnd;
    zR = yout.flx_srfs(ij).r;
    zZ = yout.flx_srfs(ij).z;
    %
    %% Build pseudo up-down symmetry in the two cases (top and bottom).
    zZr = 2*Zmag-zZ;
    iop = find(zZr > zZ);
    if (size(iop) < 2) continue; end
    clear RzR RzZ RzZr;
    RzZr(1:length(zZr)-iop(end)) = zZr(iop(end)+1:end);
    RzR(1:length(zZr)-iop(end)) = zR(iop(end)+1:end);
    RzZ(1:length(zZr)-iop(end)) = zZ(iop(end)+1:end);
    RzZ(length(zZr)-iop(end)+1:length(zZr)) = zZ(1:iop(end));
    RzZr(length(zZr)-iop(end)+1:length(zZr)) = zZr(1:iop(end));
    RzR(length(zZr)-iop(end)+1:length(zZr)) = zR(1:iop(end));
    RzZs = RzZ; RzZrs = RzZr;
    %xxx = (1:length(RzZ))./length(RzZ);
    %RzZrs = interpos(13,xxx, RzZr, xxx, 1e-06);auto
    %RzZs = interpos(13,xxx, RzZ, xxx, 1e-06);
    
    iop = find(RzZrs > RzZs);
    io1 = iop(1); if io1 == 1; io1 = iop(2) ; end;
    io1 = io1-1; io2 = iop(end);
    if io1 > 1;
        xx1 = (1:io1)./io1;
        xx2 = (1:io2-io1)./(io2-io1);
        aaa = interp1(xx1, RzZrs(io1:-1:1), xx2);
        bbb = interp1(xx2, RzZrs(io2:-1:io1+1), xx1);
        ARzZ = RzZs; BRzZ = RzZs;
        ARzZ(io1+1:io2) = aaa;
        BRzZ(1:io1) = bbb;
        xxx = 1:length(ARzZ);
        ARzZ(isnan(ARzZ)) = interp1(xxx(~isnan(ARzZ)), ARzZ(~isnan(ARzZ)), xxx(isnan(ARzZ)));
        BRzZ(isnan(BRzZ)) = interp1(xxx(~isnan(BRzZ)), BRzZ(~isnan(BRzZ)), xxx(isnan(BRzZ)));
        
        
        %% Compute Miller Geometry in the two cases
        zRN = RzR; zZN = ARzZ;
        theta= linspace(0,2*pi, length(zRN));
        [a jzmin] = min(zZN); [a jzmax] = max(zZN);
        ARmaj(ij) = (max(zRN)+min(zRN))./2;
        Ar_mid(ij) = (max(zRN)-min(zRN))./2;
        Aelong(ij) = (max(zZN)-min(zZN))./(max(zRN)-min(zRN));
        Atriup(ij) = -asin((zRN(jzmax)-ARmaj(ij))/Ar_mid(ij)); %% Miller definition
        Atridn(ij) = -asin((zRN(jzmin)-ARmaj(ij))/Ar_mid(ij)); %% Miller definition
        Atri(ij) = (Atriup(ij) + Atridn(ij))./2;
        AzRMil = (ARmaj(ij)+Ar_mid(ij)*cos(theta+Atri(ij)*sin(theta)));
        AzZMil = (max(zZN)+min(zZN))./2+Aelong(ij)*Ar_mid(ij)*sin(theta);
        
        zRN = RzR; zZN = BRzZ;
        [a jzmin] = min(zZN); [a jzmax] = max(zZN);
        BRmaj(ij) = (max(zRN)+min(zRN))./2;
        Br_mid(ij) = (max(zRN)-min(zRN))./2;
        Belong(ij) = (max(zZN)-min(zZN))./(max(zRN)-min(zRN));
        Btriup(ij) = -asin((zRN(jzmax)-BRmaj(ij))/Br_mid(ij)); %% Miller definition
        Btridn(ij) = -asin((zRN(jzmin)-BRmaj(ij))/Br_mid(ij)); %% Miller definition
        Btri(ij) = (Btriup(ij) + Btridn(ij))./2;
        theta= linspace(0,2*pi, length(zRN));
        BzRMil = (BRmaj(ij)+Br_mid(ij)*cos(theta+Btri(ij)*sin(theta)));
        BzZMil = (max(zZN)+min(zZN))./2+Belong(ij)*Br_mid(ij)*sin(theta);
        
        ABRmaj(ij) = (ARmaj(ij)+BRmaj(ij))./2;
        ABr_mid(ij) = (Ar_mid(ij)+Br_mid(ij))./2;
        ABtri(ij)  = (Atri(ij)+Btri(ij))./2;
        ABelong(ij) = (Aelong(ij)+Belong(ij))./2;
        ABzRMil = (ABRmaj(ij)+ABr_mid(ij)*cos(theta+ABtri(ij)*sin(theta)));
        ABzZMil = (max(zZ)+min(zZ))./2+ABelong(ij)*ABr_mid(ij)*sin(theta);
        
        if (optplot > 2)
            if mod(ij,8) == 0;
                figure(1001);
                plot(zR,zZ,'r-');
                axis equal; hold on;
                plot(RzR, ARzZ, 'b'); plot(RzR, BRzZ, 'k');
                op=plot(AzRMil, AzZMil, 'g');set(op,'color', [0 0.6 0]);
                op=plot(BzRMil, BzZMil, 'm');set(op,'color', [0.3 0.0 0.3]);
                op=plot(ABzRMil, ABzZMil, 'm');set(op,'color', [1 0 0]);
            end
        end
        
    end; % if io1 > 1;
    
end

clear y;
y.Rmaj     = ABRmaj(~isnan(ABRmaj));
y.r_mid    = ABr_mid(~isnan(ABRmaj));
y.tri      = ABtri(~isnan(ABRmaj));
y.elong    = ABelong(~isnan(ABRmaj));

y.UP_Rmaj  = ARmaj(~isnan(ABRmaj));
y.UP_r_mid = Ar_mid(~isnan(ABRmaj));
y.UP_tri   = Atri(~isnan(ABRmaj));
y.UP_elong = Aelong(~isnan(ABRmaj));

y.DN_Rmaj  = BRmaj(~isnan(ABRmaj));
y.DN_r_mid = Br_mid(~isnan(ABRmaj));
y.DN_tri   = Btri(~isnan(ABRmaj));
y.DN_elong = Belong(~isnan(ABRmaj));

tris = NaN*y.r_mid; elongs = NaN*y.r_mid; Rmajs = NaN*y.r_mid;
trisup = NaN*y.r_mid; trisdn = NaN*y.r_mid;
triuppri = NaN*y.r_mid; tridnpri = NaN*y.r_mid;
shift = NaN*y.r_mid; tripri = NaN*y.r_mid; akappri = NaN*y.r_mid;
zrmd = y.r_mid;
[aa bb] = interpos(13, zrmd(end-bnd:-1:bnd), y.Rmaj(end-bnd:-1:bnd), ...
    zrmd(end-bnd:-1:bnd), 1e-05, [0 0], [0 0]);
Rmajs(bnd:length(zrmd)-bnd) = aa(end:-1:1);
shift(bnd:length(zrmd)-bnd) = bb(end:-1:1);
[aa bb] = interpos(13, zrmd(end-bnd:-1:bnd), y.elong(end-bnd:-1:bnd), ...
    zrmd(end-bnd:-1:bnd), 5e-05, [0 0], [0 0]);
elongs(bnd:length(zrmd)-bnd) = aa(end:-1:1);
akappri(bnd:length(zrmd)-bnd) = bb(end:-1:1);
[aa bb] = interpos(13, zrmd(end-bnd:-1:bnd), y.tri(end-bnd:-1:bnd), ...
    zrmd(end-bnd:-1:bnd), 5e-04, [0 0], [0 0]);
tris(bnd:length(zrmd)-bnd) = aa(end:-1:1);
tripri(bnd:length(zrmd)-bnd) = bb(end:-1:1);
[aa bb] = interpos(13, zrmd(end-bnd:-1:bnd), y.UP_tri(end-bnd:-1:bnd), ...
    zrmd(end-bnd:-1:bnd), 5e-04, [0 0], [0 0]);
trisup(bnd:length(zrmd)-bnd) = aa(end:-1:1);
triuppri(bnd:length(zrmd)-bnd) = bb(end:-1:1);
[aa bb] = interpos(13, zrmd(end-bnd:-1:bnd), y.DN_tri(end-bnd:-1:bnd), ...
    zrmd(end-bnd:-1:bnd), 5e-04, [0 0], [0 0]);
trisdn(bnd:length(zrmd)-bnd) = aa(end:-1:1);
tridnpri(bnd:length(zrmd)-bnd) = bb(end:-1:1);
%
clear yGS2sm;
yGS2sm.r_mid = fliplr(zrmd); yGS2sm.Rmaj = fliplr(Rmajs); yGS2sm.shift = fliplr(shift);
yGS2sm.elong = fliplr(elongs); yGS2sm.akappri = fliplr(akappri);
yGS2sm.tri = fliplr(tris); yGS2sm.tripri = fliplr(tripri);
yGS2sm.triup = fliplr(trisup); yGS2sm.triuppri = fliplr(triuppri);
yGS2sm.tridn = fliplr(trisdn); yGS2sm.tridnpri = fliplr(tridnpri);
%
%check for NaN in the center
if (isnan(yGS2sm.Rmaj(1)) & ~isnan(yGS2sm.r_mid(1)) );
    yyy = interpos(13, yGS2sm.r_mid(2:end), yGS2sm.Rmaj(2:end), [ yGS2sm.r_mid ]);
    yGS2sm.Rmaj(1) = yyy(1);
    yyy = interpos(13, yGS2sm.r_mid(2:end), yGS2sm.shift(2:end), [ yGS2sm.r_mid ]);
    yGS2sm.shift(1) = yyy(1);
    yyy = interpos(13, yGS2sm.r_mid(2:end), yGS2sm.elong(2:end), [ yGS2sm.r_mid ]);
    yGS2sm.elong(1) = yyy(1);
    yyy = interpos(13, yGS2sm.r_mid(2:end), yGS2sm.akappri(2:end), [ yGS2sm.r_mid ]);
    yGS2sm.akappri(1) = yyy(1);
    yyy = interpos(13, yGS2sm.r_mid(2:end), yGS2sm.tri(2:end), [ yGS2sm.r_mid ]);
    yGS2sm.tri(1) = yyy(1);
    yyy = interpos(13, yGS2sm.r_mid(2:end), yGS2sm.tripri(2:end), [ yGS2sm.r_mid ]);
    yGS2sm.tripri(1) = yyy(1);
    yyy = interpos(13, yGS2sm.r_mid(2:end), yGS2sm.triup(2:end), [ yGS2sm.r_mid ]);
    yGS2sm.triup(1) = yyy(1);
    yyy = interpos(13, yGS2sm.r_mid(2:end), yGS2sm.triuppri(2:end), [ yGS2sm.r_mid ]);
    yGS2sm.triuppri(1) = yyy(1);
    yyy = interpos(13, yGS2sm.r_mid(2:end), yGS2sm.tridn(2:end), [ yGS2sm.r_mid ]);
    yGS2sm.tridn(1) = yyy(1);
    yyy = interpos(13, yGS2sm.r_mid(2:end), yGS2sm.tridnpri(2:end), [ yGS2sm.r_mid ]);
    yGS2sm.tridnpri(1) = yyy(1);
end

%
%extrapolates to full minor radius when bnd = 1;
if bnd==1
    zrrr = r2rgs2('rhot',[0 1], eqdg, shot, t1, t2, expeq, editeq)';
    yyy = interpos(13, yGS2sm.r_mid, yGS2sm.Rmaj, [zrrr(1) yGS2sm.r_mid zrrr(2)]);
    yGS2sm.Rmaj = yyy;
    yyy = interpos(13, yGS2sm.r_mid, yGS2sm.shift, [zrrr(1) yGS2sm.r_mid zrrr(2)]);
    yGS2sm.shift = yyy;
    yyy = interpos(13, yGS2sm.r_mid, yGS2sm.elong, [zrrr(1) yGS2sm.r_mid zrrr(2)]);
    yGS2sm.elong = yyy;
    yyy = interpos(13, yGS2sm.r_mid, yGS2sm.akappri, [zrrr(1) yGS2sm.r_mid zrrr(2)]);
    yGS2sm.akappri = yyy;
    yyy = interpos(13, yGS2sm.r_mid, yGS2sm.tri, [zrrr(1) yGS2sm.r_mid zrrr(2)]);
    yGS2sm.tri = yyy;
    yyy = interpos(13, yGS2sm.r_mid, yGS2sm.tripri, [zrrr(1) yGS2sm.r_mid zrrr(2)]);
    yGS2sm.tripri = yyy;
    yyy = interpos(13, yGS2sm.r_mid, yGS2sm.triup, [zrrr(1) yGS2sm.r_mid zrrr(2)]);
    yGS2sm.triup = yyy;
    yyy = interpos(13, yGS2sm.r_mid, yGS2sm.triuppri, [zrrr(1) yGS2sm.r_mid zrrr(2)]);
    yGS2sm.triuppri = yyy;
    yyy = interpos(13, yGS2sm.r_mid, yGS2sm.tridn, [zrrr(1) yGS2sm.r_mid zrrr(2)]);
    yGS2sm.tridn = yyy;
    yyy = interpos(13, yGS2sm.r_mid, yGS2sm.tridnpri, [zrrr(1) yGS2sm.r_mid zrrr(2)]);
    yGS2sm.tridnpri = yyy;
    yGS2sm.r_mid = [zrrr(1) yGS2sm.r_mid zrrr(2)]; close(gcf);
end;

 yMill=yGS2sm;

 %find points without NaNs 
 igM = find(isnan(yMill.elong)==0);
 
 % take an average about the flux surface of interest
 deps = 0.015;
 amin = max(yMill.r_mid);
 r_ov_a=epsilon*yMill.Rmaj;
 irokM = find(yMill.r_mid./yMill.Rmaj >= epsilon-deps & yMill.r_mid./yMill.Rmaj <= epsilon+deps);
 eql.shift = interp1(yMill.r_mid(igM)./yMill.Rmaj(igM), yMill.shift(igM), epsilon);
 eql.akappa = interp1(yMill.r_mid(igM)./yMill.Rmaj(igM), yMill.elong(igM), epsilon);
 eql.akappri = mean(yMill.akappri(irokM));
 eql.tri = interp1(yMill.r_mid(igM)./yMill.Rmaj(igM), yMill.tri(igM), epsilon);
 eql.tripri = max(0.01, mean(yMill.tripri(irokM)));
 %eql.beta_prime_input=0; % Needs to be calculated elsewhere

 gkwmil.GEOM.eps=epsilon;
 % Need to calculate
 % see RZsurfs to get rgs2 for PF grid and take derivative to get shat
 % remove bad
 % gkwmil.GEOM.q= 
 % gkwmil.GEOM.shat=
 
 %convert to GKW parameters
 [dum gkwmil]=gs2_gkw_gyro_miller(eql,gkwmil,'gs2');

end

