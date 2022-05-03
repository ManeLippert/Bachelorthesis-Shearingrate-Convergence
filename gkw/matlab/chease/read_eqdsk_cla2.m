function yyy = read_eqdsk_cla2(filename, pth, optplot, numcont);
%
% function yyy = read_eqdsk_cla2(filename, pth, optplot, numcont);
%
% by cla; fjc: added flux surface computation by contourc
%
% optplot = 0 : no plots
% optplot = 1 : plot eqdsk surfaces (default)
% optplot = 2 : plot eqdsk + contourc
%

if ~exist('optplot'); optplot = 1; end

if ~exist('numcont')
numcont = 200;
end
if isempty(numcont);
numcont = 200;
end


pthflnm = [pth '/' filename]
%
fptr = fopen(pthflnm);

frewind(fptr);
aaa = fgetl(fptr);
stst = fseek(fptr, -8,0);

%aaa = fscanf(fptr, '%f',1);
nw = fscanf(fptr, '%f',1)
nh = fscanf(fptr, '%f',1)
%
rwid = fscanf(fptr, '%f',1);
zhei = fscanf(fptr, '%f',1);
rcentr = fscanf(fptr, '%f',1);
rleft  = fscanf(fptr, '%f',1);
zmid = fscanf(fptr, '%f',1);
%
rmaxis = fscanf(fptr, '%f',1);
zmaxis = fscanf(fptr, '%f',1);
simag  = fscanf(fptr, '%f',1);
sibry  = fscanf(fptr, '%f',1);
bcentr = fscanf(fptr, '%f',1);
%
current = fscanf(fptr, '%f',1);
simag  = fscanf(fptr, '%f',1);
xdum   = fscanf(fptr, '%f',1);
rmaxis = fscanf(fptr, '%f',1);
xdum   = fscanf(fptr, '%f',1);
%
zmaxis = fscanf(fptr, '%f',1);
xdum   = fscanf(fptr, '%f',1);
sibry  = fscanf(fptr, '%f',1);
xdum   = fscanf(fptr, '%f',1);
xdum   = fscanf(fptr, '%f',1);
%
fpol   = fscanf(fptr, '%f',nw);
pres   = fscanf(fptr, '%f',nw);
ffprim = fscanf(fptr, '%f',nw);
pprim  = fscanf(fptr, '%f',nw);
%
clear psirz
for ij=1:nh
    psirz(1:nw, ij) = fscanf(fptr, '%f', nw);
end
%
qpsi  = fscanf(fptr, '%f',nw);
%
nbbbs  = fscanf(fptr, '%f',1);
limitr = fscanf(fptr, '%f',1);
%
aaa = fscanf(fptr, '%f', [2 nbbbs]);
rbbbs = aaa(1,:);
zbbbs = aaa(2,:);
%
aaa = fscanf(fptr, '%f', [2 limitr]);

rlim = aaa(1,:);
zlim = aaa(2,:);
fclose(fptr);
%

sefit_R=linspace(rleft,rleft+rwid, nw);
sefit_Z=linspace(zmid-zhei./2, zmid+zhei./2,nh);
spsi_bar=transpose(linspace(0,1, nw));


yyy.nw = nw;
yyy.nh = nh;
%
yyy.rwid = rwid;
yyy.zhei = zhei;
yyy.rcentr = rcentr;
yyy.rleft  = rleft;
yyy.zmid = zmid;
yyy.simag  = simag;
yyy.sibry  = sibry;
yyy.bcentr = bcentr;
yyy.current = current;
%
yyy.rmaxis=rmaxis;
yyy.zmaxis=zmaxis;
yyy.sefit_R=sefit_R;
yyy.sefit_Z=sefit_Z;
yyy.psirz = psirz;
yyy.rbbbs=rbbbs;
yyy.zbbbs=zbbbs;
yyy.rlim=rlim;
yyy.zlim=zlim;
yyy.spsi_bar = spsi_bar;
yyy.qpsi = qpsi;
yyy.fpol = fpol;
yyy.pres = pres;
yyy.pprim = pprim;
yyy.ffprim = ffprim;

[a idx_max] = max(yyy.zbbbs);
[a idx_min] = min(yyy.zbbbs);
idx_1 = min(idx_max, idx_min);
idx_2 = max(idx_max, idx_min);
%
if max(yyy.rbbbs(idx_1:idx_2)) > max(yyy.rbbbs(1:idx_1));
    %then from idx_1 to idx_2 is outerleg
    idx_outerleg = idx_1:idx_2;
    idx_innerleg = [1:idx_1-1 idx_2+1:length(yyy.rbbbs)-1];
else
    idx_outerleg = [1:idx_1-1 idx_2+1:length(yyy.rbbbs)-1];
    idx_innerleg = idx_1:idx_2;
end

[zouterleg idso] = sort(yyy.zbbbs(idx_outerleg));
[zinnerleg idsi] = sort(yyy.zbbbs(idx_innerleg));
rinnerlegz = yyy.rbbbs(idx_innerleg);
routerlegz = yyy.rbbbs(idx_outerleg);
%
routerleg = routerlegz(idso);
rinnerleg = rinnerlegz(idsi);

%zinnerleg(find(diff(zinnerleg)==0)+1)=zinnerleg(find(diff(zinnerleg)==0)+1)+1e-05;
%zouterleg(find(diff(zinnerleg)==0)+1)=zouterleg(find(diff(zinnerleg)==0)+1)+1e-05;
% rout_maxis=interp1(zouterleg, routerleg, yyy.zmaxis);
% [zinnerleg inds] = unique(zinnerleg);
% rinnerleg = rinnerleg(inds);
% rin_maxis=interp1(zinnerleg, rinnerleg, yyy.zmaxis);
% 
% yyy.aminor = (rout_maxis-rin_maxis)/2;
% yyy.Rgeo = (rout_maxis+rin_maxis)/2;

if (optplot > 0)
    figure;
    set(gcf, 'position', [750   250   400  400*zhei./rwid*0.95])
    contour(sefit_R, sefit_Z, transpose(psirz), 50);
    hold on;
    plot(rbbbs, zbbbs, 'r');
    %plotos(rlim, zlim, 'k', [2 0]);
    plot(rlim,zlim,'k-')
    plot(rmaxis, zmaxis, 'r+');
    axis([min(sefit_R) max(sefit_R) min(sefit_Z) max(sefit_Z)]);
    axis equal;
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

RR = sefit_R;
ZZ = sefit_Z;
PFL = transpose(psirz);
clear ypfm;

minpfl = min(min(PFL));
maxpfl = max(max(PFL));

%pfx = log(linspace(exp(minpfl), exp(maxpfl), numcont));
pfx = exp(linspace(log(minpfl), log(maxpfl), numcont));

CC = contourc(RR,ZZ,PFL,[pfx]);
jj = 1;
ij = 1;
clear xx;

%keyboard

while jj < length(CC(1,:));
    val = CC(1,jj);
    nbval = CC(2,jj);
    xx(ij).r = CC(1,jj+1:jj+nbval);
    xx(ij).z = CC(2,jj+1:jj+nbval);
    xx(ij).v = CC(1,jj);
    xx(ij).n = CC(2,jj);
    nb(ij) = CC(2,jj);
    zvalue(ij) = CC(1,jj);
    jj = jj+nbval+1;
    ij = ij+1;
end;

%keyboard

%remove flux surfaces outside separatix
for ij = 1:length(xx)
  dist_sq(ij) = max((xx(ij).r-rmaxis).^2 + (xx(ij).z-zmaxis).^2);
end
dist_bd = max((rbbbs-rmaxis).^2 + (zbbbs-zmaxis).^2);
% must be sorted in this order to work with build_miller_equilibrium
[dist_sq ind] = sort(dist_sq,'descend');
xx=xx(ind);
[dum ind] = find(dist_sq < dist_bd);
xx=xx(ind);
yyy.flx_srfs = xx;

% get q value for each flux surface
psigr=linspace(yyy.simag,yyy.sibry,yyy.nw);
for ij = 1:length(yyy.flx_srfs);
  yyy.flx_srfs(ij).v;
  yyy.flx_srfs(ij).q=interp1(psigr,yyy.qpsi,yyy.flx_srfs(ij).v)  ; 
end

if (optplot>1)
    figure
    for ij = 1:10:length(yyy.flx_srfs);
        plot(yyy.flx_srfs(ij).r, yyy.flx_srfs(ij).z, 'r-');
        hold on;
    end
    axis equal
end

