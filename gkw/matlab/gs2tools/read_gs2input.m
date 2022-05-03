function [eql, spc, xprkn, xkt] = read_gs2input(flnm, pth)
%
% function [eql, spc, xprkn, xkt] = read_gs2input(flnm, [pth]);
% 
% written for gs2 input using s-alpha equilibria.
%
% Input. 
%       flnm :  filename
%       pth  :  path (optional: by default looks in user's directory ~/gs2run/input)
%       if pth == 1 looks in user's directory ~/gs2run/gs2archive/input)
% Output.
%              equlibrium_option s-alpha
%       eql  : structure with s-alpha equilibrium parameters 
%              (eps, epsl, pk, shat, shift)
%              eps = r/R, used in B = B_0 / (1 + eps * cos(theta))
%              epsl = 2.0 * a/R
%              pk = 2.0 * a/(q*R)
%              shat = r/q * (dq/dr)
%              shift = alpha in the s-alpha mode
%
%              equilibrium_option eik
%       eql  : structure with eik equilibrium parameters 
%              (rhoc, Rmaj, R_geo, qinp, shift, 
%               akappa, akappri, tri, tripri, 
%               s_hat_input, beta_prime_input)
%                rhoc:
%                Rmaj:
%               R_geo:
%                qinp:
%               shift:
%              akappa:
%             akappri:
%                 tri:
%              tripri:
%         s_hat_input:
%    beta_prime_input:
%
%
%       xkt    : structure with kt grid range paramters 
%               (naky, ntheta0, aky_min, aky_max, theta0_min, theta0_max)
%
%       spc   : array of nspec structures with species input parameters
%              z=  Z_s charge of species s
%              mass =  mass_s relative to reference mass (i.e. H or D)
%              dens = n_s/n_e
%              temp = T_s/T_ref
%              tprim = a/LT_s
%              fprim = a/L_n_s
%              uprim = u_parallel gradient
%              vnewk = collision frequency (see a gs2 input)
%
%
%       xprkn : other parameters and knobs
%               gs2 beta, zeff, and knobs fphi, fapar, faperp
%
%
% CLA 10.08.03
% CLA 10.10.03 Added xkt structure in output
% 

% usrnm = find_usrnm;
% %
% if ~exist('pth')
% pth = ['/afs/ipp-garching.mpg.de/home/' usrnm(1) '/' usrnm '/gs2run/input/'];
% %pth = ['/scratch/' usrnm '/GS2scratch/output/'];
% end
% if pth == 1
% pth = ['/afs/ipp-garching.mpg.de/home/' usrnm(1) '/' usrnm '/gs2run/gs2archive/input/'];
% end

pthflnm = [pth flnm];
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Sets names of variables it looks for
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Local equilibrium parameters
%% Only in case equilibrium option is s-alpha
sa_pareq{1} = 'eps';
sa_pareq{2} = 'epsl';
sa_pareq{3} = 'pk';
sa_pareq{4} = 'shat';
sa_pareq{5} = 'shift';
%% Only in case equilibrium option is eik
eik_pareq{1} = 'rhoc';
eik_pareq{2} = 'Rmaj';
eik_pareq{3} = 'R_geo';
eik_pareq{4} = 'qinp';
eik_pareq{5} = 'shift';
eik_pareq{6} = 'akappa';
eik_pareq{7} = 'akappri';
eik_pareq{8} = 'tri';
eik_pareq{9} = 'tripri';
%% additional &theta_grid_eik_knobs
eik_pareq{10} = 's_hat_input';
eik_pareq{11} = 'beta_prime_input';
%
%Parameters for each species
parsp{1}  = 'z';
parsp{2}  = 'mass';
parsp{3}  = 'dens';
parsp{4}  = 'temp';
parsp{5}  = 'tprim';
parsp{6}  = 'fprim';
parsp{7}  = 'uprim';
parsp{8}  = 'upara';
parsp{9}  = 'vnewk';
%
%Other parameters and knobs
parkn{1}  = 'beta';
parkn{2}  = 'zeff';
parkn{3}  = 'fphi';
parkn{4}  = 'fapar';
parkn{5}  = 'faperp';
%
%Parameters for ktheta range
parkt{1}  = 'naky';
parkt{2}  = 'ntheta0';
parkt{3}  = 'aky_min';
parkt{4}  = 'aky_max';
parkt{5}  = 'theta0_min';
parkt{6}  = 'theta0_max';

frdi = fopen(pthflnm, 'r');

frewind(frdi);
yk = 0;
while yk == 0
sss = fscanf(frdi,'%s', 1);
if strcmp(sss, '&theta_grid_knobs'), yk = 1; end;
end
yk = 0;
while yk == 0
sss = fscanf(frdi,'%s', 1);
if strncmp(sss, 'equilibrium_option=',19), yk = 1; end;
end
equil_option = sss(21:end-1);
%
if strcmp(equil_option, 's-alpha')
     pareq = sa_pareq;
elseif strcmp(equil_option, 'eik')
     pareq = eik_pareq;
%  disp('This matlab input-output interface is forseen only for equilibrium_option=s-alpha.')
%     disp('The equilibrium structure eql in output is not used in the actual GS2 input.')
%     disp('Structures on species parameters, aky values and knobs are valid.')
%%%  return  % runs unofficially also with 'eik' equilibrium_option
end
%
frewind(frdi);
yk = 0;
while yk == 0
sss = fscanf(frdi,'%s', 1);
if strcmp(sss, '&theta_grid_parameters'), yk = 1; end;
end
%
for je = 1:min(length(pareq),9)
yk = 0;
while yk == 0
sss = fscanf(frdi,'%s', 1);
%if strcmp(sss, pareq{je}), %% modified to include reading of R_geo= in 'eik'
if (strcmp(sss, [pareq{je} '=']) | strcmp(sss, [pareq{je}])),
yk = 1; end;
end
sss=fscanf(frdi,'%c', 2);
eval(['eql.' pareq{je} '= fscanf(frdi,''%e'', 1);']);
if sss(2) ~= '='
eval(['eql.' pareq{je} '= eql.' pareq{je} ' + ' sss ';']);
end
end
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
yk = 0;
while yk == 0
sss = fscanf(frdi,'%s', 1);
if strcmp(sss, '&parameters'), yk = 1; end
end
yk = 0;
while yk == 0
sss = fscanf(frdi,'%s', 1);
if (strcmp(sss, [parkn{1} '=']) | strcmp(sss, [parkn{1}])), yk = 1; end;
end
if strcmp(sss, [parkn{1} '='])
eval(['xprkn.' parkn{1} '= fscanf(frdi,''%e'', 1);']);
else
sss = fscanf(frdi,'%s', 1);
eval(['xprkn.' parkn{1} '= fscanf(frdi,''%e'', 1);']);
end
yk = 0;
while yk == 0
sss = fscanf(frdi,'%s', 1);
if (strcmp(sss, [parkn{2} '=']) | strcmp(sss, [parkn{2}])), yk = 1; end;
end
if strcmp(sss, [parkn{2} '='])
eval(['xprkn.' parkn{2} '= fscanf(frdi,''%e'', 1);']);
else
sss = fscanf(frdi,'%s', 1);
eval(['xprkn.' parkn{2} '= fscanf(frdi,''%e'', 1);']);
end
%
%%%% &theta_grid_eik_knobs     ! Only read if equilibrium_option='eik'
%
if strcmp(equil_option, 'eik')
yk = 0;
while yk == 0
sss = fscanf(frdi,'%s', 1);
if strcmp(sss, 'bishop>1,'), yk = 1; end
end
for je = 10:length(pareq)
yk = 0;
while yk == 0
sss = fscanf(frdi,'%s', 1);
%if strcmp(sss, pareq{je}), %% modified to include reading of R_geo= in 'eik'
if (strcmp(sss, [pareq{je} '=']) | strcmp(sss, [pareq{je}])),
yk = 1; end;
end
sss=fscanf(frdi,'%c', 2);
eval(['eql.' pareq{je} '= fscanf(frdi,''%e'', 1);']);
if sss(2) ~= '='
eval(['eql.' pareq{je} '= eql.' pareq{je} ' + ' sss ';']);
end
end
end
%
%%%%%
%
yk = 0;
while yk == 0
sss = fscanf(frdi,'%s', 1);
if strcmp(sss, '&knobs'), yk = 1; end
end
for jp = 3:length(parkn)
yk = 0;
while yk == 0
sss = fscanf(frdi,'%s', 1);
if strcmp(sss, [parkn{jp} '=']), yk = 1; end
end
eval(['xprkn.' parkn{jp} '= fscanf(frdi,''%e'', 1);']);
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% kt parameters
frewind(frdi)
yk = 0;
while yk == 0
sss = fscanf(frdi,'%s', 1);
if strcmp(sss, '&kt_grids_range_parameters'), yk = 1; end
end
for je = 1:length(parkt)
yk = 0;
while yk == 0
sss = fscanf(frdi,'%s', 1);
%if strcmp(sss, [parkt{je} '=']), yk = 1; end;
if strncmp(sss, [parkt{je} '='],length(parkt{je})+1), yk = 1; 
if length(sss) > length(parkt{je})+1; yk = 10; end;
end;end;
if yk == 1
eval(['xkt.' parkt{je} '= fscanf(frdi,''%e'', 1);']);
elseif yk == 10
eval(['xkt.' parkt{je} '= str2num(sss(length(parkt{je})+2:end));']);
end;end;
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% species parameters
yk = 0;
while yk == 0
sss = fscanf(frdi,'%s', 1);
if strcmp(sss, 'nspec='), yk = 1; end;
end
nbspecies = fscanf(frdi,'%e', 1);
%
for js = 1:nbspecies
%
%js
yk = 0;
while yk == 0
sss = fscanf(frdi,'%s', 1);
if strcmp(sss, ['&species_parameters_' num2str(js)]), yk = 1; end;
end
%
for jp = 1:length(parsp)
yk = 0;
%jp
while yk == 0
sss = fscanf(frdi,'%s', 1);
if strcmp(sss, [parsp{jp} '=']), yk = 1; end;
end
eval(['spc(js).' parsp{jp} '= fscanf(frdi,''%e'', 1);']);
end
%
end
fclose(frdi);
