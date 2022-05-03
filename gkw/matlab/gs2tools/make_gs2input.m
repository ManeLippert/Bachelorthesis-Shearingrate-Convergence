function y = make_gs2input(filename, xpareq, xkt, nbspecies, xspc, xprkn, pth, reffile);
%
% function y = make_gs2input(filename,xpareq,xkt,nbspecies,xspc,xprkn,[pth], [reffile])
%
% You can create xpareq, xkt, xspc and xprkn structures by editing
% your local file gs2_standard_set.m, in which you can save the input values
% or by creating an "experimental" standard set with function
% read_arp_2_stset from Peeters' idl GUI output in a matlab gs2 standard set.
%
% To create all the structures, before running the function make_gs2input,
% you run first the matlab script
% >> gs2_standard_set
% or
% >> stset_17170_3_5  (here 17170 is the shot number and 3_5 = 3.5 s the time)
%
% These matlab scripts prepare your matlab workspace, then you can call
% y = make_gs2input(filename, xpareq, xkt, nbspecies, xspc, xprkn)
%
% If you want to give a specific reference file to give different
% input parameters which are not included in the standard set structures
% (for instance nstep, delt, nperiod, phiinit or any other), you
% can create your reference file in your /u/ 'usr' /gs2run/input
% directory, and give the name of this reference file  at the
% variable 'reffile' in input of this function.
% Please give separately both the path and the name of the
% reference file in input in the two input variables [pth]
% and [reffile].
%
% GS2 input file is saved in directory /u/ 'usr'/gs2run/input
%
% CLA 10.08.03 / Update for Reference file in input 20.07.04
%

%[a zhomedir] = unix('echo $HOME');
%homedir = zhomedir(1:end-1)
%%
%if ~exist('pth')
%pth = [homedir '/gs2run/input/'];
%end
%if isempty(pth)
%pth = [homedir '/gs2run/input/'];
%end

usrnm = find_usrnm;
%
if ~exist('pth')
    pth = ['/afs/ipp-garching.mpg.de/home/' usrnm(1) '/' usrnm '/gs2run/input/'];
end
if isempty(pth);
    pth = ['/afs/ipp-garching.mpg.de/home/' usrnm(1) '/' usrnm '/gs2run/input/'];
end

if ~exist('reffile')
    reffile = 'abcde';
end
if isempty(reffile)
    reffile = 'abcde';
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
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
parsp{8}  = 'vnewk';
%
comsp{1}='  ! Z_s charge of species s';
comsp{2}='  ! mass_s relative to reference mass (i.e. H or D)';
comsp{3}='  ! n_s/n_e';
comsp{4}='  ! T_s/T_ref';
comsp{5}='  ! a/LT_s';
comsp{6}='  ! a/L_n_s';
comsp{7}='  ! u_parallel gradient';
comsp{8}='  ! collisionality for this species, see below';
%
%
flnm = filename;
pthflnm = [pth flnm];
if strcmp(reffile, 'abcde');
    %fprd=fopen([homedir '/gs2run/input/hmsc_start'],'r');
    %fprd=fopen(['/afs/ipp-garching.mpg.de/home/' usrnm(1) '/' usrnm '/gs2run/input/hmsc_start'],'r');
    %fprd=fopen([ '/gs2run/input/hmsc_start'],'r');
    fprd=fopen('~cla/gs2run/input/waltz_ref','r');
else
    fprd=fopen([pth reffile],'r');
end
%
fptr = fopen(pthflnm, 'w');

frewind(fprd);
frewind(fptr);

yk = 0;
while yk == 0
    sss = fgetl(fprd);
    fprintf(fptr,'%s\n',sss);
    if strcmp(sss(1:end), '&theta_grid_knobs'), yk = 1; end
end
sss = fgetl(fprd);
fprintf(fptr,'%s\n',sss);
sss = fgetl(fprd);
if (isfield(xpareq,'eps')) %% s-alpha equilibrium input structure
    if (strncmp(sss, 'equil',5)); %% reference file with s-alpha equilibrium
        %% leaves this way
        fprintf(fptr,'%s\n',sss);
        ss1 = fgetl(fprd);
        fprintf(fptr,'%s\n',ss1);
    else %% reference file with eik equilibrium
        fprintf(fptr,'%s\n',sss(2:end));
        ss1 = fgetl(fprd);
        fprintf(fptr,'%s\n',['!' ss1]);
    end
elseif (isfield(xpareq,'rhoc')) %% eik equilibrium input structure
    if (strncmp(sss, 'equil',5)); %% reference file with s-alpha equilibrium
        %% comment this
        fprintf(fptr,'%s\n',['!' sss]);
        ss1 = fgetl(fprd);
        fprintf(fptr,'%s\n',ss1(2:end));
    else %% reference file with eik equilibrium
        %% leave in this form
        fprintf(fptr,'%s\n',sss);
        ss1 = fgetl(fprd);
        fprintf(fptr,'%s\n',ss1);
    end
end
%
yk = 0;
while yk == 0
    sss = fgetl(fprd);
    fprintf(fptr,'%s\n',sss);
    if strcmp(sss(1:end), '&theta_grid_parameters'), yk = 1; end
end
%
for jj = 1:4
    sss = fgetl(fprd);
    fprintf(fptr,'%s\n',sss);
end
%
%% Writes equilibrium data !!!!
%
if isfield(xpareq, 'eps');
    % s-alpha equilibrium structure
    sss = fgetl(fprd);
    is = 1; while(~strcmp(sss(is),'!')), is=is+1; end;
    fprintf(fptr,'%s',' eps =');
    fprintf(fptr,' %f', xpareq.eps);
    fprintf(fptr,'%s\n', sss(is-6: end));
    %
    sss = fgetl(fprd);
    fprintf(fptr,'%s\n',sss);
    %
    sss = fgetl(fprd);
    is = 1; while(~strcmp(sss(is),'!')), is=is+1; end;
    fprintf(fptr,'%s',' epsl =');
    fprintf(fptr,' %f', xpareq.epsl);
    fprintf(fptr,'%s\n', sss(is-6: end));
    %
    sss = fgetl(fprd);
    is = 1; while(~strcmp(sss(is),'!')), is=is+1; end;
    fprintf(fptr,'%s',' pk =');
    fprintf(fptr,' %f', xpareq.pk);
    fprintf(fptr,'%s\n', sss(is-4: end));
    %
    sss = fgetl(fprd);
    is = 1; while(~strcmp(sss(is),'!')), is=is+1; end;
    fprintf(fptr,'%s',' shat =');
    fprintf(fptr,' %f', xpareq.shat);
    fprintf(fptr,'%s\n', sss(is-5: end));
    %
    sss = fgetl(fprd);
    is = 1; while(~strcmp(sss(is),'!')), is=is+1; end;
    fprintf(fptr,'%s',' shift =');
    fprintf(fptr,' %f', xpareq.shift);
    fprintf(fptr,'%s\n', sss(is-6: end));
    %
    %%%%
    %
elseif isfield(xpareq, 'rhoc');
    % eik equilibrium structure
    yk = 0;
    while yk == 0
        sss = fgetl(fprd);
        fprintf(fptr,'%s\n',sss);
        if strncmp(sss, '! If equilibrium_option=''eik''',29), yk = 1; end
    end
    sss = fgetl(fprd);
    fprintf(fptr,'%s\n',sss);
    sss = fgetl(fprd);
    fprintf(fptr,'%s\n',sss);
    sss = fgetl(fprd);
    fprintf(fptr,'%s\n',sss);
    %
    sss = fgetl(fprd);
    is = 1; while(~strcmp(sss(is),'!')), is=is+1; end;
    fprintf(fptr,'%s',' rhoc =');
    fprintf(fptr,' %f', xpareq.rhoc);
    fprintf(fptr,'%s\n', sss(is-2: end));
    %
    sss = fgetl(fprd);
    is = 1; while(~strcmp(sss(is),'!')), is=is+1; end;
    fprintf(fptr,'%s',' Rmaj =');
    fprintf(fptr,' %f', xpareq.Rmaj);
    fprintf(fptr,'%s\n', sss(is-2: end));
    %
    sss = fgetl(fprd);
    is = 1; while(~strcmp(sss(is),'!')), is=is+1; end;
    fprintf(fptr,'%s',' R_geo=');
    fprintf(fptr,' %f', xpareq.R_geo);
    fprintf(fptr,'%s\n', sss(is-2: end));
    %
    sss = fgetl(fprd);
    fprintf(fptr,'%s\n',sss);
    sss = fgetl(fprd);
    fprintf(fptr,'%s\n',sss);
    %
    sss = fgetl(fprd);
    is = 1; while(~strcmp(sss(is),'!')), is=is+1; end;
    fprintf(fptr,'%s',' qinp =');
    fprintf(fptr,' %f', xpareq.qinp);
    fprintf(fptr,'%s\n', sss(is-4: end));
    %
    sss = fgetl(fprd);
    is = 1; while(~strcmp(sss(is),'!')), is=is+1; end;
    igf = 0;
    if is <= 5; %%% '! shift' in reference input file for s-alpha equilibrium
        igf = 1;
        while(~strcmp(sss(is+5),'!')), is=is+1; end;
        is=is-5;
    end
    fprintf(fptr,'%s',' shift =');
    fprintf(fptr,' %f', xpareq.shift);
    fprintf(fptr,'%s\n', sss(is-2+igf*10: end));
    %
    sss = fgetl(fprd);
    is = 1; while(~strcmp(sss(is),'!')), is=is+1; end;
    fprintf(fptr,'%s',' akappa =');
    fprintf(fptr,' %f', xpareq.akappa);
    fprintf(fptr,'%s\n', sss(is-2: end));
    %
    sss = fgetl(fprd);
    is = 1; while(~strcmp(sss(is),'!')), is=is+1; end;
    fprintf(fptr,'%s',' akappri =');
    fprintf(fptr,' %f', xpareq.akappri);
    fprintf(fptr,'%s\n', sss(is-2: end));
    %
    sss = fgetl(fprd);
    is = 1; while(~strcmp(sss(is),'!')), is=is+1; end;
    fprintf(fptr,'%s',' tri =');
    fprintf(fptr,' %f', xpareq.tri);
    fprintf(fptr,'%s\n', sss(is-2: end));
    %
    sss = fgetl(fprd);
    fprintf(fptr,'%s\n',sss);
    %
    sss = fgetl(fprd);
    is = 1; while(~strcmp(sss(is),'!')), is=is+1; end;
    fprintf(fptr,'%s',' tripri =');
    fprintf(fptr,' %f', xpareq.tripri);
    fprintf(fptr,'%s\n', sss(is-2: end));
end
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
if isfield(xprkn,'nperiod');
    yk = 0;
    while yk == 0
        sss = fgetl(fprd);
        fprintf(fptr,'%s\n',sss);
        if strncmp(sss, ' ntheta=',8), yk = 1; end
    end
    sss = fgetl(fprd);
    is = 1; while(~strcmp(sss(is),'!')), is=is+1; end;
    fprintf(fptr,'%s',' nperiod=');
    fprintf(fptr,' %i', xprkn.nperiod);
    fprintf(fptr,'   %s\n', sss(is: end));
end
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
yk = 0;
while yk == 0
    sss = fgetl(fprd);
    fprintf(fptr,'%s\n',sss);
    if strcmp(sss, '&parameters'), yk = 1; end
end
sss = fgetl(fprd);
is = 1; while(~strcmp(sss(is),'!')), is=is+1; end;
fprintf(fptr,'%s',' beta=');
fprintf(fptr,' %f', xprkn.beta);
fprintf(fptr,'%s\n', sss(is-3: end));
for jj = 1:15
    sss = fgetl(fprd);
    fprintf(fptr,'%s\n',sss);
end
sss = fgetl(fprd);
is = 1; while(~strcmp(sss(is),'!')), is=is+1; end;
fprintf(fptr,'%s',' zeff =');
fprintf(fptr,' %f', xprkn.zeff);
fprintf(fptr,'%s\n', sss(is-3: end));
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%
%%%% &theta_grid_eik_knobs     ! Only read if equilibrium_option='eik'
%
if (isfield(xpareq,'rhoc')) %% eik equilibrium input structure
    yk = 0;
    while yk == 0
        sss = fgetl(fprd);
        fprintf(fptr,'%s\n',sss);
        if strncmp(sss, '&theta_grid_eik_knobs',21), yk = 1; end
    end
    for ij = 1:30
        sss = fgetl(fprd);
        fprintf(fptr,'%s\n',sss);
    end
    sss = fgetl(fprd);
    is = 1; while(~strcmp(sss(is),'!')), is=is+1; end;
    fprintf(fptr,'%s',' s_hat_input = ');
    fprintf(fptr,' %f', xpareq.s_hat_input);
    fprintf(fptr,'%s\n', sss(is-2: end));
    sss = fgetl(fprd);
    is = 1; while(~strcmp(sss(is),'!')), is=is+1; end;
    fprintf(fptr,'%s',' beta_prime_input = ');
    fprintf(fptr,' %f', xpareq.beta_prime_input);
    fprintf(fptr,'%s\n', sss(is-2: end));
    
end

yk = 0;
while yk == 0
    sss = fgetl(fprd);
    fprintf(fptr,'%s\n',sss);
    if strcmp(sss(1:end), '&kt_grids_range_parameters'), yk = 1; end
end
%
sss = fgetl(fprd);
is = 1; while(~strcmp(sss(is),'!')), is=is+1; end;
fprintf(fptr,'%s',' naky=');
fprintf(fptr,' %i', xkt.naky);
fprintf(fptr,'%s\n', sss(is-5: end));
sss = fgetl(fprd);
is = 1; while(~strcmp(sss(is),'!')), is=is+1; end;
fprintf(fptr,'%s',' ntheta0=');
fprintf(fptr,' %i', xkt.ntheta0);
fprintf(fptr,'%s\n', sss(is-3: end));
sss = fgetl(fprd);
fprintf(fptr,'%s',' aky_min=');
fprintf(fptr,' %f\n', xkt.aky_min);
sss = fgetl(fprd);
fprintf(fptr,'%s',' aky_max=');
fprintf(fptr,' %f\n', xkt.aky_max);
sss = fgetl(fprd);
fprintf(fptr,'%s',' theta0_min=');
fprintf(fptr,' %f\n', xkt.theta0_max);
sss = fgetl(fprd);
fprintf(fptr,'%s',' theta0_max=');
fprintf(fptr,' %f\n', xkt.theta0_max);
%
%%%%%%%%%%%%%%%%%%%%%
%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
yk = 0;
while yk == 0
    sss = fgetl(fprd);
    fprintf(fptr,'%s\n',sss);
    if strcmp(sss, '&knobs'), yk = 1; end
end
sss = fgetl(fprd);
is = 1; while(~strcmp(sss(is),'!')), is=is+1; end;
fprintf(fptr,'%s',' fphi=');
fprintf(fptr,' %f', xprkn.fphi);
fprintf(fptr,'%s\n', sss(is-3: end));
%
sss = fgetl(fprd);
is = 1; while(~strcmp(sss(is),'!')), is=is+1; end;
fprintf(fptr,'%s',' fapar=');
fprintf(fptr,' %f', xprkn.fapar);
fprintf(fptr,'%s\n', sss(is-3: end));
%
sss = fgetl(fprd);
is = 1; while(~strcmp(sss(is),'!')), is=is+1; end;
fprintf(fptr,'%s',' faperp=');
fprintf(fptr,' %f', xprkn.faperp);
fprintf(fptr,'%s\n', sss(is-3: end));
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
yk = 0;
while yk == 0
    sss = fgetl(fprd);
    fprintf(fptr,'%s\n',sss);
    if strcmp(sss(1:end), '&species_knobs'), yk = 1; end
end
%
sss = fgetl(fprd);
is = 1; while(~strcmp(sss(is),'!')), is=is+1; end;
fprintf(fptr,'%s',' nspec=');
fprintf(fptr,' %i', nbspecies);
fprintf(fptr,'%s\n', sss(is-3: end));
%
sss = fgetl(fprd);
fprintf(fptr,'%s\n',sss);
%
for js = 1:nbspecies
    %
    sss = fgetl(fprd);
    if (js > 1)
        while (~strncmp(sss, ['&species_parameters_' num2str(js)] ,21))
            fprintf(fptr,'%s\n',sss);
            sss = fgetl(fprd);
        end
    end
    %
    fprintf(fptr, '%s\n',['&species_parameters_' num2str(js)]);
    %
    for jp = 1:length(parsp)
        sss = fgetl(fprd);
        fprintf(fptr, ' %s',[parsp{jp} '= ']);
        %fprintf(fptr, '%s', sss(1:7));
        eval(['fprintf(fptr, '' %e'',[xspc(js).' parsp{jp} ']);']);
        fprintf(fptr, '%s\n',comsp{jp});
    end
end
%
for j = 1:2000
    sss = fgetl(fprd);
    if (sss == -1),
        break;
    else
        fprintf(fptr,'%s\n',sss);
    end
end
%
fclose(fptr);
fclose(fprd);
y = 1;
