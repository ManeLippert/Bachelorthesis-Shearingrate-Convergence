function outcffsa = read_gkwcffsa(sim_name, pth);
%
% function outcffsa = read_gkwcffsa(sim_name, pth);
%
% Note about output:
%
% e_0, e_1, e_2, e_3, e_4, e_eta = averages included in Memo
% Vcal_spgr = Vcal with each species parameters
% Vcal_iongr_uprim_only = Vcal with uprim from D and each species parameters for the other terms
%

if ~exist('pth');
[a pth] = unix('echo $GKWMAT_DIR');
end
if isempty(pth);
clear pth
[a pth] = unix('echo $GKWMAT_DIR');
end
pthsim = [pth(1:end-1) '/other/'];

clear flnm
flnm = [pthsim sim_name '/cffsa.dat'];
fptr = fopen(flnm,'r');
if fptr == -1
error(['No such simulation in ' pth]);
disp(['These are the simulations available']);
[a avsims] = unix(['ls ' pthsim]);
disp(avsims)
ystr = -1;
return;
end;
fclose(fptr);

clear outcffsa;

%fptr = fopen(flnm); 
%frewind(fptr);
%aaa = fscanf(fptr, '%s',1);
%while ~isempty(aaa);
%outcffsa.(aaa) = fscanf(fptr, '%f');
%aaa = fscanf(fptr, '%s',1);
%end;

fptr = fopen(flnm);
frewind(fptr);
aaa = fscanf(fptr, '%s',1);
while ~isempty(aaa);
  strmem = aaa;
  fpos=ftell(fptr);
  [zzz, nel] = fscanf(fptr, '%f');
  fseek(fptr,fpos,'bof');
  outcffsa.(aaa) = fscanf(fptr, '%f',nel)';
  aaa = fscanf(fptr, '%s',1);
end;
fclose(fptr);

