function [out, outnosqueeze, outpar] = read_gkwparallel(sim_name,proj,eiv);
%
% read the gkw parallel.dat into a structure 
% with data reshaped into multidimensional arrays
%
% function [out, outnosqueeze, outpar] = read_gkwparallel(sim_name, proj,eiv);
%
% out is the data of parallel.dat reshaped into 2d arrays when possible
% outnosqueeze is the data of parallel.dat reshaped into 4d arrays
% outpar is the raw data of parallel.dat not reshaped
%
% eiv is an optional integer specifying which eigenfunction 
% to read in case of runs using the eigenvalue solver
%
% Uses folder structure of gkwnlin and requires gkwpath to be set
%
% CLA with modifications by FJC, 11.11.13
%

if ~exist('eiv','var')
  pthsim = gkwpath('parallel',proj);
  flnm = [pthsim sim_name]
else
  pthsim = gkwpath('other',proj)
  flnm = [pthsim sim_name '/parallel' num2str(eiv,'%03i') '.dat']
end

try
  strin = read_gkwinput(sim_name, proj);
catch
  disp(['No such simulation in project' proj]);
  disp(['These are the simulations available']);
  list_gkwscan(proj)  
end
    
fptr = fopen(flnm,'r');
if fptr == -1
    disp(['No such output' flnm]);  
    disp(['These are the runs available']);  
    [a avsims] = unix(['ls ' pthsim]);
    disp(avsims)
    ystr = -1;
    return;
end;
fclose(fptr);

aaa=load(flnm);
nbspecies = strin.GRIDSIZE.number_of_species;

outpar.s_grid = aaa(1:end/nbspecies,1);
outpar.PhiR = aaa(1:end/nbspecies,2);
outpar.PhiI = aaa(1:end/nbspecies,3);
outpar.AparR = aaa(1:end/nbspecies,4);
outpar.AparI = aaa(1:end/nbspecies,5);
for ij = 1:nbspecies;
    outpar.densR(:,ij)  = aaa((ij-1)*end./nbspecies+1 : ij*end/nbspecies,6);
    outpar.densI(:,ij)  = aaa((ij-1)*end./nbspecies+1:ij*end/nbspecies,7);
    outpar.TparR(:,ij)  = aaa((ij-1)*end./nbspecies+1:ij*end/nbspecies,8);
    outpar.TparI(:,ij)  = aaa((ij-1)*end./nbspecies+1:ij*end/nbspecies,9);
    outpar.TperpR(:,ij) = aaa((ij-1)*end./nbspecies+1:ij*end/nbspecies,10);
    outpar.TperpI(:,ij) = aaa((ij-1)*end./nbspecies+1:ij*end/nbspecies,11);
    outpar.UparR(:,ij)  = aaa((ij-1)*end./nbspecies+1:ij*end/nbspecies,12);
    outpar.UparI(:,ij)  = aaa((ij-1)*end./nbspecies+1:ij*end/nbspecies,13);
end
outpar.BparR  = aaa(1:end/nbspecies,14);
outpar.BparI  = aaa(1:end/nbspecies,15);
outpar.Phisq = outpar.PhiR.^2 + outpar.PhiI.^2;
outpar.input = strin;

n_s_grid = outpar.input.GRIDSIZE.n_s_grid; out.n_s_grid = n_s_grid;
nmod = outpar.input.GRIDSIZE.nmod; out.nmod = nmod;
nx = outpar.input.GRIDSIZE.nx; out.nx = nx;
%
out.s_grid = permute(reshape(transpose(outpar.s_grid), n_s_grid, nx, nmod), [3 2 1]);
out.PhiR = permute(reshape(transpose(outpar.PhiR), n_s_grid, nx, nmod), [3 2 1]);
out.PhiI = permute(reshape(transpose(outpar.PhiI), n_s_grid, nx, nmod), [3 2 1]);
out.AparR = permute(reshape(transpose(outpar.AparR), n_s_grid, nx, nmod), [3 2 1]);
out.AparI = permute(reshape(transpose(outpar.AparI), n_s_grid, nx, nmod), [3 2 1]);
for js = 1:nbspecies;
    out.densR(js, :, :, :) = permute(reshape(transpose(outpar.densR(:,js)), n_s_grid, nx, nmod), [3 2 1]);
    out.densI(js, :, :, :) = permute(reshape(transpose(outpar.densI(:,js)), n_s_grid, nx, nmod), [3 2 1]);
    out.TparR(js, :, :, :) = permute(reshape(transpose(outpar.TparR(:,js)), n_s_grid, nx, nmod), [3 2 1]);
    out.TparI(js, :, :, :) = permute(reshape(transpose(outpar.TparI(:,js)), n_s_grid, nx, nmod), [3 2 1]);
    out.TperpR(js, :, :, :) = permute(reshape(transpose(outpar.TperpR(:,js)), n_s_grid, nx, nmod), [3 2 1]);
    out.TperpI(js, :, :, :) = permute(reshape(transpose(outpar.TperpI(:,js)), n_s_grid, nx, nmod), [3 2 1]);
    out.UparR(js, :, :, :) = permute(reshape(transpose(outpar.UparR(:,js)), n_s_grid, nx, nmod), [3 2 1]);
    out.UparI(js, :, :, :) = permute(reshape(transpose(outpar.UparI(:,js)), n_s_grid, nx, nmod), [3 2 1]);
end
out.BparR = permute(reshape(transpose(outpar.BparR), n_s_grid, nx, nmod), [3 2 1]);
out.BparI = permute(reshape(transpose(outpar.BparI), n_s_grid, nx, nmod), [3 2 1]);
out.Phisq = permute(reshape(transpose(outpar.Phisq), n_s_grid, nx, nmod), [3 2 1]);

outnosqueeze = out;

out.s_grid = squeeze(out.s_grid);
out.PhiR   = squeeze(out.PhiR);
out.PhiI   = squeeze(out.PhiI);
out.AparR  = squeeze(out.AparR);
out.AparI  = squeeze(out.AparI);
out.densR  = squeeze(out.densR);
out.densI  = squeeze(out.densI);
out.TparR  = squeeze(out.TparR);
out.TparI  = squeeze(out.TparI);
out.TperpR = squeeze(out.TperpR);
out.TperpI = squeeze(out.TperpI);
out.UparR  = squeeze(out.UparR);
out.UparI  = squeeze(out.UparI);
out.BparR  = squeeze(out.BparR);
out.BparI  = squeeze(out.BparI);
out.Phisq  = squeeze(out.Phisq);

return
