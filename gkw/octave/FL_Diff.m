%*********************************************************************************************************
%This program calculates the diffusion coefficient for a single fieldline and averages over all fieldlines
%The output is a plot of the average against the number of turns.
%It was written by Benjamin Schobert.
%*********************************************************************************************************

%Inputs:                Poincare1.mat (as created by FL_Tracing.m)

%Outputs:               Diffusion1.mat (contains all diffusion coefficients and the average)
%                       figure(1) (plot of average against number of turns)


clear;
clf(figure(2));


%Set to 1 for saving all plots as eps files with latex code
save_plot = 0;


poincare = load('Poincare1.mat');
ly = poincare.ly;
diff_psi_n = poincare.diff_psi_n;               %This array contains all differences for psi for each turn
diff_zeta_n = poincare.diff_zeta_n.*ly;
size = size(diff_psi_n);
num_fl = size(1);                               %number of field lines
nturns = size(2);                               %number of toriodal turns

psi_dico_n(num_fl,nturns) = 0;                  %stores all the diffusion coefficients
zeta_dico_n(num_fl,nturns) = 0;
psi_dico_av(nturns) = 0;                        %stores the average of the diffusion coefficients
zeta_dico_av(nturns) = 0;

psi_dico_n = (diff_psi_n .* diff_psi_n);
zeta_dico_n = (diff_zeta_n .* diff_zeta_n);

for k=1:nturns
  psi_dico_n(:,k) = psi_dico_n(:,k)/(k);
  zeta_dico_n(:,k) = zeta_dico_n(:,k)/(k);
end

psi_dico_av = mean(psi_dico_n);                 %average for each turn
zeta_dico_av = mean(zeta_dico_n);

figure(2);                                      %plot the average
hold on;
plot(psi_dico_av,'r');
plot(zeta_dico_av,'k');
%plot(zeta_dico_n([55,72],:)');
%plot(zeta_dico_n');
%plot(psi_dico_n','k');
legend ({'$\psi$', '$\zeta$'});
xlabel('number of turns')
ylabel('normalized diffusion coefficient');

save('Diffusion1.mat','psi_dico_n','psi_dico_av','zeta_dico_n','zeta_dico_av');
if(save_plot)
  print('epslatex', 'figure2.tex');               %create eps file with plot and latex file with labels
endif

figure(4);                                      %plot the average
hold on;
plot(psi_dico_av,'r');
%plot(zeta_dico_av,'k');
%plot(zeta_dico_n([55,72],:)');
%plot(zeta_dico_n');
%plot(psi_dico_n','k');
legend ({'$\psi$'});
xlabel('number of turns')
ylabel('normalized diffusion coefficient');

save('Diffusion1.mat','psi_dico_n','psi_dico_av','zeta_dico_n','zeta_dico_av');
if(save_plot)
  print('epslatex', 'figure4.tex');               %create eps file with plot and latex file with labels
endif