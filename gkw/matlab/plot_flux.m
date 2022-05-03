% Plots the fluxes of a certain structure 
%
% Usage : plot_flux(in) 
% where 'in' is the structure that contains the fluxes 

function plot_flux(in) 

  % Get the screensize and generate the figure 
  scrsz = get(0,'ScreenSize');
  figure('Position',[1 3*scrsz(4)/4 3*scrsz(3)/4 scrsz(4)/2])

  % The symbols of the lines use (append if necessary)
  symbol = ['r-' 'b-' 'g-' 'm-' 'k-' 'r--' 'b--' 'g--' 'm--' 'k--' 'r.-' 'b.-' 'g.-' 'm.-' 'k.-' ]; 
  
  % Plot the heat fluxes (Electro-static is always assumed to exist)
  subplot(2,2,1); 
  nspec = size(in.qfluxes,2); 
  for i = 1: nspec 
    plot(in.time(:,1),in.qfluxes(:,i),symbol(2*i-1:2*i)); 
    npnt = size(in.qfluxes,1); 
    label = ['ES' int2str(i)]; 
    text(in.time(npnt,1),in.qfluxes(npnt,i),label); 
    hold on;
  end 
  % Would be easier in future just to check for existance of file
  if isfield(in.CONTROL,'nlapar'); 
    if (in.CONTROL.nlapar(1)==1) 
      for i = 1: nspec 
        plot(in.time(:,1),in.qfluxem(:,i),symbol(2*i-1+2*nspec:2*i+2*nspec)); 
        npnt = size(in.qfluxem,1); 
        label = ['EM' int2str(i)]; 
        text(in.time(npnt,1),in.qfluxem(npnt,i),label); 
      end 
    end 
  end 
  if isfield(in.CONTROL,'nlbpar'); 
    if (in.CONTROL.nlbpar(1)==1) 
      for i = 1: nspec 
        plot(in.time(:,1),in.qfluxbpar(:,i),symbol(2*i-1+4*nspec:2*i+4*nspec)); 
        npnt = size(in.qfluxbpar,1); 
        label = ['BP' int2str(i)]; 
        text(in.time(npnt,1),in.qfluxbpar(npnt,i),label); 
      end 
    end    
  end 
  title('Heat fluxes'); 
  xlabel('t v_{thref} / R_{ref}'); 
  ylabel('Q')


  % Plot the particle fluxes (Electro-static is always assumed to exist)
  subplot(2,2,2)
  nspec = size(in.pfluxes,2); 
  for i = 1: nspec 
    plot(in.time(:,1),in.pfluxes(:,i),symbol(2*i-1:2*i)); 
    npnt = size(in.pfluxes,1); 
    label = ['ES' int2str(i)]; 
    text(in.time(npnt,1),in.pfluxes(npnt,i),label); 
    hold on;
  end 
  if isfield(in.CONTROL,'nlapar'); 
    if (in.CONTROL.nlapar(1)==1) 
      for i = 1: nspec 
        plot(in.time(:,1),in.pfluxem(:,i),symbol(2*i-1+2*nspec:2*i+2*nspec)); 
        npnt = size(in.pfluxem,1); 
        label = ['EM' int2str(i)]; 
        text(in.time(npnt,1),in.pfluxem(npnt,i),label); 
      end 
    end 
  end 
  if isfield(in.CONTROL,'nlbpar'); 
    if (in.CONTROL.nlbpar(1)==1) 
      for i = 1: nspec 
        plot(in.time(:,1),in.pfluxbpar(:,i),symbol(2*i-1+4*nspec:2*i+4*nspec)); 
        npnt = size(in.pfluxbpar,1); 
        label = ['BP' int2str(i)]; 
        text(in.time(npnt,1),in.pfluxbpar(npnt,i),label); 
      end 
    end    
  end 
  title('Particle fluxes'); 
  xlabel('t v_{thref} / R_{ref}'); 
  ylabel('\Gamma')



  % Plot the toroidal momentum fluxes (Electro-static is always assumed to exist)
  subplot(2,2,3)
  nspec = size(in.vfluxes,2); 
  for i = 1: nspec 
    plot(in.time(:,1),in.vfluxes(:,i),symbol(2*i-1:2*i)); 
    npnt = size(in.vfluxes,1); 
    label = ['ES' int2str(i)]; 
    text(in.time(npnt,1),in.vfluxes(npnt,i),label); 
    hold on;
  end 
  if isfield(in.CONTROL,'nlapar'); 
    if (in.CONTROL.nlapar(1)==1) 
      for i = 1: nspec 
        plot(in.time(:,1),in.vfluxem(:,i),symbol(2*i-1+2*nspec:2*i+2*nspec)); 
        npnt = size(in.vfluxem,1); 
        label = ['EM' int2str(i)]; 
        text(in.time(npnt,1),in.vfluxem(npnt,i),label); 
      end 
    end 
  end 
  if isfield(in.CONTROL,'nlbpar'); 
    if (in.CONTROL.nlbpar(1)==1) 
      for i = 1: nspec 
        plot(in.time(:,1),in.vfluxbpar(:,i),symbol(2*i-1+4*nspec:2*i+4*nspec)); 
        npnt = size(in.vfluxbpar,1); 
        label = ['BP' int2str(i)]; 
        text(in.time(npnt,1),in.vfluxbpar(npnt,i),label); 
      end 
    end    
  end 
  title('Momentum fluxes'); 
  xlabel('t v_{thref} / R_{ref}'); 
  ylabel('\Gamma_\phi')
