#!/usr/bin/env octave

% This is a standalone GNU Octave script.

try
libgkwoct;
catch
disp('')
disp('ERROR: Could not source libgkwoct. Please put the line ''addpath("~/gkw/octave");'' (with correct path to the $GKW_HOME/octave folder) into the file ~/.octaverc .')
disp('')
end

graphics_toolkit("gnuplot")
set(0, 'defaultfigurevisible', 'off');

set(0,'defaulttextinterpreter','Tex')
set(0, 'defaultlinelinewidth', 6);
%set(0, 'defaultaxesygrid','on');


pkg load hdf5oct


krho = h5read('gkwdata.h5','/grid/krho')(1,:);
number_of_species = h5readatt('gkwdata.h5','/','grid.number_of_species');
degr_grid = h5read('gkwdata.h5','/diagnostic/diagnos_cross_phase/cross_phase_degr_grid');

% the same tokens as in diagnos_generic
field_tokens = {'phi','Apar','Bpar'}
% 'fdis' not needed
moment_tokens= {'phi_ga','Apar_ga','Bpar_ga','dens','T','vpar','vparsq','vperpsq','Tpar','Tperp','dens_ps', 'dens_tr', 'dens_ga', 'dens_polar'}

for tok1 = [field_tokens, moment_tokens]
  is_field1 = find(ismember(field_tokens, tok1{1}))
  if(is_field1)
    m1 = 1;
  else
    m1 = number_of_species;
  end
  for im1 = 1:m1
    if(m1 == 1)
      species1 = '';
    else
      species1 = sprintf("_sp%02i", im1);
    end

    for tok2 = [field_tokens, moment_tokens]
      is_field2 = find(ismember(field_tokens, tok2{1}))
      if(is_field2)
        m2 = 1;
      else
        m2 = number_of_species;
      end
      for im2 = 1:m2

        if(m2 == 1)
          species2 = '';
        else
          species2 = sprintf("_sp%02i", im2);
        end
        %sprintf('%s%4.4d%s%4.4d','PhiTD',i,'_',j));

        eim_number = 0;
        while(true)
          if(eim_number == 0)
            eim = '';
          else
            eim = sprintf('_eim%01i', eim_number);
          end
          luname = sprintf('cross_phase_pdf_%s%s_%s%s%s', tok1{1}, species1, tok2{1}, species2, eim)

          try
            data = h5read('gkwdata.h5', ['/diagnostic/diagnos_cross_phase/', luname]);
          catch
            break
          end
          % test if axes labelling is correct:
          %data(1,:) = 0.2;
          %data(:,1) = 0.2;
          %data(:,100) = 1;
          
          % want a plot with degrees on the x axis and modes on the y axis
          data = transpose(data);

          figure()

          ## imagesc(degr_grid,krho,data)
          ## colormap(flip(bone(200)))
          ## colorbar();
          hold on
          hold all

          contour(degr_grid,krho(2:end),data(2:end,:),'k')

          set(gca,'YDir','normal')
          set(gca(),'xtick',-180:45:180)
          xlabel("deg.")
          ylabel("k_\\theta\\rho")
          title(["weighted PDF of phase difference ",tok1{1},' vs. ',tok2{1}] )

          % highlight the maximum, note that the data was transposed!
          [m, imax] = max(data,[],2);
          plot(degr_grid(imax),krho,'ro')
          
          grid minor on
          print_all_formats(luname)
          hold off

          % check how well it sums to 1:

          cross_phase_nclasses = size(data,1);
          delta_angle = pi / (cross_phase_nclasses - 1);
          s = 0;
          for i = 1:cross_phase_nclasses
            if(i == 1 || i == cross_phase_nclasses)
              s += data(:,i) * delta_angle*0.5;
            else
              s += data(:,i) * delta_angle;
            end
          end
          disp("Sum to zero:")
          s

          eim_number++;
        end

      end
    end
  end
end
