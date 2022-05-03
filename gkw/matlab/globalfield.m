% function [thetazero, theta, r, q, s, x, y, ind, zeta, z2, rzc, izc, rzc2, izc2] = GetDataForPhi4(parallel, prof_back, nx, ns, rhostar, n, field = 1, write = '')
%
% ----------------------------------------------------------------------
%
% This function will extract from the data of a run the grids as well as
% the perturbed potential.
% The x,y coordinates are calculated from these, also the potential for
% constant phi.
%
% \attention This is an OCTAVE script, _not_ a matlab script, as matlab
%   does not know default parameters.
%
% Input parameters:
%   parallel:  Array created from loading the parallel.dat file.
%   prof_back: Array created from loading the prof_back.dat file.
%   nx:        Number of points in x direction.
%   ns:        Number of points in s (parallel) direction.
%   rhostar:   The value of rhostar used for the simulation.
%   n:         The value  of the modenumber used for the simulation.
%   field:     The field that should be extracted from the file.
%              1 = perturbed potential (default).
%              2 = || component of vector potential.
%              3 = perturbed density.
%              4 = perturbed parallel temperature.
%              5 = perturbed perpendicular temperature.
%              6 = perturbed parallel flow velocity.
%              7 = perturbed parallel (compressional) magnetic field.
%              \note As real and imaginary part are extracted the field
%                    number and the column number differ.
%   write:   A string, if this is not empty (the default) then is it used
%            as a name to save the fields and grids to.
%
%
% Output parameters:
%   \note : All arrays (fields/grids) that are returned, will have a size of ns times nx.
%   thetazero: The bloch shift parameter. Calculation doen't work yet.
%   theta:     Grid with the theta coordinate. 
%   r:         Grid with the r/epsilon/psi coordinate
%   q:         Array with the value of the safety factor for each point.
%   s:         Grid with the s-coordinate.
%   x:         Grid with the x-coordinate.
%   y:         Grid with the y-coordinate.
%   ind:       Grid with the indices that correspond to the radial direction.
%   zeta:      Grid  with the zeta-coordinate/difference.
%   z2:        Array with the rescaled absolute values of the field.
%   rzc:       Array with the real part of the perturbed field.
%   izc:       Array with the imaginary part of the perturbed field.
%   rzc2:      Array with the real part of the perturbed rescaled field.
%   izc2:      Array with the imgainary part of the perturbed rescaled field.
%
%   \note For working with (parallel.dat) files that contain multiple
%         species we do not use the complete data, but only the first
%         ns*nx entries, that correspond to the first species.
%   \todo Make used species selectable.
%   \todo Determine as many parameters as possible itself.
%   \todo Use a method for determining thetazero that works.
% ----------------------------------------------------------------------

function [thetazero, theta, r, q, s, x, y, ind, zeta, z2, rzc, izc, rzc2, izc2] = GetDataForPhi4(parallel, prof_back, nx, ns, rhostar, n, field = 1, write = '')
  kzeta = 2*pi*rhostar*n;

  r = prof_back(:,1);
  r = repmat(r.', ns, 1);% r-coordinate does not depend on s, so simply repeat the vector ns times.

  q = 0.854+17.037*r.^2;

  % Get the s-coordinate from the data. Also read the data, in this case it is the square of the value of the perturbed potential.
  s = reshape(parallel(1:ns*nx, 1), ns, nx);

  % Get the real/imaginary part of the perturbed potential.
  rzc = reshape(parallel(1:ns*nx, 2*field + 0), ns, nx);
  izc = reshape(parallel(1:ns*nx, 2*field + 1), ns, nx);
  % Build the complex potential, and calculate the absolute value.
  zc  = complex(rzc, izc);
  z   = abs(zc);

  % For x(radial)-direction we use at first a simple index.
  ind = reshape(floor((0:(ns*nx-1))/ns), ns, nx);

  for ii = 1:ns
    for jj = 1:nx
      % Compute the angle from the s-coordinate.
      % Note: This assumes circular geometry.
      theta(ii, jj) = fzero(@(x) (x + r(ii, jj)*sin(x) - 2*pi*s(ii, jj)), 2*pi*s(ii, jj));
    end
  end

  % Compute x- and y-(cartesian)coordinate.
  x = r.*cos(theta);
  y = r.*sin(theta);

  zeta = (q/pi).* atan(sqrt((1-r)./(1+r)).*tan(theta/2));
  %mesh(r, theta, zeta)
  z2 = z.*abs(sin(kzeta.*zeta));
  zc2 = zc.*exp(complex(0,1)*kzeta*zeta/rhostar);
  izc2 = imag(zc2);
  rzc2 = real(zc2);
  %mesh(x, y, z2);

  if(1 == strcmp('',write))
    % do nothing
  else
    % Store the data to file, such that it can be plotted with gnuplot.
    % The data is written in columns, where the meaning of the columns is
    % index r-coord s-coord x-coord y-coord z-coord mod_z-coord
    %
    % Also, before the r-coordinate (or the index) is increased, there is
    % a blank line (required by gnuplot).
    fileID = fopen(write, 'w');
    % Comment for indicating the columns.
    fprintf(fileID, '# index        eps             s               x               y               Re(zc)          Im(zc)          Re(zc2)         Im(zc2)\n');
    for ii = 1:nx
      fprintf(fileID, '%12.8e %12.8e %12.8e %12.8e %12.8e %12.8e %12.8e %12.8e %12.8e\n', [transpose(ind(:,ii)); transpose(r(:,ii)); transpose(s(:,ii)); transpose(x(:,ii)); transpose(y(:,ii)); transpose(rzc(:,ii)); transpose(izc(:,ii)); transpose(rzc2(:,ii)); transpose(izc2(:,ii));]);
      fprintf(fileID, '\n');
    end
    fclose(fileID);
  end

  % Compute thetazero.
  imean = 0.0;
  jmean = 0.0;
  zsum = 0.0;
  for ii = 1:ns
    for jj = 1:nx
      zsum = zsum + z(ii,jj);
      imean = imean + z(ii,jj)*ii;
      jmean = jmean + z(ii,jj)*jj;
    end
  end
  imean = imean/zsum;
  jmean = jmean/zsum;
  thetazero = fzero(@(x) (x + r(floor(imean), floor(jmean))*sin(x) - 2*pi*s(floor(imean), floor(jmean))), 2*pi*s(floor(imean), floor(jmean)));
endfunction
