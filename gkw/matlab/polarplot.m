function [ x,y,rzc2,izc2 ] = polarplot( parallel,prof_back,n,nx,ns,rhostar,field,interpnum)
% polarplot 
% This function will extract from the data of a run the grids as well as
% the perturbed potential.
% The x,y coordinates are calculated from these, also the potential for
% constant phi and plot poloidal slices of the real and imag parts
%%
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
%   interpnum  The interger which determines how many times the parallel
%              direction is interpolated.....8-16 is recommended
%
%
% Output parameters:
%   x:         Grid with the x-coordinate.
%   y:         Grid with the y-coordinate.
%   rzc2:      Array with the real part of the perturbed rescaled field.
%   izc2:      Array with the imgainary part of the perturbed rescaled field.
% ----------------------------------------------------------------------
    kzeta = 2*pi*rhostar*n;
    r = prof_back(:,1);
    r = repmat(r.', ns, 1);% r-coordinate does not depend on s, so simply repeat the vector ns times.
    q = prof_back(:,2);
    q = repmat(q.',ns, 1);
    % Get the s-coordinate from the data. Also read the data, in this case it is the square of the value of the perturbed potential.
    s = reshape(parallel(1:ns*nx, 1), ns, nx);
    rzc = reshape(parallel(1:ns*nx, 2*field + 0), ns, nx);
    izc = reshape(parallel(1:ns*nx, 2*field + 1), ns, nx);

    ds = parallel(2,1)-parallel(1,1);
    r(end+1,:)=r(1,:);
    q(end+1,:)=q(1,:);
    s(end+1,:)=s(end,:) + ds;
    rzc(end+1,:)=rzc(1,:);
    izc(end+1,:)=izc(1,:);
    ns=ns+1;%Since we have interpolated a bit
    zc  = complex(rzc, izc);
    zc(end,:) = zc(end,:);%.*exp(-1i*kzeta.*q(1,:));
    
    for ii = 1:nx
        r_i(:,ii) = interp(r(:,ii),interpnum);
        q_i(:,ii) = interp(q(:,ii),interpnum);
        s_i(:,ii) = interp(s(:,ii),interpnum);
        zc_i(:,ii) = interp(zc(:,ii),interpnum);
    end

    % Build the complex potential, and calculate the absolute value.
    
    z   = abs(zc_i);
    % For x(radial)-direction we use at first a simple index.
    %ind = reshape(floor((0:(ns*nx-1))/ns), ns, nx);
    for ii = 1:ns*interpnum
      for jj = 1:nx
        % Compute the angle from the s-coordinate.
        % Note: This assumes circular geometry.
        theta(ii, jj) = fzero(@(x) (x + r_i(ii, jj)*sin(x) - 2*pi*s_i(ii, jj)), 2*pi*s_i(ii, jj));
      end
    end
    % Compute x- and y-(cartesian)coordinate.
    x = r_i.*cos(theta);
    y = r_i.*sin(theta);
    zeta = (q_i/pi).* atan(sqrt((1-r_i)./(1+r_i)).*tan(theta/2));

    %mesh(r, theta, zeta)
    z2 = z.*abs(sin(kzeta.*zeta));
    zc2 = zc_i.*exp(complex(0,1)*kzeta*zeta/rhostar);
    izc2 = imag(zc2);
    rzc2 = real(zc2);
    figure;pcolor(x,y,rzc2);shading flat
    figure;pcolor(x,y,izc2);shading flat
end

