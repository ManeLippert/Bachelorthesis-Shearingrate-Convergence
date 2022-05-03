#!/usr/bin/env octave

%% This is a GNU Octave script.

% this is important:
1;

function init()
  graphics_toolkit("gnuplot")
  set(0, 'defaultfigurevisible', 'off');
  ## c = get_color_palette("srg2");
  ## assignin('caller','c',c);
  ## set(0, 'defaultaxescolororder', c);
  set(0, 'defaultlinelinewidth', 6);
  set(0, 'defaultaxesygrid','on');
  set(0, 'defaultaxesfontweight', 'bold');
  set(0, 'defaultaxeslinewidth', 3);
  set(0, 'defaultaxesfontsize',10);
  % setting the font is a workaround to prevent an error that may
  % happen when printing plot to PNG files:
  set(0, 'defaultaxesfontname', 'Helvetica')
  set(0, 'defaulttextinterpreter','Tex')
  set(0,'defaultfigurepapertype', '<custom>')
  set(0,'defaultfigurepaperunits','points');
  set(0,'defaultfigurepaperposition',[0 0 660 450])
%get (0, "factory") %returns a list of factory defaults. get (0,
%"default") %returns a list of user-defined default values
end



function ret = from_inputdat(pattern,filename)
  %% Use this function like this:
  %%    eval(from_inputdat("vpmax"))
  %%    eval(from_inputdat("mumax"))
  %% Note that this function will grep from input.out, not input.dat!
  %% So logical values are either 'T' or 'F' in this file.
  %% A useful trick is to do it like in this example:
  %%    f = false;
  %%    t = true;
  %%    eval(from_inputdat("spectral_radius"))
  %% Lists of floats can be evaluated like this
  %%  eval([strrep(from_inputdat("parallel_output_timestamps", "input.dat"),'=','=[') , ']'])
  %%
  %% If the optional parameter filename is given,
  %% then given file is read instead of input.out.

  infilename="input.out";
  if(nargin == 2)
      infilename=filename;
  end
  
  commandline = ["grep -w --ignore-case '",pattern,"' ",infilename];
  [status, line_from_inputdat] = system (commandline);
  if(status != 0)
    printf (["problem with ",commandline],"\n");
    exit
  end
  line_from_inputdat = tolower(line_from_inputdat);

  %% remove Fortran comments
  if(index(line_from_inputdat, "!") > 1) 
    line_from_inputdat = line_from_inputdat(1:index(line_from_inputdat, "!")-1);
  end
  
  %% remove a trailing comma
  line_from_inputdat = strrep(line_from_inputdat, ',', '');

  ret = line_from_inputdat;
end


function print_all_formats(basename,extraoptions)
  %% from the docs: -Sxsize,ysize Plot size in pixels for EMF, GIF,
  %% JPEG, PBM, PNG and SVG. For PS, EPS, PDF, and other vector
  %% formats the plot size is in points.
         
  %% suffix, foldername, options..	 
  outformat{1} = {"png","./","-dpng","-color"};
  %outformat{end+1} = {"eps","eps/","-depsc","-color","-tight"};
  %outformat{end+1} = {"eps","","-depslatex","-color","-tight","-r300"}; % I failed to obtain a plot with the right size with epslatex
  %outformat{end+1} = {"tex","tex/","-dtex","-color","-tight"};
  outformat{end+1} = {"tex","tikz/","-dtikz","-color","-tight"};
  %outformat{end+1} = {"eps","epsmono/","-depsc","-mono","-tight"};
  %outformat{end+1} = {"epsc","epsc2/","-depsc2","-tight"};
  %outformat{end+1} = {"pdf","pdf/","-dpdf"}
  %outformat{end+1} = {"pdf","pdfwrite/","-dpdfwrite"}

  for i = 1:length(outformat)
      folder = outformat{i}{2} 
    if(length(folder) > 0 && ~mkdir(folder))
      error(["could not create folder '",folder,"/'"])
    endif
    
    printcommand =["print ", outformat{i}{2},basename,'.',outformat{i}{1}," "];
    for j = 3:length(outformat{i})
        printcommand = [printcommand," ", outformat{i}{j}];
    end
    %(2:length(outformat{i})-1) = outformat{i}(3:end)
    
    if(nargin() == 2)
      for j = 1:length(extraoptions)
        printcommand = [printcommand, " ", extraoptions{j}];
      end
    end
    printcommand
    eval(printcommand)

    if(strcmp(outformat{i}{3},"-depslatex"))
      %% this is a workaround. I would like to avoid that the text
      %% file contains 'epslatex/foo.eps'
      mkdir("epslatex")
      rename([basename,'.eps'],['epslatex/',basename,'.eps'])
      rename([basename,'.tex'],['epslatex/',basename,'.tex'])
    end
  endfor
  
endfunction



function hcb = give_imagesc_nancolor(a,colormap_,nancolor)	 
  %% IMAGESC with NaNs assigning a specific color to NaNs
  


  %% find minimum and maximum
  amin=min(a(:));
  amax=max(a(:));
  maxbound = max(abs(amin), abs(amax));
  %% size of colormap
  n = size(colormap_,1);
  %% color step
  %dmap=(amax-amin)/n;
  dmap=(2 * maxbound)/n;

  %% add nan color to colormap
  colormap([nancolor; colormap_]);
  
  %% changing color limits
  %caxis([amin-dmap amax]);
  caxis([-maxbound-dmap maxbound]);
  %% place a colorbar
  hcb = colorbar();
  %%  change Y limit for colorbar to avoid showing NaN color
  %clim(hcb, [amin amax]) %function undefined in octave?

end

function puttext(string, location)
  %% this function allows to write a string on a plot.
  %% instead of coordinates like the usual text(...) functions
  %% it takes a location like in
  %% legend("the foobar function", "location","NorthEast");

  switch(location)
    case 'North'
      t = text(0.5, 0.95, string,
	       'horizontalalignment','center',
	       'units','normalized'
	       );
    case 'NorthEast'
      t = text(0.95, 0.95, string,
	       'horizontalalignment','right',
	       'units','normalized'
	       );
    case 'East'
      t = text(0.95, 0.5, string,
	       'horizontalalignment','right',
	       'units','normalized'
	       );
    case 'SouthEast'
      t = text(0.95, 0.05, string,
	       'horizontalalignment','right',
	       'units','normalized'
	       );
    case 'South'
      t = text(0.5, 0.05, string,
	       'horizontalalignment','center',
	       'units','normalized'
	       );
    case 'SouthWest'
      t = text(0.05, 0.05, string,
	       'horizontalalignment','left',
	       'units','normalized'
	       );
    case 'West'
      t = text(0.05, 0.5, string,
	       'horizontalalignment','left',
	       'units','normalized'
	       );
    case 'NorthWest'
      t = text(0.05, 0.95, string,
	       'horizontalalignment','left',
	       'units','normalized'
	       );
    otherwise
      %South
      t = text(0.5, 0.05, string,
	       'horizontalalignment','center',
	       'units','normalized'
	       );
  endswitch

endfunction

function ret = dlmread_datbin(luname)
  filename = [luname,'.dat'];
  [info, err, msg] = stat(filename);
  if(err ~= -1)
    ret = dlmread(filename);
    return;
  end

  filename = [luname,'.bin'];
  [info, err, msg] = stat(filename);
  if(err ~= -1)
    ret = dlmread(filename);
    return;
  end

  % this option should not exist with newer GKW versions:
  filename = luname
  ret = dlmread(filename);
end
