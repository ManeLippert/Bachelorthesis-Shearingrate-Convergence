#!/usr/bin/octave --silent

%% This is a GNU Octave script.
%% It can be run with something like ./plot_entropy_balance.m



graphics_toolkit("gnuplot")
set(0, 'defaultfigurevisible', 'off');



%%%%%%%%%%%%%% USEFUL DEFINITIONS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function ret = get_color_palette(name)
  if (strcmp(name,"srg2"))
    ret=[hex2dec(reshape("66a61e",3,2))',
	 hex2dec(reshape("e3221c",3,2))',
	 hex2dec(reshape("d95f01",3,2))',
	 hex2dec(reshape("a6761d",3,2))',
	 hex2dec(reshape("1e78b4",3,2))',
	 hex2dec(reshape("a6cee3",3,2))',
	 hex2dec(reshape("199e77",3,2))',
	 hex2dec(reshape("e72b8a",3,2))',
	 hex2dec(reshape("666666",3,2))',
	 hex2dec(reshape("e6ab02",3,2))',
	 hex2dec(reshape("c2ca00",3,2))'
	 ]/255.0;
  else
    ret=[0.9, 0.0, 0.0,
	 0.9, 0.9, 0.0,
	 0.0, 0.9, 0.0,
	 0.9, 0.0, 0.9,
	 0.0, 0.9, 0.9];
  endif
endfunction

function ret = get_repeated_color_palette(name,ngraphs)
  %% this returns a palette with *twice* the first ngraphs colors from
  %% the palette with the given name
  c = get_color_palette(name);
  ret = [c(1:ngraphs,:),
	 c(1:ngraphs,:)];
endfunction

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
  endif
  
  commandline = ["grep -w --ignore-case '",pattern,"' ",infilename]
  [status, line_from_inputdat] = system (commandline)
  if(status != 0)
    printf (["problem with ",commandline],"\n");
    exit
  endif
  line_from_inputdat = tolower(line_from_inputdat);

  %% remove Fortran comments
  if(index(line_from_inputdat, "!") > 1) 
    line_from_inputdat = line_from_inputdat(1:index(line_from_inputdat, "!")-1);
  endif
  
  %% remove a trailing comma
  line_from_inputdat = strrep(line_from_inputdat, ',', '');

  ret = line_from_inputdat;
endfunction


function loglog_minusaware(x,y,palettename, lastgraphdashed)
  figure();
  
  if(nargin>=3)
    c = get_repeated_color_palette(palettename,columns(y));
  else
    c = get_repeated_color_palette("srg2",columns(y));
  endif

  set(gca(), 'colororder', c);

  if(nargin==4)
    loglog(x,max(y(:,1:end-1),0),'linestyle','-'
	   ,x,max(y(:,end),0),'linestyle',':'%,'color',c(1,end),
	   ,x,abs(min(y(:,1:end-1),0)),'linestyle','--'
	   ,x,abs(min(y(:,end),0)),'linestyle','--'%,'color',c(1,end)
	   )
  else
    loglog(x,max(y,0),'linestyle','-'
	   ,x,abs(min(y,0)),'linestyle','--'
	   )
  endif
  puttext("solid: > 0\ndashed: < 0","SouthWest")
endfunction

function semilogy_minusaware(x,y, palettename, lastgraphdashed)
  figure();

  if(nargin>=3)
    c = get_repeated_color_palette(palettename,columns(y));
  else
    c = get_repeated_color_palette("srg2",columns(y));
  endif
  set(gca(), 'colororder', c);
  
  if(nargin == 4)
    semilogy(x,max(y(:,1:end-1),0),'linestyle','-',
   	     x,max(y(:,end),0),'linestyle',':',
	     x,abs(min(y(:,1:end-1),0)),'linestyle','--',
	     x,abs(min(y(:,end),0)),'linestyle','--'
   	     )
  else
    semilogy(x,max(y,0),'linestyle','-',
   	     x,abs(min(y,0)),'linestyle','--'
   	     )
  endif
  puttext("solid: > 0\ndashed: < 0","SouthWest")
  
endfunction

function semilogyerr_minusaware(x,y,stddev, palettename, lastgraphdashed)
  figure();

  if(nargin>=4)
    c = get_repeated_color_palette(palettename,columns(y));
  else
    c = get_repeated_color_palette("srg2",columns(y));
  endif

  set(gca(), 'colororder', c);
  
  if(nargin == 5)
    semilogy(x,max(y(:,1:end-1),0),'linestyle','-',
   	     x,abs(min(y(:,end),0)),'linestyle','--',
	     x,max(y(:,1:end-1),0),'linestyle','-',
	     x,abs(min(y(:,end),0)),'linestyle','--'
   	     )
  else
    semilogy(x,max(y,0),'linestyle','-',
   	     x,abs(min(y,0)),'linestyle','--'
   	     )
  endif
  puttext("bars:\n1 standard deviation",'NorthEast');
  puttext("solid: > 0\ndashed: < 0","SouthWest")
  
endfunction

function print_all_formats(basename)
  SUFFIX{1} = ".eps";
  PRINTOPTS{1} = "-depsc";
  SUFFIX{2} = ".png";
  PRINTOPTS{2} = "-dpng";

  folder{1} = "eps/";
  folder{2} = "./";

  for i = 1:length(SUFFIX)
    if(~mkdir(folder{i}))
      error(["could not create folder '",folder{i},"'"])
    endif
    print([folder{i},basename,SUFFIX{i}], PRINTOPTS{i})
  endfor
  
endfunction

function ret = get_corresponding_index(time,skip_start_intervall)
  %% arguments:
  %%   time: a 2D array, the first column containing the time values.
  %%         For GKW, take time = dlmread("time.dat");
  %%   skip_start_intervall: a time intervall in time units.
  %% This function returns the (row-)index of the first sample after the time skip_start_intervall.
  
  
  %% skip something at the beginning
  if(issorted(time(:,1)))
    idx = lookup(time(:,1), skip_start_intervall);
    ret = idx + 1
    %% clip time series like this:
    %% clipped_time = time(idx+1:end,:);
    printf(["The first sample after ",num2str(skip_start_intervall)," [time units] is index ",int2str(idx),"\n"])
  else
    error("Something is not sorted!")
  endif

endfunction


function clipped_data = clip_beginning(data, time, skip_start_intervall)
  %% arguments:
  %%   data: a 2D array
  %%   time: a 2D array, the first column containing the time values.
  %%         For GKW, take time = dlmread("time.dat");
  %%   skip_start_intervall: a time intervall in time units.
  %% This function returns a part of the array data. It cuts away so many rows of data, that
  %% the signal begins at the first sample after skip_start_intervall.
  printf(["Skip ", num2str(skip_start_intervall)," [time units] at the beginning.\n"])
  
  clipped_data = data(get_corresponding_index(time,skip_start_intervall):end,:);
endfunction

%%%%%%%%%%% SETTINGS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
palettename="srg2"
set(0, 'defaultaxescolororder', get_color_palette(palettename));
set(0, 'defaultlinelinewidth', 6);


%%%%%%%%%%%%%%%%%% LOAD DATA %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
global time ;
global krho ;

global dt_entr ;
global dt_entr_field ;
global entr ;
global entr_field ;
global sources ;
global num_source01 ;
global num_source02 ;
global num_source03 ;
global dealiasing ;
global coll ;
global temp_src ;
global outflow ;

global ene_f ;
global ene_e ;
global ene_e2 ;


global total ;
global total_norm ;
global quotient ;
global mode_group_labels;

%% sum over all modes
global SUM_ALL_MODE=0;

%% sum over all modes except the zeromode
global SUM_ALL_EXCEPT_ZERO_MODE=3;

%% use all modes, load the spectrum
global ALL_INDIVIDUALLY_MODE=4;


function load_time()
  global time;
  %% the first column of time.dat contains the time
  time = dlmread("time.dat");
endfunction

function ret = remove_parseval_factor2(data)
  
  if(columns(data) > 1)
    %% the first column is the zero mode and is left as it is.
    ret = data;
    %% the other fourierkoeffs have received a factor 2 before (for simple integration via Parseval's theorem)
    %% in the diagnostic.
    %% this function removes this factor 2.
    ret(:,2:end) = data(:,2:end) * 0.5;
  else
    ret = data * 0.5;
  endif
endfunction

function load_energetics_diagn(MODE, igroup_or_imod, groups)
  %% load the output files of the energetics diagnostic of GKW

  global time
  
  global dt_entr
  global dt_entr_field
  global entr
  global entr_field
  global sources
  global num_source01
  global num_source02
  global num_source03
  global dealiasing
  global coll
  global temp_src
  global outflow

  global ene_f ;
  global ene_e ;
  global ene_e2 ;

  global total
  global total_norm
  global quotient
  global quotient_with_outflow

  global SUM_ALL_MODE
  global mode_group_labels;
  global SUM_ALL_EXCEPT_ZERO_MODE
  global ALL_INDIVIDUALLY_MODE

  if(nargin() < 1)
    error("not enough input arguments");
  endif

  _dt_entr = dlmread("dt_entr.kyspec");
  _dt_entr_field = dlmread("dt_entr_field.kyspec");
  _entr = dlmread("entr.kyspec");
  _entr_field = dlmread("entr_field.kyspec");
  _src01 = dlmread("entr_source01.kyspec");
  _src02 = dlmread("entr_source02.kyspec");
  _src03 = dlmread("entr_source03.kyspec");
  _src04 = dlmread("entr_source04.kyspec");
  _src05 = dlmread("entr_source05.kyspec");
  _src06 = dlmread("entr_source06.kyspec");  
  _sources = _src01 + _src02 + _src03 + _src04 + _src05 + _src06;
  _num_source01 = dlmread("entr_num_dis.kyspec");
  _num_source02 = dlmread("entr_num_vp.kyspec");
  _num_source03 = dlmread("entr_num_perp.kyspec");

  _dealiasing = _dt_entr*0;
  _dealiasing(:,columns(_dealiasing)) = dlmread("entr_src_dealiasing.dat");
  _coll = dlmread("entr_coll.kyspec");
  _temp_src = dlmread("entr_temp_src.kyspec");
  _outflow = dlmread("entr_outflow.kyspec");

  %% the outflow can be *very* small. In a logarithmic plot, this will lead to
  %% unpractical ranges. Therefore suppress (=do not plot) very small values of
  %% the outflow dissipation.
  tinyval=1e-10;
  _outflow(_outflow<tinyval) = 0;

  _ene_f = dlmread("ene_f.kyspec");
  _ene_e = dlmread("ene_e.kyspec");
  _ene_e2 = dlmread("ene_e2.kyspec");

  %% the entropy contribution from dealiasing is just a scalar				
  if (MODE == SUM_ALL_MODE)
    %% Sum over all modes (i.e. columns).
    %% To do the sum properly according to Parsevals theorem and the FFTW data structure
    %% the following lines assume that the data contains a factor 2 on all but the zero-mode column
    dt_entr = sum(_dt_entr,2);
    dt_entr_field = sum(_dt_entr_field,2);
    entr = sum(_entr,2);
    entr_field = sum(_entr_field,2);
    sources = sum(_sources,2);
    num_source01 = sum(_num_source01,2);
    num_source02 = sum(_num_source02,2);
    num_source03 = sum(_num_source03,2);
    coll = sum(_coll,2);
    temp_src = sum(_temp_src,2);
    outflow = sum(_outflow,2);

    dealiasing = sum(_dealiasing,2);

    ene_f = sum(_ene_f,2);
    ene_e = sum(_ene_e,2);
    ene_e2 = sum(_ene_e,2);

    mode_group_labels = {"sum of modes"};

  elseif (MODE == SUM_ALL_EXCEPT_ZERO_MODE)
    %% Sum over all modes (i.e. columns), except the zeromode (i.e the first column)
    dt_entr = sum(_dt_entr(:,2:end),2);
    dt_entr_field = sum(_dt_entr_field(:,2:end),2);
    entr = sum(_entr(:,2:end),2);
    entr_field = sum(_entr_field(:,2:end),2);
    sources = sum(_sources(:,2:end),2);
    num_source01 = sum(_num_source01(:,2:end),2);
    num_source02 = sum(_num_source02(:,2:end),2);
    num_source03 = sum(_num_source03(:,2:end),2);
    coll = sum(_coll(:,2:end),2);
    temp_src = sum(_temp_src(:,2:end),2);
    outflow = sum(_outflow(:,2:end),2);

    dealiasing = sum(_dealiasing(:,2:end),2);

    ene_f = sum(_ene_f(:,2:end),2);
    ene_e = sum(_ene_e(:,2:end),2);
    ene_e2 = sum(_ene_e(:,2:end),2);

    mode_group_labels = {"sum of modes"};
    
    elseif(MODE == ALL_INDIVIDUALLY_MODE)

      %% take all modes
      for i=1:columns(_dt_entr)
	mode_group_labels{i}=["mode ",int2str(i)];
      endfor
      
      dt_entr = _dt_entr;
      dt_entr_field = _dt_entr_field;
      entr = _entr;
      entr_field = _entr_field;
      sources = _sources;
      num_source01 = _num_source01;
      num_source02 = _num_source02;
      num_source03 = _num_source03;
      coll = _coll;
      temp_src = _temp_src;
      outflow = _outflow;
      dealiasing = _dealiasing;

      ene_f = _ene_f;
      ene_e = _ene_e;
      ene_e2 = _ene_e2;
    
  else
    disp("error loading energetics data.")
    exit 1
  endif
    
  total=dt_entr+dt_entr_field+sources+num_source01+num_source02+num_source03+coll+temp_src+dealiasing+outflow;
  
  total_norm=total./(entr+entr_field);

  quotient = (dt_entr+dt_entr_field)./(sources + num_source01 + num_source02 + num_source03 + coll + temp_src + dealiasing + outflow);
endfunction

load_time()
load_energetics_diagn(SUM_ALL_MODE)


%%%%%%%%%%%%%%% CLIP TRANSIENT PHASE %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%% PLOT INTEGRALS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%% EVERYTHING IN ONE PLOT, MANY GRAPHS %%%%%%

y = [dt_entr, dt_entr_field, sources, num_source01, num_source02, num_source03, coll, dealiasing, outflow, total];
semilogy_minusaware(time(:,1)
		    ,y
		    ,"srg2"
		    ,"lastgraphdashed"
		    );

xlabel("time");
legend("dSdt",
       "dWdt",
       "fluxes*gradients",
       "parall. diss.",
       "vel. parall. diss.",
       "perp. diss.",
       "collisional dissipation",
       %"temperature source",
       "dealiasing dissipation",
       "outflow dissipation",
       "total",
       "location","north");
print_all_formats("balancelog")

%%%%%% EVERYTHING IN ONE PLOT,  TWO GRAPHS %%%%%%%%
y = [(dt_entr + dt_entr_field), (sources + num_source01 + num_source02 + num_source03 + coll + temp_src + dealiasing + outflow), total ];

%% do not plot the time=0 sample in a loglog plot.
loglog_minusaware(time(2:end,1),
		  y(2:end,:),
		  palettename,
		  "lastgraphdashed"
	 );

xlabel("time");
legend("dSdt+dWdt",
       "fluxInGrads+num.dissip.+coll.+dealias.+outflow",
       "total",
       "location","SouthEast");

print_all_formats("balancelog2")



%%%%%%%%%%%%% LOAD DATA, NOW SPECTRALLY %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

eval(from_inputdat("nmod"))
if(nmod == 1)
  disp("run has only one mode. done")
  exit 0
endif

function load_wavenumbers()
  
  global krho;

  %% (q/2πψ)kζ ρref = k⊥ (kψ = 0, s = 0)ρref
  krho = dlmread("krho");

endfunction

load_energetics_diagn(ALL_INDIVIDUALLY_MODE)
load_wavenumbers()

SKIP_START_INTERVALL = 0.4 * time(rows(time),1);

entr = clip_beginning(entr, time, SKIP_START_INTERVALL);
entr = remove_parseval_factor2(entr);
entr_mean = mean(entr,1)
entr_std = std(entr)

entr_field = clip_beginning(entr_field, time, SKIP_START_INTERVALL);
entr_field = remove_parseval_factor2(entr_field);
entr_field_mean = mean(entr_field,1)
entr_field_std = std(entr_field)

total = clip_beginning(total, time, SKIP_START_INTERVALL);
total = remove_parseval_factor2(total);
total_mean = mean(total,1)
total_std = std(total)

dt_entr = clip_beginning(dt_entr, time, SKIP_START_INTERVALL);
dt_entr = remove_parseval_factor2(dt_entr);
dt_entr_mean = mean(dt_entr,1)
dt_entr_std = std(dt_entr)

dt_entr_field = clip_beginning(dt_entr_field, time, SKIP_START_INTERVALL);
dt_entr_field = remove_parseval_factor2(dt_entr_field);
dt_entr_field_mean = mean(dt_entr_field,1)
dt_entr_field_std = std(dt_entr_field)

sources = clip_beginning(sources, time, SKIP_START_INTERVALL);
sources = remove_parseval_factor2(sources);
sources_mean = -mean(sources,1)
sources_std = std(sources)

num_sources = clip_beginning(num_source01 + num_source02 + num_source03, time, SKIP_START_INTERVALL);
num_sources = remove_parseval_factor2(num_sources);
num_sources_mean = -mean(num_sources,1)
num_sources_std = std(num_sources)

num_source01 = clip_beginning(num_source01, time, SKIP_START_INTERVALL);
num_source01 = remove_parseval_factor2(num_source01);
num_source01_mean = -mean(num_source01,1)
num_source01_std = std(num_source01)

num_source02 = clip_beginning(num_source02, time, SKIP_START_INTERVALL);
num_source02 = remove_parseval_factor2(num_source02);
num_source02_mean = -mean(num_source02,1)
num_source02_std = std(num_source02)

num_source03 = clip_beginning(num_source03, time, SKIP_START_INTERVALL);
num_source03 = remove_parseval_factor2(num_source03);
num_source03_mean = -mean(num_source03,1)
num_source03_std = std(num_source03)

coll = clip_beginning(coll, time, SKIP_START_INTERVALL);
coll = remove_parseval_factor2(coll);
coll_mean = mean(coll,1)
coll_std = std(coll)

time = clip_beginning(time, time, SKIP_START_INTERVALL);

abscissa = krho(:,1)
abscissa_label = "bi-normal wavenumber";

%%eval(from_inputdat("krhomax"))
%% You have to edit this line probably, according to your krhomax:
global xlimits = [(-0.5*(abscissa(2) - abscissa(1))) 2.8+(0.5*(abscissa(2) - abscissa(1)))];

%%%%%%%%%%%%%%%%%% PLOT SPECTRUM %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function spectral_bar_plot(x,y,stddev)
  global xlimits;
  global mode_group_labels;

  figure();

  %% workaround: clip the bars to the ylimits (otherwise they may not be drawn...)
  # std_uy=min(stddev,ylimits(2)-y)
  # std_ly=min(stddev,y-ylimits(1))

  # err_uy=min(err,ylimits(2)-y)
  # err_ly=min(err,y-ylimits(1))
 
  groups = columns(y);
  hold on
  barwidth=1;
  bar(x, diag(y), barwidth, "stacked")
  ylimits = get(gca(), "ylim");
  legend(mode_group_labels, 'Location', 'SouthEast');

  errorbar(x,y,stddev,'~.k')
  puttext("bars:\n1 standard deviation",'NorthEast');
  hold off
  set(gca(), "xlim", xlimits);
  set(gca(), "ylim", ylimits);
endfunction

spectral_bar_plot(abscissa, entr_mean + entr_field_mean, entr_std + entr_field_std)
xlabel(abscissa_label)
ylabel("entropy (S+W) time-averaged")
print_all_formats("spec_entr_mean")

%%%%%%%% and now relative %%%%%%%%%%%%%%%m

spectral_bar_plot(abscissa,
		  mean((dt_entr+dt_entr_field)./(entr + entr_field)),
		  std((dt_entr+dt_entr_field)./(entr + entr_field))
		  );
xlabel(abscissa_label)
ylabel("d/dt (S + W)  / entropy in mode, time averaged")
print_all_formats("spec_rel_dt_entr_mean")

spectral_bar_plot(abscissa,
		  mean(total./(entr + entr_field)),
		  std(total./(entr + entr_field))
		  );
xlabel(abscissa_label)
ylabel("equation (left side - right side) time-averaged  / entropy in mode")
print_all_formats("spec_rel_total_mean")

spectral_bar_plot(abscissa,
		  mean(sources./(entr + entr_field)),
		  std(sources./(entr + entr_field))
		  );
xlabel(abscissa_label)
ylabel("fluxes * gradients  / entropy in mode, time averaged")
print_all_formats("spec_rel_sources_mean")

spectral_bar_plot(abscissa,
		  mean(num_sources ./(entr + entr_field)),
		  std(num_sources ./(entr + entr_field))
		  );
xlabel(abscissa_label)
ylabel("total numerical dissipation  / entropy in mode, time-averaged")
print_all_formats("spec_rel_num_sources_mean")

spectral_bar_plot(abscissa,
		  mean(coll ./(entr + entr_field)),
		  std(coll ./(entr + entr_field))
		  );
xlabel(abscissa_label)
ylabel("collisional entropy contribution  / entropy in mode, time-averaged")
print_all_formats("spec_rel_coll_mean")


%%%%%%% SPECTRAL LOGLOG LINE PLOTS %%%%%%%%%%%%%%%%%%%%%%%

function spectral_loglog_line_plot(x,y,stddev)
  loglog(x,max(y,0),'marker','+',%'color',c(1,:),
	 x,abs(min(y,0)),'marker','*'%,'color',c(2,:)
	 )
  legend({"positive sign","negative sign"}, 'Location', 'SouthEast');
  %%set(gca(), "xlim", [(-0.5*(x(2) - x(1))) 2.8]);
endfunction


spectral_loglog_line_plot(abscissa,entr_mean + entr_field_mean,std(entr+entr_field))
xlabel(abscissa_label)
ylabel("entropy S+W time-averaged")
print_all_formats("loglog_spec_entr_mean")

%%%%%%%%%%% SEMILOGARITHMIC PLOTS %%%%%%%%%%%%%%%%%%%%%%%%%

%% in those one can also see the zero-mode...

semilogyerr_minusaware(abscissa,entr_mean,entr_std)
xlabel(abscissa_label)
ylabel("entropy S time-averaged")
print_all_formats("semilog_spec_entr_mean")
