% load the metrics from various files and plot all the terms for comparison

% uncomment the line corresponding to the eps value you want to plot (and comment the other one)
flnm={'GeomComp_sa_1','GeomComp_circ_1','GeomComp_ch_1'}; %eps=0.0507
%flnm={'GeomComp_sa_2','GeomComp_circ_2','GeomComp_ch_2'}; %eps=0.1498


c={'b','m','r'}; % s-alpha, circ, chease

for ii = 1:length(flnm)

  % read data
  flpth=pwd;
  flpth=[flpth '/'];
  if unix(['test -e ' flpth flnm{ii}])==0
	fid = fopen([flpth flnm{ii}], 'r');
  else	
	error(['The file ' flpth flnm{ii} ' does not exist' ])
  end
  frewind(fid);
  sss='';
  % scalars
  for jj=1:13
	sss=deblank(fgets(fid));
	eval(['G.' lower(sss) '=fscanf(fid,''%f'',1);'])
	fgets(fid);
  end
  % 1_D quantities
  for jj=1:26
	sss=deblank(fgets(fid));
	eval(['G.' lower(sss) '=fscanf(fid,''%f'',G.ns);'])
	fgets(fid);
  end
  fclose(fid);


  % plot data
  F = fieldnames(G);
  for jj=1:length(F)
    if eval(['length(G.' F{jj} ')'])==G.ns
       figure(jj)
       hold on
       eval(['plot(G.s_grid,G.' F{jj} ',c{ii})']);
       title(F{jj})
    end
  end

  clear G

end


