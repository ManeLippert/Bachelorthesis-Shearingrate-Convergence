%This script makes a matlab function that plots a set of time traces from
%GKW output onto one plot.  Requires gkwpath to be configured.
%
%Usage: time_traces('proj','filename','data',column,optplot,optline)
%
%proj is a string with your project name
%Filename to select which runs of project to plot (accepts wildcards).
%Data file to plot (default is fluxes) - as goes into gkwpath input
%Column to plot (default is 2)
%optplot selects the time vs timestep plot  (optional)
%optline is the plot line type e.g. '-')   (optional)
%
%Example: time_traces('my_proj','runs*','fluxes',2)

function[]=time_traces(proj,filename,data,col,optplot,line)
 
if ~exist('proj','var')
	proj='default';
    disp('You must provide the project name')
    %No default value
end
if ~exist('filename','var')
    %list all
    filename='';
end
if ~exist('data','var')
    data='fluxes'
end
if ~exist('optplot','var')
    optplot=0;
end
if ~exist('line','var')
    line='-';
end


count=0;            %Counter for output files
no_files=0;         %Counter for input files

% if ~exist('names')
%  %files=dir([gkwpath('root',proj) 'good/' filename]);
  files=dir([gkwpath('time',proj) filename]);
  total = size(files)
% else
%  files=names 
%  [files(1:size(files')).isdir]=deal(0);
%  total = size(files')
% end

if optplot ==0; h=figure; end
hold all

for i = 1:total(1)
   
   if(files(i).isdir==0)
        no_files=no_files+1;
        
        %input=read_gkwinput(files(i).name,proj);
        %gkwpath(data,proj)
        disp([gkwpath(data,proj) files(i).name]);
        if (exist([gkwpath(data,proj) files(i).name],'file')==2)
            files(i).name;
            count = count+1;
            time=load([gkwpath('time',proj) files(i).name]);
            datafile=load([gkwpath(data,proj) files(i).name]); 
                        
            len=min(length(time),length(datafile));
            
            %These lines renormalise to gyrobohm if desired
            %time(:,1)=time(:,1)/(sqrt(2)*3);
            %data(:,col)=data(:,col)/9.545;
            
            hold all;
            
            %Add the plot
            if (length(time) > 1)
                %if (~isnan(datafile(end,col)) && datafile(end,col) < 1e3) 
                  if optplot==0
                      figure(h) 
                      subplot(1,2,1)
                  end
                  %size(time)
                  %size(datafile)
                  plot(time(1:len,1),datafile(1:len,col),line,'DisplayName',files(i).name);
                  hold all
                %end
            end
            
            if optplot==0
                subplot(1,2,2)
                plot(time(:,1),[1:1:length(time)],line);
                xlabel('time');
                ylabel('timestep');
                hold all
            end
            
        else
            disp(['Skipped missing', gkwpath(data,proj) files(i).name])
        end
        %cd ..
   end
end

if (no_files==0)
    disp('No input files found')
    dir([gkwpath('input',proj) filename])
end    
if (count==0)
    disp('No output files found')
    dir([gkwpath(data,proj) filename])
else    
disp(sprintf('Loaded %i output files for %i input files', count, no_files))
end

if optplot==0; subplot(1,2,1); end
xlabel('t (v_{thi} / R)');
ylabel([data 'column' col]);

end
