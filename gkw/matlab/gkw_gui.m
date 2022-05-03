function varargout = gkw_gui(varargin)
% GKW_GUI M-file for gkw_gui.fig
%      GKW_GUI, by itself, creates a new GKW_GUI or raises the existing
%      singleton*.
%
%      H = GKW_GUI returns the handle to a new GKW_GUI or the handle to
%      the existing singleton*.
%
%      GKW_GUI('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in GKW_GUI.M with the given input arguments.
%
%      GKW_GUI('Property','Value',...) creates a new GKW_GUI or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before gkw_gui_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to gkw_gui_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help gkw_gui

% Last Modified by GUIDE v2.5 13-Feb-2010 22:26:22

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @gkw_gui_OpeningFcn, ...
                   'gui_OutputFcn',  @gkw_gui_OutputFcn, ...
                   'gui_LayoutFcn',  [] , ...
                   'gui_Callback',   []);
if nargin && ischar(varargin{1})
    gui_State.gui_Callback = str2func(varargin{1});
end

if nargout
    [varargout{1:nargout}] = gui_mainfcn(gui_State, varargin{:});
else
    gui_mainfcn(gui_State, varargin{:});
end
% End initialization code - DO NOT EDIT



% --- Executes just before gkw_gui is made visible.
function gkw_gui_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to gkw_gui (see VARARGIN)

% Global variables 
global gkw_gui_directory_name
global gkw_gui_filename
global gkw_gui_flux_select
global gkw_gui_par_select
global gkw_gui_current_plot

% Choose default command line output for gkw_gui
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

%Initialize
gkw_gui_directory_name = pwd;
gkw_gui_filename       = '';

%Check if an input.dat exists in the parrent directory 
name = [gkw_gui_directory_name '/input.dat'];
alf = dir(name);
if (isempty(alf)) 
  set(handles.Filename,'String','*** Select file ***');
else 
  read_all_data; 
  set(handles.Filename,'String','*** Single dir ***');
end 
set(handles.Directory_files,'String',gkw_gui_directory_name);

% Set all the selection switches to on 
set(handles.Select_species_1,'Value',1);
set(handles.Select_species_2,'Value',1);
set(handles.Select_species_3,'Value',1);
set(handles.Select_species_all,'Value',1);
set(handles.Particle_flux_selection,'Value',1);
set(handles.Heat_flux_selection,'Value',1);
set(handles.Mom_flux_selection,'Value',1);
set(handles.Select_dens,'Value',0);
set(handles.Select_phi,'Value',1);
set(handles.Select_Apar,'Value',0);
set(handles.Select_Bpar,'Value',0);
set(handles.Select_Tpar,'Value',0);
set(handles.Select_Tperp,'Value',0);
set(handles.Select_vpar,'Value',0);

% Set all flux select values to 1. 
for i =1 : 3 
  gkw_gui_flux_select(i) = 1; 
end 
for i = 1 : 7 
  gkw_gui_par_select(i) = 0; 
end
gkw_gui_par_select(1) = 1; 

% at pressent there is no plot selected 
gkw_gui_current_plot = 'none';

%For manipulation of the figure
set(hObject,'toolbar','figure');

% UIWAIT makes gkw_gui wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = gkw_gui_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


% --- Executes on button press in Growth_rate_freq.
function Growth_rate_freq_Callback(hObject, eventdata, handles)
% hObject    handle to Growth_rate_freq (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

%used global variables 
global gkw_gui_current_plot 

gkw_gui_current_plot = 'time';
axes(handles.Figure1);
redo_the_plot;


function Directory_files_Callback(hObject, eventdata, handles)
% hObject    handle to Directory_files (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of Directory_files as text
%        str2double(get(hObject,'String')) returns contents of Directory_files as a double

%used global variables 
global gkw_gui_directory_name

%Get the name 
gkw_gui_directory_name = get(hObject,'String');

% Check if there is an input.dat 
name = [gkw_gui_directory_name '/input.dat'];
alf = dir(name);
if (isempty(alf)) 
  set(handles.Filename,'String','*** Select a file ***');
else 
  read_all_data(handles); 
  set(handles.Filename,'String','*** Single dir ***');
end 

% --- Executes during object creation, after setting all properties.
function Directory_files_CreateFcn(hObject, eventdata, handles)
% hObject    handle to Directory_files (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function Filename_Callback(hObject, eventdata, handles)
% hObject    handle to Filename (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of Filename as text
%        str2double(get(hObject,'String')) returns contents of Filename as a double

%used globales
global gkw_gui_filename
gkw_gui_filename = get(hObject,'String'); 
read_all_data(handles); 


% --- Executes during object creation, after setting all properties.
function Filename_CreateFcn(hObject, eventdata, handles)
% hObject    handle to Filename (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in File_select.
function File_select_Callback(hObject, eventdata, handles)
% hObject    handle to File_select (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

global gkw_gui_directory_name
global gkw_gui_filename

name = [gkw_gui_directory_name '/input'];
alf = dir(name); 
if (isempty(alf)) 
  disp('Could not find files in the input directory. I was looking in:');
  disp(name);
else 
  current_directory = pwd; 
  cd(name); 
  [file path filt] = uigetfile('*','Choose a GKW input file'); 
  gkw_gui_filename = file; 
  set(handles.Filename,'String',gkw_gui_filename);
  cd(current_directory);
  read_all_data(handles);
end 


% --- Executes on button press in Fluxes.
function Fluxes_Callback(hObject, eventdata, handles)
% hObject    handle to Fluxes (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

%used global variables 
global gkw_gui_current_plot 

%Set the plotting to fluxes
gkw_gui_current_plot = 'fluxes';

%Redo the plotting 
axis(handles.Figure1); 
redo_the_plot; 


% --------------------------------------------------------------------
function Fluxes_select_ButtonDownFcn(hObject, eventdata, handles)
% hObject    handle to Fluxes_select (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% --- Executes on button press in Particle_flux_selection.
function Particle_flux_selection_Callback(hObject, eventdata, handles)
% hObject    handle to Particle_flux_selection (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of Particle_flux_selection
global gkw_gui_flux_select

if (get(hObject,'Value') == get(hObject,'Max'))
   % Checkbox is checked-take approriate action
   gkw_gui_flux_select(1) = 1; 
else
   % Checkbox is not checked-take approriate action
   gkw_gui_flux_select(1) = 0;
end

axes(handles.Figure1); 
redo_the_plot; 


% --- Executes on button press in Heat_flux_selection.
function Heat_flux_selection_Callback(hObject, eventdata, handles)
% hObject    handle to Heat_flux_selection (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of Heat_flux_selection

global gkw_gui_flux_select

if (get(hObject,'Value') == get(hObject,'Max'))
   % Checkbox is checked-take approriate action
   gkw_gui_flux_select(2) = 1;
else
   % Checkbox is not checked-take approriate action
   gkw_gui_flux_select(2) = 0;
end
axes(handles.Figure1); 
redo_the_plot; 

% --- Executes on button press in Mom_flux_selection.
function Mom_flux_selection_Callback(hObject, eventdata, handles)
% hObject    handle to Mom_flux_selection (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of Mom_flux_selection

global gkw_gui_flux_select

if (get(hObject,'Value') == get(hObject,'Max'))
   % Checkbox is checked-take approriate action
   gkw_gui_flux_select(3) = 1; 
else
   % Checkbox is not checked-take approriate action
   gkw_gui_flux_select(3) = 0;
end
axes(handles.Figure1)
redo_the_plot; 

% --- Executes on button press in Select_species_1.
function Select_species_1_Callback(hObject, eventdata, handles)
% hObject    handle to Select_species_1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

global gkw_gui_nspecies
global gkw_gui_species_select

% Hint: get(hObject,'Value') returns toggle state of Select_species_1
if (gkw_gui_nspecies >=1)
  if (get(hObject,'Value') == 1) 
    gkw_gui_species_select(1) = 1; 
  else 
    gkw_gui_species_select(1) = 0;
    set(handles.Select_species_all,'Value',0);
  end
  axes(handles.Figure1); 
  redo_the_plot; 
end  
 
    

% --- Executes on button press in Select_species_2.
function Select_species_2_Callback(hObject, eventdata, handles)
% hObject    handle to Select_species_2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of Select_species_2

global gkw_gui_nspecies
global gkw_gui_species_select

if (gkw_gui_nspecies >=2) 
  if (get(hObject,'Value') == 1) 
    gkw_gui_species_select(2) = 1;
  else
    gkw_gui_species_select(2) = 0;
    set(handles.Select_species_all,'Value',0);
  end
  axes(handles.Figure1);
  redo_the_plot;
end 

% --- Executes on button press in Select_species_3.
function Select_species_3_Callback(hObject, eventdata, handles)
% hObject    handle to Select_species_3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of Select_species_3
global gkw_gui_nspecies
global gkw_gui_species_select 

if (gkw_gui_nspecies >=3) 
  if (get(hObject,'Value') == 1) 
    gkw_gui_species_select(3) = 1; 
  else
    gkw_gui_species_select(3) = 0; 
    set(handles.Select_species_all,'Value',0); 
  end
  axes(handles.Figure1);
  redo_the_plot;
end

% --- Executes on button press in Select_species_all.
function Select_species_all_Callback(hObject, eventdata, handles)
% hObject    handle to Select_species_all (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of Select_species_all
global gkw_gui_nspecies
global gkw_gui_species_select

if (get(hObject,'Value') == 1) 
  for i = 1 : gkw_gui_nspecies 
    gkw_gui_species_select(i) = 1; 
  end
  set(handles.Select_species_1,'Value',1); 
  set(handles.Select_species_2,'Value',1); 
  set(handles.Select_species_3,'Value',1); 
else 
  for i = 4 : gkw_gui_nspecies 
    gkw_gui_species_select(i) = 0; 
  end 
end 
axes(handles.Figure1);
redo_the_plot;

% --- Executes on button press in Eigenmode.
function Eigenmode_Callback(hObject, eventdata, handles)
% hObject    handle to Eigenmode (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDAT)

global gkw_gui_current_plot 

gkw_gui_current_plot = 'parallel';
axes(handles.Figure1);
redo_the_plot; 



% --- Executes on button press in Select_dens.
function Select_dens_Callback(hObject, eventdata, handles)
% hObject    handle to Select_dens (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of Select_dens
global gkw_gui_par_select

if (get(hObject,'Value') == 1) 
  gkw_gui_par_select(4) = 1;
else 
  gkw_gui_par_select(4) = 0; 
end
axes(handles.Figure1);
redo_the_plot; 


% --- Executes on button press in Select_phi.
function Select_phi_Callback(hObject, eventdata, handles)
% hObject    handle to Select_phi (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of Select_phi
global gkw_gui_par_select

if (get(hObject,'Value') == 1)
  gkw_gui_par_select(1) = 1;
else 
  gkw_gui_par_select(1) = 0; 
end
axes(handles.Figure1);
redo_the_plot; 

% --- Executes on button press in Select_Apar.
function Select_Apar_Callback(hObject, eventdata, handles)
% hObject    handle to Select_Apar (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of Select_Apar
global gkw_gui_par_select

if (get(hObject,'Value') == 1)
  gkw_gui_par_select(2) = 1;
else 
  gkw_gui_par_select(2) = 0; 
end
axes(handles.Figure1);
redo_the_plot; 

% --- Executes on button press in Select_Bpar.
function Select_Bpar_Callback(hObject, eventdata, handles)
% hObject    handle to Select_Bpar (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of Select_Bpar
global gkw_gui_par_select

if (get(hObject,'Value') == 1)
  gkw_gui_par_select(3) = 1;
else 
  gkw_gui_par_select(3) = 0; 
end
axes(handles.Figure1);
redo_the_plot; 

% --- Executes on button press in Select_Tpar.
function Select_Tpar_Callback(hObject, eventdata, handles)
% hObject    handle to Select_Tpar (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of Select_Tpar
global gkw_gui_par_select

if (get(hObject,'Value') == 1) 
  gkw_gui_par_select(5) = 1;
else 
  gkw_gui_par_select(5) = 0; 
end
axes(handles.Figure1);
redo_the_plot; 

% --- Executes on button press in Select_Tperp.
function Select_Tperp_Callback(hObject, eventdata, handles)
% hObject    handle to Select_Tperp (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of Select_Tperp
global gkw_gui_par_select

if (get(hObject,'Value') == 1) 
  gkw_gui_par_select(6) = 1;
else 
  gkw_gui_par_select(6) = 0; 
end
axes(handles.Figure1);
redo_the_plot; 

% --- Executes on button press in Select_vpar.
function Select_vpar_Callback(hObject, eventdata, handles)
% hObject    handle to Select_vpar (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of Select_vpar
global gkw_gui_par_select

if (get(hObject,'Value') == 1) 
  gkw_gui_par_select(7) = 1;
else 
  gkw_gui_par_select(7) = 0; 
end
axes(handles.Figure1);
redo_the_plot; 


%-----------------------------------------------------------------------
% This function reads all the data 
%-----------------------------------------------------------------------
function read_all_data(handles) 

global gkw_gui_directory_name
global gkw_gui_filename 
global gkw_gui_time 
global gkw_gui_fluxes
global gkw_gui_parallel 
global gkw_gui_nspecies
global gkw_gui_species_select

% Check first what directory structure is used 
name = [gkw_gui_directory_name '/input/' gkw_gui_filename]; 
alf = dir(name);
if (isempty(alf))
  name = [gkw_gui_directory_name '/input.dat'];
  alf = dir(name); 
  if (isempty(alf)) 
    disp('Found no valid imput name. Can not read data');
    one_directory = -1; 
  else 
    one_directory = 1; 
  end
else
  one_directory = 0; 
end 

if (one_directory >= 0) 

  if (one_directory ==0) 
    name = [gkw_gui_directory_name '/time/' gkw_gui_filename];  
  else 
    name = [gkw_gui_directory_name '/time.dat'];
  end   
  alf = dir(name); 
  if (isempty(alf)) 
    disp('Could not load the data of growth rate / freq. Looked for:');
    disp(name);
  else 
    gkw_gui_time = load(name); 
  end 

  if (one_directory ==0) 
    name = [gkw_gui_directory_name '/fluxes/' gkw_gui_filename]; 
  else 
    name = [gkw_gui_directory_name '/fluxes.dat'];
  end
  alf = dir(name); 
  if (isempty(alf)) 
    disp('Could not load the data of the fluxes. Looked for:');
    disp(name);
  else
    gkw_gui_fluxes = load(name); 
  end 

  if (one_directory == 0) 
    name = [gkw_gui_directory_name '/parallel/' gkw_gui_filename];
  else 
    name = [gkw_gui_directory_name '/parallel.dat'];
  end 
  alf = dir(name); 
  if (isempty(alf))
    disp('Could not load the parallel mode structure. Looked for');
    disp(name);
  else 
    gkw_gui_parallel = load(name); 
  end 
  
  %Determine the number of species 
  gkw_gui_nspecies = floor(size(gkw_gui_fluxes,2) / 3);
  if (gkw_gui_nspecies == 0) 
    disp('Found zero species. Output can be incorrect');
  end
  
  if  (get(handles.Select_species_1,'Value')==1)
    gkw_gui_species_select(1) = 1;
  else
    gkw_gui_species_select(1) = 0;
  end; 
  if  (get(handles.Select_species_2,'Value')==1)
    gkw_gui_species_select(2) = 1;
  else
    gkw_gui_species_select(2) = 0;
  end; 
  if  (get(handles.Select_species_3,'Value')==1)
    gkw_gui_species_select(3) = 1;
  else
    gkw_gui_species_select(3) = 0;
  end; 
  if  (get(handles.Select_species_all,'Value')==1)
    for i = 1: gkw_gui_nspecies
      gkw_gui_species_select(i) = 1;
    end
  else
    for i = 1: gkw_gui_nspecies
      gkw_gui_species_select(i) = 0;
    end;
  end; 
  
end  

%-----------------------------------------------------------------------
% The general plotting function (calls the appropriate routine 
%-----------------------------------------------------------------------
function redo_the_plot 

global gkw_gui_current_plot 

% Select the right figure

switch gkw_gui_current_plot
case('none') 
  cla;
case('fluxes')
  plot_the_fluxes 
case('time')
  plot_the_time
case('parallel')
  plot_the_mode 
end 

%-----------------------------------------------------------------------
% Function that plots the fluxes 
%-----------------------------------------------------------------------
function plot_the_fluxes

global gkw_gui_nspecies
global gkw_gui_fluxes
global gkw_gui_flux_select
global gkw_gui_species_select
global gkw_gui_time

cla; 
for i = 1 : 3 
  switch i
  case(1) 
    c = 'b-';
  case(2) 
    c = 'g-';
  case(3)
    c = 'r-';
  end
  for j = 1 : gkw_gui_nspecies 
    if (gkw_gui_flux_select(i)*gkw_gui_species_select(j)) 
      plot(gkw_gui_time(:,1),gkw_gui_fluxes(:,3*(j-1)+i),c); 
      hold on;
      ns = size(gkw_gui_time,1);
      is = floor(0.8*ns);
      text(gkw_gui_time(is,1),gkw_gui_fluxes(is,3*(j-1)+i),int2str(j));
    end 
  end 
end 
xlabel('Time [R/v_{thref}]'); 
ylabel('Fluxes [GKW units]'); 
hold off;

%------------------------------------------------------------------------
% Function that plots the growth rate and frequency
%------------------------------------------------------------------------
function plot_the_time 

global gkw_gui_time 
global growth_rate 
global frequency 

% Initialize 
growth_rate = 0;
frequency = 0;

% clean the plot 
cla; 

% size of the array 
ns = size(gkw_gui_time,2); 
ntp = size(gkw_gui_time,1);

if (ns>=2) 
  plot(gkw_gui_time(:,1),gkw_gui_time(:,2),'r-'); 
  growth_rate = gkw_gui_time(ntp,2);
  ym1 = max(gkw_gui_time(floor(ntp/2):ntp,2));
  ym2 = min(gkw_gui_time(floor(ntp/2):ntp,2));
  hold on; 
  if (ns>=3) 
    plot(gkw_gui_time(:,1),gkw_gui_time(:,3),'b-');
    frequency = gkw_gui_time(ntp,3);
    yz1 = max(gkw_gui_time(floor(ntp/2):ntp,3));
    ym1 = max(ym1,yz1);
    yz2 = min(gkw_gui_time(floor(ntp/2):ntp,3));
    ym2 = max(ym2,yz2);
  end
  % Try to get the axis more reasonable 
  ymax = max([0 1.2*growth_rate 1.2*frequency ym1]);
  ymin = min([0 1.2*growth_rate 1.2*frequency ym2]); 
  axis([0 gkw_gui_time(ntp,1) ymin ymax]); 
  xloc = 0.1*gkw_gui_time(ntp,1); 
  yloc = ymin + 0.9 *(ymax-ymin); 
  text(xloc,yloc,['Growth rate :' num2str(growth_rate)]);
  yloc = ymin + 0.85*(ymax-ymin); 
  text(xloc,yloc,['Frequency   :' num2str(frequency)]);
  xlabel('Time [R/v_{thref}]');
  ylabel('Growth rate / Freq. [v_{thref}/R]')
  hold off; 
end

%-------------------------------------------------------------------------
% Function that plots the mode structure
%-------------------------------------------------------------------------
function plot_the_mode 

global gkw_gui_nspecies
global gkw_gui_parallel 
global gkw_gui_par_select
global gkw_gui_species_select

cla; 

% The length of the array with one species 
nlen = size(gkw_gui_parallel,1) / gkw_gui_nspecies; 

% unrole s 
sadd = 0; 
sss(1) = gkw_gui_parallel(1,1);
for i = 1: nlen - 1
  if (abs(gkw_gui_parallel(i+1,1)-gkw_gui_parallel(i,1)) > 0.5)
    sadd = sadd + 1; 
  end
  sss(i+1) = gkw_gui_parallel(i+1,1) + sadd; 
end
val = 0.5*(sss(1) + sss(nlen));
for i = 1:nlen
  sss(i) = sss(i) - val;
end

% Loop over parallel mode structures 
for i = 1 : 7

  % switched on ?? 
  if (gkw_gui_par_select(i))
      
    %set the colour 
    switch(i)
    case(1)
      %potential in red 
      c = 'r';
      i1 = 2; 
      i2 = 3; 
    case(2)
      c = 'b';
      i1 = 4;
      i2 = 5; 
    case(3)
      c = 'g';
      i1 = 0; 
      i2 = 0; 
    case(4)
      c = 'k';
      i1 = 6;
      i2 = 7;
    case(5)
      c = 'c';
      i1 = 8;
      i2 = 9;
    case(6)
      c= 'y';
      i1 = 10;
      i2 = 11;
    case(7)
      c = 'm';
      i1 = 12;
      i2 = 13;
    end    
  
    if (i < 4) 
      nb = 1;
      ne = nlen;
      if (i1 ~= 0)
        plot(sss,gkw_gui_parallel(nb:ne,i1),[c '-']);
        hold on;
        plot(sss,gkw_gui_parallel(nb:ne,i2),[c '--']);
      end    
    else
      for j = 1: gkw_gui_nspecies 
        if (gkw_gui_species_select(j))
          nb = (j-1)*nlen+1;
          ne = j*nlen;
          plot(sss,gkw_gui_parallel(nb:ne,i1),[c '-']);
          hold on;
          plot(sss,gkw_gui_parallel(nb:ne,i2),[c '--']);
          text(sss(floor(nlen/2)),gkw_gui_parallel(floor((ne+nb)/2),i1),int2str(j))
          text(sss(floor(nlen/2)),gkw_gui_parallel(floor((nb+ne)/2),i2),int2str(j))
        end
      end
    end
  end
end
xlabel('s [GKW units]');
ylabel('Eigenfunction [GKW units]');
hold off;
