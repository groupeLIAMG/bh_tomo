function varargout = bh_tomo(varargin)
% BH_TOMO M-file for bh_tomo.fig
%      BH_TOMO, by itself, creates a new BH_TOMO or raises the existing
%      singleton*.
%
%      H = BH_TOMO returns the handle to a new BH_TOMO or the handle to
%      the existing singleton*.
%
%      BH_TOMO('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in BH_TOMO.M with the given input arguments.
%
%      BH_TOMO('Property','Value',...) creates a new BH_TOMO or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before bh_tomo_OpeningFunction gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to bh_tomo_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Copyright (C) 2005 Bernard Giroux
%
% This program is free software; you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation; either version 3 of the License, or
% (at your option) any later version.
%
% This program is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details.
% 
% You should have received a copy of the GNU General Public License
% along with this program.  If not, see <http://www.gnu.org/licenses/>.
% 
% 

% Edit the above text to modify the response to help bh_tomo

% Last Modified by GUIDE v2.5 27-May-2015 22:46:47

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @bh_tomo_OpeningFcn, ...
                   'gui_OutputFcn',  @bh_tomo_OutputFcn, ...
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

% --- Executes just before bh_tomo is made visible.
function bh_tomo_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to bh_tomo (see VARARGIN)

% Choose default command line output for bh_tomo
handles.output = hObject;
handles.init_data = 0;
% Update handles structure
guidata(hObject, handles);

str = get_str_locale();
set_String_locale(handles, str)

db_file = [];
setappdata(handles.fig_bh_tomo,'db_file',db_file)

%add_path;
%get(handles.FileMenuItem,'Parent')

% UIWAIT makes bh_tomo wait for user response (see UIRESUME)
% uiwait(handles.fig_bh_tomo);

% --- Outputs from this function are returned to the command line.
function varargout = bh_tomo_OutputFcn(hObject, eventdata, handles)
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;

function pointe_tt_Callback(hObject, eventdata, handles)
db_file = getappdata(handles.fig_bh_tomo,'db_file');
bh_tomo_tt( 'UserData', db_file )

function pointe_amp_Callback(hObject, eventdata, handles)
db_file = getappdata(handles.fig_bh_tomo,'db_file');
bh_tomo_amp( 'UserData', db_file )

function inversion_Callback(hObject, eventdata, handles)
db_file = getappdata(handles.fig_bh_tomo,'db_file');
bh_tomo_inv( 'UserData', db_file )

function fitCovar_Callback(hObject, eventdata, handles)
db_file = getappdata(handles.fig_bh_tomo,'db_file');
bh_tomo_fitCovar( 'UserData', db_file )

function gere_db_Callback(hObject, eventdata, handles)
db_file = getappdata(handles.fig_bh_tomo,'db_file');
bh_tomo_db( 'UserData', db_file );

function interpretation_Callback(hObject, eventdata, handles)
db_file = getappdata(handles.fig_bh_tomo,'db_file');
bh_tomo_interp( 'UserData', db_file );


function pick_Callback(hObject, eventdata, handles)
db_file = getappdata(handles.fig_bh_tomo,'db_file');
bh_tomo_pick( 'UserData', db_file )


function phase_pick_Callback(hObject, eventdata, handles)
db_file = getappdata(handles.fig_bh_tomo,'db_file');
bh_tomo_phase_pick( 'UserData', db_file )

%%%%%added 2012-10-   YH
function pushbutton_timelapse_Callback(hObject, eventdata, handles)
db_file = getappdata(handles.fig_bh_tomo,'db_file');
bh_tomo_timelapse( 'UserData', db_file )


function fitCovar3d_Callback(hObject, eventdata, handles)
db_file = getappdata(handles.fig_bh_tomo,'db_file');
bh_tomo_fitCovar3d( 'UserData', db_file )

function inversion3d_Callback(hObject, eventdata, handles)
db_file = getappdata(handles.fig_bh_tomo,'db_file');
bh_tomo_inv3d( 'UserData', db_file )


function pushbutton_timelapse3d_Callback(hObject, eventdata, handles)
db_file = getappdata(handles.fig_bh_tomo,'db_file');
bh_tomo_timelapse3d( 'UserData', db_file )

%%%%%%%%%%%%%%%%%%%%%%%
function set_String_locale(handles, str)
set(handles.titre,'String',str.s123)
set(handles.gere_db,'String',str.s124)
set(handles.pick,'String',str.s273)
set(handles.phase_pick,'String',str.s282)
set(handles.pointe_tt,'String',str.s125)
set(handles.pointe_amp,'String',str.s126)
set(handles.fitCovar,'String',str.s93)
set(handles.inversion,'String',str.s33)
set(handles.FileMenuItem,'Label',str.s25)
set(handles.OpenMenuItem,'Label',[str.s127,' ...'])
set(handles.QuitMenuItem,'Label',str.s31)
set(handles.HelpMenuItem,'Label',str.s32)

function FileMenuItem_Callback(hObject, eventdata, handles)


function QuitMenuItem_Callback(hObject, eventdata, handles)
delete(handles.fig_bh_tomo)
 

function OpenMenuItem_Callback(hObject, eventdata, handles)
loc_str = get_str_locale();
[file, rep] = uigetfile('*.mat',loc_str.s127);
if isequal(file,0)
	return
end
db_file = [rep,file];
name = getDBname(db_file);
set(handles.text_db_name,'String',name)
setappdata(handles.fig_bh_tomo,'db_file',db_file)

%%%%%%%%%%%%added for changing names of variables 2012-11-15  Yh
load(db_file);
% for i = 1:length(air)
%     air(i).name = air(i).nom;
%     air(i).tt_done = air(i).tt_fait;
%     
% end
% air = rmfield(air,'nom');
% air = rmfield(air,'tt_fait');
% for i = 1:length(panels)
% 
%     for j = 1:length(panels(i).inv_res)
% % if ~isempty(panels(i).amp_covar)
%         panels(i).inv_res(j).covar.model = panels(i).inv_res(j).covar.modele;
%     panels(i).inv_res(j).covar.nugget_t = panels(i).inv_res(j).covar.pepite_t;
%     panels(i).inv_res(j).covar.nugget_l = panels(i).inv_res(j).covar.pepite_l;
%     
%     panels(i).inv_res(j).covar = rmfield(panels(i).inv_res(j).covar,'modele');
%     panels(i).inv_res(j).covar = rmfield(panels(i).inv_res(j).covar,'pepite_t');
%     panels(i).inv_res(j).covar = rmfield(panels(i).inv_res(j).covar,'pepite_l');
%     end
% end
% for i = 1:length(mogs)
%     mogs(i).amp_name_Ldc = mogs(i).amp_nom_Ldc;
% end
%  mogs = rmfield(mogs,'amp_nom_Ldc');
% for i = 1:length(boreholes)
%     boreholes(i).name = boreholes(i).nom;
%     boreholes(i).Z_water = boreholes(i).Z_eau;
% end
% boreholes = rmfield(boreholes,'nom');
% boreholes = rmfield(boreholes,'Z_eau');
%  save(db_file);

verify_struct_fields(panels,mogs,boreholes);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


function HelpMenuItem_Callback(hObject, eventdata, handles)
dos('Help.chm');


% --- Executes on button press in pushbutton_nanofluid.
function pushbutton_nanofluid_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton_nanofluid (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
db_file = getappdata(handles.fig_bh_tomo,'db_file');
bh_tomo_nanofluid( 'UserData', db_file )


% --- Executes on button press in pushbutton_tltomo2d.
function pushbutton_tltomo2d_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton_tltomo2d (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
db_file = getappdata(handles.fig_bh_tomo,'db_file');
bh_tomo_tltomo2d( 'UserData', db_file )


% --- Executes on button press in pushbutton_tltomo3d.
function pushbutton_tltomo3d_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton_tltomo3d (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
db_file = getappdata(handles.fig_bh_tomo,'db_file');
bh_tomo_tltomo3d( 'UserData', db_file )
