function varargout = choisirMOG(varargin)
% CHOISIRMOG M-file for choisirMOG.fig
%      CHOISIRMOG, by itself, creates a new CHOISIRMOG or raises the existing
%      singleton*.
%
%      H = CHOISIRMOG returns the handle to a new CHOISIRMOG or the handle to
%      the existing singleton*.
%
%      CHOISIRMOG('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in CHOISIRMOG.M with the given input arguments.
%
%      CHOISIRMOG('Property','Value',...) creates a new CHOISIRMOG or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before choisirMOG_OpeningFunction gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to choisirMOG_OpeningFcn via varargin.
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

% Edit the above text to modify the response to help choisirMOG

% Last Modified by GUIDE v2.5 17-Jan-2011 14:49:56

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @choisirMOG_OpeningFcn, ...
                   'gui_OutputFcn',  @choisirMOG_OutputFcn, ...
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


% --- Executes just before choisirMOG is made visible.
function choisirMOG_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to choisirMOG (see VARARGIN)

% Choose default command line output for choisirMOG
handles.output = 0;
handles.second_output = '';
handles.third_output = hObject;

db_file = get(handles.choix_MOG,'UserData');
if ~isempty(db_file)
	handles.second_output = db_file;
	file = getDBname(db_file);
	set(handles.text_db,'String',file)
	load(db_file,'names_mog')
	set(handles.popupmenu_mogs,'String',names_mog)
end

% Update handles structure
guidata(hObject, handles);

str = get_str_locale();
set_String_locale(handles, str)

% UIWAIT makes choisirMOG wait for user response (see UIRESUME)
uiwait(handles.choix_MOG);


% --- Outputs from this function are returned to the command line.
function varargout = choisirMOG_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;
varargout{2} = handles.second_output;
varargout{3} = handles.third_output;


function pushbutton_db_Callback(hObject, eventdata, handles)
[file, rep] = uigetfile('*.mat','Ouvrir DB');
if isequal(file,0)
	return
end
handles.second_output = [rep,file];
set(handles.text_db,'String',file)
load([rep,file],'names_mog')
set(handles.popupmenu_mogs,'String',names_mog)
guidata(hObject, handles);

function popupmenu_mogs_Callback(hObject, eventdata, handles)


function popupmenu_mogs_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


function pushbutton_OK_Callback(hObject, eventdata, handles)
handles.output = get(handles.popupmenu_mogs,'Value');
% Update handles structure
guidata(hObject, handles);
uiresume;


function pushbutton_annuler_Callback(hObject, eventdata, handles)
uiresume;

function set_String_locale(handles, str)
set(handles.pushbutton_db,'String',str.s127)
set(handles.pushbutton_annuler,'String',str.s91)
