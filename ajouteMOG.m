function varargout = ajouteMOG(varargin)
% AJOUTEMOG M-file for ajouteMOG.fig
%      AJOUTEMOG, by itself, creates a new AJOUTEMOG or raises the existing
%      singleton*.
%
%      H = AJOUTEMOG returns the handle to a new AJOUTEMOG or the handle to
%      the existing singleton*.
%
%      AJOUTEMOG('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in AJOUTEMOG.M with the given input arguments.
%
%      AJOUTEMOG('Property','Value',...) creates a new AJOUTEMOG or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before ajouteMOG_OpeningFunction gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to ajouteMOG_OpeningFcn via varargin.
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


% Edit the above text to modify the response to help ajouteMOG

% Last Modified by GUIDE v2.5 06-Jul-2005 13:13:05

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @ajouteMOG_OpeningFcn, ...
                   'gui_OutputFcn',  @ajouteMOG_OutputFcn, ...
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


% --- Executes just before ajouteMOG is made visible.
function ajouteMOG_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to ajouteMOG (see VARARGIN)

% Choose default command line output for ajouteMOG
handles.output = 0;
handles.second_output = hObject;

% Update handles structure
guidata(hObject, handles);

if nargin ~= 4
	uiwait(warndlg('Oh oh!'))
	delete(hObject)
else
	str = varargin{1};
	set(handles.popupmenu_MOG,'String',str)
end

str = get_str_locale();
set(handles.text_selection,    'String', str.s182)
set(handles.pushbutton_ajoute, 'String', str.s128)
set(handles.pushbutton_annule, 'String', str.s91)



% UIWAIT makes ajouteMOG wait for user response (see UIRESUME)
uiwait(handles.fig_ajouteMOG);


% --- Outputs from this function are returned to the command line.
function varargout = ajouteMOG_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;
varargout{2} = handles.second_output;


function pushbutton_ajoute_Callback(hObject, eventdata, handles)
handles.output = get(handles.popupmenu_MOG,'Value');

% Update handles structure
guidata(hObject, handles);
uiresume;


function pushbutton_annule_Callback(hObject, eventdata, handles)
uiresume;


function popupmenu_MOG_Callback(hObject, eventdata, handles)


function popupmenu_MOG_CreateFcn(hObject, eventdata, handles)

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


function edit_poids_Callback(hObject, eventdata, handles)


function edit_poids_CreateFcn(hObject, eventdata, handles)

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


