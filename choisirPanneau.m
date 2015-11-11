function varargout = choisirPanneau(varargin)
% CHOISIRPANNEAU M-file for choisirPanneau.fig
%      CHOISIRPANNEAU, by itself, creates a new CHOISIRPANNEAU or raises the existing
%      singleton*.
%
%      H = CHOISIRPANNEAU returns the handle to a new CHOISIRPANNEAU or the handle to
%      the existing singleton*.
%
%      CHOISIRPANNEAU('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in CHOISIRPANNEAU.M with the given input arguments.
%
%      CHOISIRPANNEAU('Property','Value',...) creates a new CHOISIRPANNEAU or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before choisirPanneau_OpeningFunction gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to choisirPanneau_OpeningFcn via varargin.
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

% Edit the above text to modify the response to help choisirPanneau

% Last Modified by GUIDE v2.5 24-May-2013 09:51:26

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @choisirPanneau_OpeningFcn, ...
                   'gui_OutputFcn',  @choisirPanneau_OutputFcn, ...
                   'gui_LayoutFcn',  [] , ...
                   'gui_Callback',   []);
if nargin && ischar(varargin{1})
    if ~isempty(varargin{1}), gui_State.gui_Callback = str2func(varargin{1}); end
end

if nargout
    [varargout{1:nargout}] = gui_mainfcn(gui_State, varargin{:});
else
    gui_mainfcn(gui_State, varargin{:});
end
% End initialization code - DO NOT EDIT


% --- Executes just before choisirPanneau is made visible.
function choisirPanneau_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to choisirPanneau (see VARARGIN)

% Choose default command line output for choisirPanneau
handles.output = 0;
handles.second_output = '';
handles.third_output = hObject;
%%%%YH
if nargin > 5
    handles.model = varargin{3};
else
    handles.model = '2d';
end
guidata(hObject,handles);
%%%%%%%%%%%%%%
db_file = get(handles.choix_panneau, 'UserData');
if ~isempty(db_file)
	handles.second_output = db_file;
	tmp = findstr('/',db_file);
	if ~isempty(tmp)
		file = db_file((tmp(length(tmp))+1):length(db_file));
	else
		tmp = findstr('\',db_file); %windoze
		if ~isempty(tmp)
			file = db_file((tmp(length(tmp))+1):length(db_file));
		else
			file = dbfile;
		end
	end
	set(handles.text_db,'String',file)
	load(db_file,'panels')
	noms_panels = {};
	for n=1:length(panels)
		noms_panels{n} = char( panels(n).name );
    end
    %%%%% add model3d in the list  YH
    if strfind(handles.model,'3d')
    load(db_file,'model3d')
    noms_panels = {};
        for n=1:length(model3d)
            noms_panels{n} = char( model3d(n).name );
        end
    end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%
	set(handles.popupmenu_pannneaux,'String',noms_panels)
end

% Update handles structure
guidata(hObject, handles);

str = get_str_locale();
set_String_locale(handles, str)

% UIWAIT makes choisirPanneau wait for user response (see UIRESUME)
uiwait(handles.choix_panneau);


% --- Outputs from this function are returned to the command line.
function varargout = choisirPanneau_OutputFcn(hObject, eventdata, handles) 
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
load([rep,file],'panels')
noms_panels = {};
for n=1:length(panels)
	noms_panels{n} = char( panels(n).name );
end
%%%%% add model3d in the list  YH
if strfind(handles.model,'3d')
    load([rep,file],'model3d')
    noms_panels = {};
    for n=1:length(model3d)
        noms_panels{n} = char( model3d(n).name );
    end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%
set(handles.popupmenu_pannneaux,'String',noms_panels)
guidata(hObject, handles);

function popupmenu_pannneaux_Callback(hObject, eventdata, handles)


function popupmenu_pannneaux_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


function pushbutton_OK_Callback(hObject, eventdata, handles)
handles.output = get(handles.popupmenu_pannneaux,'Value');
% Update handles structure
guidata(hObject, handles);
uiresume;


function pushbutton_annuler_Callback(hObject, eventdata, handles)
uiresume;

function set_String_locale(handles, str)
set(handles.pushbutton_db,'String',str.s127)
set(handles.pushbutton_annuler,'String',str.s91)


% --- Executes on key press with focus on choix_panneau or any of its controls.
function choix_panneau_WindowKeyPressFcn(hObject, eventdata, handles)
% hObject    handle to choix_panneau (see GCBO)
% eventdata  structure with the following fields (see FIGURE)
%	Key: name of the key that was pressed, in lower case
%	Character: character interpretation of the key(s) that was pressed
%	Modifier: name(s) of the modifier key(s) (i.e., control, shift) pressed
% handles    structure with handles and user data (see GUIDATA)
currChar = get(hObject,'CurrentCharacter');
if isequal(currChar,char(13)) %char(13) == enter key
   pushbutton_OK_Callback(hObject, eventdata, handles)
end
