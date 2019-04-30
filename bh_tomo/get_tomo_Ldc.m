function varargout = get_tomo_Ldc(varargin)
% GET_TOMO_LDC M-file for get_tomo_Ldc.fig
%      GET_TOMO_LDC, by itself, creates a new GET_TOMO_LDC or raises the existing
%      singleton*.
%
%      H = GET_TOMO_LDC returns the handle to a new GET_TOMO_LDC or the handle to
%      the existing singleton*.
%
%      GET_TOMO_LDC('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in GET_TOMO_LDC.M with the given input arguments.
%
%      GET_TOMO_LDC('Property','Value',...) creates a new GET_TOMO_LDC or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before get_tomo_Ldc_OpeningFunction gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to get_tomo_Ldc_OpeningFcn via varargin.
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

% Edit the above text to modify the response to help get_tomo_Ldc

% Last Modified by GUIDE v2.5 17-Dec-2012 14:54:48

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @get_tomo_Ldc_OpeningFcn, ...
                   'gui_OutputFcn',  @get_tomo_Ldc_OutputFcn, ...
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


% --- Executes just before get_tomo_Ldc is made visible.
function get_tomo_Ldc_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to get_tomo_Ldc (see VARARGIN)

% Choose default command line output for get_tomo_Ldc
handles.output = [];
handles.second_output = {};
handles.third_output = hObject;

% Update handles structure
guidata(hObject, handles);

str = get_str_locale();
set_String_locale(handles, str)

db_file = '';
if nargin >=4
	db_file = char( varargin{1} );
end
if ~isempty(db_file)
	tmp = findstr('/',db_file);
	if ~isempty(tmp)
		file = db_file((tmp(length(tmp))+1):length(db_file));
	else
		tmp = findstr('\',db_file); %windoze
		if ~isempty(tmp)
			file = db_file((tmp(length(tmp))+1):length(db_file));
		else
			file = db_file;
		end
	end
	set(handles.text_db,'String',file)

%     load(db_file,'models')
%     noms_models = {};
%     for n=1:length(models)
%         noms_models{n} = char( models(n).name );
%     end
%     set(handles.popupmenu_pannneaux,'String',noms_models)
% 
%     if isfield(models(1), 'inv_res')
%         noms_inv_res = {};
%         for n=1:length(models(1).inv_res)
%             noms_inv_res{n} = char( models(1).inv_res(n).name );
%         end
%         set(handles.popupmenu_inv_res,'String',noms_inv_res)        
%     end
%%%%%%%%%%%%%%%%%%%%%%  YH
    noms_models = {};
    noms_inv_res = {};
  	load(db_file,'models')
    nb_models = length(models);
    if nb_models > 0
        for n=1:nb_models
            noms_models{n} = char(models(n).name );
        end
        if isfield(models(1), 'inv_res') || isprop(models(1), 'inv_res')
            for n=1:length(models(1).inv_res)
                noms_inv_res{n} = char(models(1).inv_res(n).name );
            end
        end
    end
    if isempty(noms_models{1})
        set(handles.popupmenu_pannneaux,'Value',1,'String','-');
    else
        set(handles.popupmenu_pannneaux,'Value',1,'String',noms_models);
    end
    
    if isempty(noms_inv_res{1})
        set(handles.popupmenu_inv_res,'Value',1,'String','-'); 
    else
        set(handles.popupmenu_inv_res,'Value',1,'String',noms_inv_res);        
    end    

end
%%%%%%%%%%%%%%%%%%%%%%%%%%%
setappdata(handles.fig_get_Ldc,'db_file',db_file)

% UIWAIT makes get_tomo_Ldc wait for user response (see UIRESUME)
uiwait(handles.fig_get_Ldc);

% --- Outputs from this function are returned to the command line.
function varargout = get_tomo_Ldc_OutputFcn(hObject, eventdata, handles) 
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
set(handles.text_db,'String',file)
db_file = [rep,file];

% load([rep,file],'models')
% noms_models = {};
% for n=1:length(models)
%     noms_models{n} = char( models(n).name );
% end
% set(handles.popupmenu_pannneaux,'String',noms_models)
% if isfield(models(1), 'inv_res')
%     noms_inv_res = {};
%     for n=1:length(models(1).inv_res)
%         noms_inv_res{n} = char( models(1).inv_res(n).name );
%     end
%     set(handles.popupmenu_inv_res,'String',noms_inv_res)
% else
%     %warning
% end
%%%%%%%%%%%%%%%%%%%% YH     
noms_models = {};
noms_inv_res = {};
load(db_file,'models')
nb_models = length(models);
if nb_models > 0
    for n=1:nb_models
        noms_models{n} = char(models(n).name );
    end
    if isfield(models(1), 'inv_res')        
        for n=1:length(models(1).inv_res)
            noms_inv_res{n} = char(models(1).inv_res(n).name );
        end
    end
end
if isempty(noms_models)
    set(handles.popupmenu_pannneaux,'Value',1,'String','-');
else
    set(handles.popupmenu_pannneaux,'Value',1,'String',noms_models);
end

if isempty(noms_inv_res)
    set(handles.popupmenu_inv_res,'Value',1,'String','-'); 
else
    set(handles.popupmenu_inv_res,'Value',1,'String',noms_inv_res);        
end    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%
setappdata(handles.fig_get_Ldc,'db_file',db_file)
guidata(hObject, handles);

function popupmenu_pannneaux_Callback(hObject, eventdata, handles)
no = get(hObject,'Value');
db_file = getappdata(handles.fig_get_Ldc,'db_file');

% load(db_file,'models')
% if isfield(models(no), 'inv_res')
%     noms_inv_res = {};
%     for n=1:length(models(no).inv_res)
%         noms_inv_res{n} = char( models(no).inv_res(n).name );
%     end
%     set(handles.popupmenu_inv_res,'String',noms_inv_res)
% else
%     str = get_str_locale();
%     uiwait(warndlg(str.s183))
% end
%%%%%%%%%%%%%%%%%%%% YH
 
load(db_file,'models')
nb_models = length(models);
noms_inv_res = {};
if nb_models > 0
    if no <= nb_models
        if isfield(models(no), 'inv_res')
            for n=1:length(models(no).inv_res)
                noms_inv_res{n} = char( models(no).inv_res(n).name );
            end
        else
            str = get_str_locale();
            uiwait(warndlg(str.s183))
        end
    end
end
if isempty(noms_inv_res)
    set(handles.popupmenu_inv_res,'Value',1,'String','-'); 
else
    set(handles.popupmenu_inv_res,'Value',1,'String',noms_inv_res);        
end 
%%%%%%%%%%%%%%%%%%%%%

function popupmenu_pannneaux_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function pushbutton_OK_Callback(hObject, eventdata, handles)
no1 = get(handles.popupmenu_pannneaux,'Value');
no2 = get(handles.popupmenu_inv_res,'Value');
db_file = getappdata(handles.fig_get_Ldc,'db_file');

% load(db_file,'models')
% handles.output = models(no1).inv_res(no2).tomo;
% handles.second_output = {models(no1).name; models(no1).inv_res(no2).name};

%%%%%%%%%%%%%%%%%YH
load(db_file,'models')
nb_models = length(models);
if nb_models > 0
    if no1 <= nb_models
        if isfield(models(no1), 'inv_res')
            if isempty(models(no1).inv_res)
                handles.output = [];
                handles.second_output = {models(no1).name; '?'};
            else
                handles.output = models(no1).inv_res(no2).tomo;
                handles.second_output = {models(no1).name; models(no1).inv_res(no2).name};
            end
        end
    end
end
%%%%%%%%%%%%%
% Update handles structure
guidata(hObject, handles);
uiresume;


function pushbutton_annuler_Callback(hObject, eventdata, handles)
uiresume;

function set_String_locale(handles, str)
set(handles.pushbutton_db,'String',str.s127)
set(handles.pushbutton_annuler,'String',str.s91)


function popupmenu_inv_res_Callback(hObject, eventdata, handles)


function popupmenu_inv_res_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
