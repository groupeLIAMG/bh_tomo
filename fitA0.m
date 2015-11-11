function varargout = fitA0(varargin)
% FITA0 M-file for fitA0.fig
%      FITA0, by itself, creates a new FITA0 or raises the existing
%      singleton*.
%
%      H = FITA0 returns the handle to a new FITA0 or the handle to
%      the existing singleton*.
%
%      FITA0('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in FITA0.M with the given input arguments.
%
%      FITA0('Property','Value',...) creates a new FITA0 or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before fitA0_OpeningFunction gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to fitA0_OpeningFcn via varargin.
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

% Edit the above text to modify the response to help fitA0

% Last Modified by GUIDE v2.5 13-May-2005 16:18:47

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @fitA0_OpeningFcn, ...
                   'gui_OutputFcn',  @fitA0_OutputFcn, ...
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



% --- Executes just before fitA0 is made visible.
function fitA0_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to fitA0 (see VARARGIN)

% Choose default command line output for fitA0
handles.output = 0;
handles.second_output = hObject;
handles.third_output = 0;

% Update handles structure
guidata(hObject, handles);

h2 = varargin{1};
h1.tau = h2.tau;
h1.m = h2.m;
h1.A0initial = h2.m(2);
h1.mInitial = h2.m(1);
h2.sigma_s2 = median(h2.var);
set(handles.edit_A0,'String',num2str(h1.m(2)))
set(handles.edit_pente,'String',num2str(h1.m(1)))
setappdata(handles.fig_A0, 'h1', h1)
setappdata(handles.fig_A0, 'h2', h2)

str = get_str_locale();
set_String_locale(handles, str)
update_plots(handles)

% UIWAIT makes fitA0 wait for user response (see UIRESUME)
uiwait(handles.fig_A0);



% --- Outputs from this function are returned to the command line.
function varargout = fitA0_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;
varargout{2} = handles.second_output;
varargout{3} = handles.third_output;

function update_plots(handles)
h1 = getappdata(handles.fig_A0, 'h1');
h2 = getappdata(handles.fig_A0, 'h2');

axes(handles.axes1)
plot(h2.lrai, h2.rawdata, 'o')
ylabel(h2.label{1})
xlabel(h2.label{6})

axes(handles.axes2)
plot(h2.lrai, h2.var, 'o')
ylabel(h2.label{2})
xlabel(h2.label{6})

axes(handles.axes3)
plot(h2.lrai, h2.data, 'o')
xlim=get(gca,'XLim');
hold on
plot(xlim,[xlim' [1 1]']*h1.m,'g-')
hold off
ylabel(h2.label{3})
xlabel(h2.label{6})

axes(handles.axes4)
plot(h2.lrai, (h1.tau)./h2.lrai, 'o')
ylabel(h2.label{4})
xlabel(h2.label{6})

axes(handles.axes5)
plot(h2.lrai, h1.tau, 'o')
ylabel(h2.label{5})
xlabel(h2.label{6})


function edit_A0_Callback(hObject, eventdata, handles)
h1 = getappdata(handles.fig_A0, 'h1');
h2 = getappdata(handles.fig_A0, 'h2');
if strcmp(h2.type,'amp')==1
	h1.m(2) = str2double(get(hObject,'String'));
	h1.tau = h1.m(2) - h2.data;
else
	h1.m(2) = 1e6*h2.f_Tx*str2double(get(hObject,'String'))/h2.sigma_s2;
	data = h2.f_Tx * h2.rawdata / h2.sigma_s2;
	h1.tau = h1.m(2) - data;
end
setappdata(handles.fig_A0, 'h1', h1)
update_plots(handles)


function edit_A0_CreateFcn(hObject, eventdata, handles)

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


function edit_pente_Callback(hObject, eventdata, handles)
h1 = getappdata(handles.fig_A0, 'h1');
h1.m(1) = str2double(get(hObject,'String'));
setappdata(handles.fig_A0, 'h1', h1)
update_plots(handles)


function edit_pente_CreateFcn(hObject, eventdata, handles)

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


function pushbutton_OK_Callback(hObject, eventdata, handles)
h1 = getappdata(handles.fig_A0, 'h1');
handles.output = h1.m(2);
handles.third_output = h1.m(1);


% Update handles structure
guidata(hObject, handles);
uiresume;


function pushbutton_cancel_Callback(hObject, eventdata, handles)
h1 = getappdata(handles.fig_A0, 'h1');
handles.output = h1.A0initial;
handles.third_output = h1.mInitial;

% Update handles structure
guidata(hObject, handles);
uiresume;


function set_String_locale(handles, str)
set(handles.uipanel1,              'Title',  str.s89);
set(handles.text2,                 'String', str.s90);
set(handles.pushbutton_cancel,     'String', str.s91);
