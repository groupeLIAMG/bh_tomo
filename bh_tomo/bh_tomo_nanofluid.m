function varargout = bh_tomo_nanofluid(varargin)
% BH_TOMO_NANOFLUID MATLAB code for bh_tomo_nanofluid.fig
%      BH_TOMO_NANOFLUID, by itself, creates a new BH_TOMO_NANOFLUID or raises the existing
%      singleton*.
%
%      H = BH_TOMO_NANOFLUID returns the handle to a new BH_TOMO_NANOFLUID or the handle to
%      the existing singleton*.
%
%      BH_TOMO_NANOFLUID('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in BH_TOMO_NANOFLUID.M with the given input arguments.
%
%      BH_TOMO_NANOFLUID('Property','Value',...) creates a new BH_TOMO_NANOFLUID or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before bh_tomo_nanofluid_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to bh_tomo_nanofluid_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help bh_tomo_nanofluid

% Last Modified by GUIDE v2.5 16-Sep-2014 17:27:58

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @bh_tomo_nanofluid_OpeningFcn, ...
                   'gui_OutputFcn',  @bh_tomo_nanofluid_OutputFcn, ...
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


% --- Executes just before bh_tomo_nanofluid is made visible.
function bh_tomo_nanofluid_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to bh_tomo_nanofluid (see VARARGIN)

% Choose default command line output for bh_tomo_nanofluid
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes bh_tomo_nanofluid wait for user response (see UIRESUME)
% uiwait(handles.fig_nanofluid);
h.db_file = get(handles.fig_nanofluid,'UserData');
h.S = [];
h.ivel1 = [];
h.ivel2 = [];
h.iatt1 = [];
h.iatt2 = [];

h.epsilon_m = 8;
h.a = 0.62;
h.m = 1.72;
h.sigma_nf = 0.01;
h.epsilon_nf = 75;
h.mu_nf = 3;
h.epsilon_f1 = 2.25;
h.sigma_f1 = 0.001;
h.omega = 2*pi*1e6*50;
h.S0 = 0.5;

set(handles.edit_epsilon_m,'String', num2str(h.epsilon_m ) );
set(handles.edit_a,'String', num2str(h.a) );
set(handles.edit_m,'String', num2str(h.m) );
set(handles.edit_sigma_nf,'String', num2str(h.sigma_nf) );
set(handles.edit_epsilon_nf,'String', num2str(h.epsilon_nf) );
set(handles.edit_mu_nf,'String', num2str(h.mu_nf) );
set(handles.edit_epsilon_f1,'String', num2str(h.epsilon_f1) );
set(handles.edit_sigma_f1,'String', num2str(h.sigma_f1) );
set(handles.edit_frequency,'String', num2str(h.omega/(2*pi*1e6)) );
set(handles.edit_S0,'String', num2str(h.S0) );

setappdata(handles.fig_nanofluid, 'h', h);


% --- Outputs from this function are returned to the command line.
function varargout = bh_tomo_nanofluid_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


% --- Executes on button press in pushbutton_compute.
function pushbutton_compute_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton_compute (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% check consistency of input data
ivel1 = get(handles.popupmenu_baseline_vel,'Value');
iatt1 = get(handles.popupmenu_baseline_att,'Value');
ivel2 = get(handles.popupmenu_repeat_vel,'Value');
iatt2 = get(handles.popupmenu_repeat_att,'Value');

if ivel1 == ivel2 || iatt1 == iatt2
    msgbox('Repeat data should not equal baseline data','Error','modal');
    return
end
h = getappdata(handles.fig_nanofluid, 'h');

data_v1 = 1./h.panels(h.no_panel).inv_res(ivel1).tomo.s;
data_v2 = 1./h.panels(h.no_panel).inv_res(ivel2).tomo.s;
data_a1 = h.panels(h.no_panel).inv_res(iatt1).tomo.s;
data_a2 = h.panels(h.no_panel).inv_res(iatt2).tomo.s;

h.epsilon_m = str2double( get(handles.edit_epsilon_m,'String') );
h.a = str2double( get(handles.edit_a,'String') );
h.m = str2double( get(handles.edit_m,'String') );
h.sigma_nf = str2double( get(handles.edit_sigma_nf,'String') );
h.epsilon_nf = str2double( get(handles.edit_epsilon_nf,'String') );
h.mu_nf = str2double( get(handles.edit_mu_nf,'String') );
h.epsilon_f1 = str2double( get(handles.edit_epsilon_f1,'String') );
h.sigma_f1 = str2double( get(handles.edit_sigma_f1,'String') );
h.omega = 2*pi*1e6*str2double( get(handles.edit_frequency,'String') );
h.S0 = str2double( get(handles.edit_S0,'String') );

h.S = nan(size(data_v1));
h.ivel1 = ivel1;
h.ivel2 = ivel2;
h.iatt1 = iatt1;
h.iatt2 = iatt2;

hh = waitbar(0,'1','Name','Computing saturation ...','CreateCancelBtn','setappdata(gcf,''canceling'',1)');
setappdata(gcf,'canceling',0)
for n=130%1:length(h.S)
    % Check for Cancel button press
    drawnow
    if getappdata(gcf,'canceling')
        break
    end
    
    waitbar(n/length(h.S), hh, sprintf('Node %d of %d',n, length(h.S)))
    h.S(n) = findS(data_v1(n),data_a1(n),data_v2(n),data_a2(n),...
        h.epsilon_m,h.a,h.m,h.sigma_nf,h.epsilon_nf,h.mu_nf,h.epsilon_f1,...
        h.sigma_f1,h.omega,h.S0);
end
delete (hh)
setappdata(handles.fig_nanofluid,'h',h);
update_fig(handles);

% --------------------------------------------------------------------
function File_menu_Callback(hObject, eventdata, handles)
% hObject    handle to File_menu (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --------------------------------------------------------------------
function Open_menu_item_Callback(hObject, eventdata, handles)
% hObject    handle to Open_menu_item (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
h = getappdata(handles.fig_nanofluid, 'h');
[h.no_panel,h.db_file,hh] = choisirPanneau('UserData', h.db_file);
if ishandle(hh), delete(hh), end
if isequal(h.no_panel,0)
	return
end
load(h.db_file,'panels')
if isempty(panels)
    msgbox('Please define panels and carry out inversion','warning','modal');
    return;
elseif isempty(panels(h.no_panel).inv_res)
    msgbox('Please define panels and carry out inversion','warning','modal');
    return;
else
    h.panels = panels;
    setappdata(handles.fig_nanofluid,'h',h);
    name_panels = {};
    for n = 1:length(panels(h.no_panel).inv_res)
        name_panels{n} = char( panels(h.no_panel).inv_res(n).name );
    end
    set(handles.popupmenu_baseline_vel,'String',name_panels);
    set(handles.popupmenu_baseline_att,'String',name_panels);
    set(handles.popupmenu_repeat_vel,'String',name_panels);
    set(handles.popupmenu_repeat_att,'String',name_panels);

    set(handles.popupmenu_baseline_vel,'Value',1);
    set(handles.popupmenu_baseline_att,'Value',1);
    set(handles.popupmenu_repeat_vel,'Value',1);
    set(handles.popupmenu_repeat_att,'Value',1);
end

if isfield(h.panels(h.no_panel),'saturation')
    name_maps = {};
    for n = 1:length(panels(h.no_panel).saturation)
        name_maps{n} = char( panels(h.no_panel).saturation(n).name );
    end
    set(handles.popupmenu_Smaps,'String',name_maps);
    set(handles.popupmenu_Smaps,'Value',1);
    
    sat = panels(h.no_panel).saturation(1);
    
    set(handles.popupmenu_baseline_vel,'Value',sat.ivel1);
    set(handles.popupmenu_baseline_att,'Value',sat.iatt1);
    set(handles.popupmenu_repeat_vel,'Value',sat.ivel2);
    set(handles.popupmenu_repeat_att,'Value',sat.iatt2);

    set(handles.edit_epsilon_m,'String', num2str(sat.epsilon_m ) );
    set(handles.edit_a,'String', num2str(sat.a) );
    set(handles.edit_m,'String', num2str(sat.m) );
    set(handles.edit_sigma_nf,'String', num2str(sat.sigma_nf) );
    set(handles.edit_epsilon_nf,'String', num2str(sat.epsilon_nf) );
    set(handles.edit_mu_nf,'String', num2str(sat.mu_nf) );
    set(handles.edit_epsilon_f1,'String', num2str(sat.epsilon_f1) );
    set(handles.edit_sigma_f1,'String', num2str(sat.sigma_f1) );
    set(handles.edit_frequency,'String', num2str(sat.omega/(2*pi*1e6)) );
    set(handles.edit_S0,'String', num2str(sat.S0) );

    h.S = sat.S;
    h.ivel1 = sat.ivel1;
    h.ivel2 = sat.ivel2;
    h.iatt1 = sat.iatt1;
    h.iatt2 = sat.iatt2;

    h.epsilon_m = sat.epsilon_m;
    h.a = sat.a;
    h.m = sat.m;
    h.sigma_nf = sat.sigma_nf;
    h.epsilon_nf = sat.epsilon_nf;
    h.mu_nf = sat.mu_nf;
    h.epsilon_f1 = sat.epsilon_f1;
    h.sigma_f1 = sat.sigma_f1;
    h.omega = sat.omega;
    h.S0 = sat.S0;

    setappdata(handles.fig_nanofluid, 'h', h)
end

% --------------------------------------------------------------------
function exit_menu_item_Callback(hObject, eventdata, handles)
% hObject    handle to exit_menu_item (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
delete(handles.fig_nanofluid)


% --- Executes on selection change in popupmenu2.
function popupmenu2_Callback(hObject, eventdata, handles)
% hObject    handle to popupmenu2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns popupmenu2 contents as cell array
%        contents{get(hObject,'Value')} returns selected item from popupmenu2


% --- Executes during object creation, after setting all properties.
function popupmenu2_CreateFcn(hObject, eventdata, handles)
% hObject    handle to popupmenu2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in popupmenu3.
function popupmenu3_Callback(hObject, eventdata, handles)
% hObject    handle to popupmenu3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns popupmenu3 contents as cell array
%        contents{get(hObject,'Value')} returns selected item from popupmenu3


% --- Executes during object creation, after setting all properties.
function popupmenu3_CreateFcn(hObject, eventdata, handles)
% hObject    handle to popupmenu3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit_epsilon_nf_Callback(hObject, eventdata, handles)
% hObject    handle to edit_epsilon_nf (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_epsilon_nf as text
%        str2double(get(hObject,'String')) returns contents of edit_epsilon_nf as a double


% --- Executes during object creation, after setting all properties.
function edit_epsilon_nf_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_epsilon_nf (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit_mu_nf_Callback(hObject, eventdata, handles)
% hObject    handle to edit_mu_nf (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_mu_nf as text
%        str2double(get(hObject,'String')) returns contents of edit_mu_nf as a double


% --- Executes during object creation, after setting all properties.
function edit_mu_nf_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_mu_nf (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit_sigma_nf_Callback(hObject, eventdata, handles)
% hObject    handle to edit_sigma_nf (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_sigma_nf as text
%        str2double(get(hObject,'String')) returns contents of edit_sigma_nf as a double


% --- Executes during object creation, after setting all properties.
function edit_sigma_nf_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_sigma_nf (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit_epsilon_m_Callback(hObject, eventdata, handles)
% hObject    handle to edit_epsilon_m (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_epsilon_m as text
%        str2double(get(hObject,'String')) returns contents of edit_epsilon_m as a double


% --- Executes during object creation, after setting all properties.
function edit_epsilon_m_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_epsilon_m (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit_a_Callback(hObject, eventdata, handles)
% hObject    handle to edit_a (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_a as text
%        str2double(get(hObject,'String')) returns contents of edit_a as a double


% --- Executes during object creation, after setting all properties.
function edit_a_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_a (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit_m_Callback(hObject, eventdata, handles)
% hObject    handle to edit_m (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_m as text
%        str2double(get(hObject,'String')) returns contents of edit_m as a double


% --- Executes during object creation, after setting all properties.
function edit_m_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_m (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit_frequency_Callback(hObject, eventdata, handles)
% hObject    handle to edit_frequency (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_frequency as text
%        str2double(get(hObject,'String')) returns contents of edit_frequency as a double


% --- Executes during object creation, after setting all properties.
function edit_frequency_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_frequency (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit_S0_Callback(hObject, eventdata, handles)
% hObject    handle to edit_S0 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_S0 as text
%        str2double(get(hObject,'String')) returns contents of edit_S0 as a double


% --- Executes during object creation, after setting all properties.
function edit_S0_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_S0 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit_epsilon_f1_Callback(hObject, eventdata, handles)
% hObject    handle to edit_epsilon_f1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_epsilon_f1 as text
%        str2double(get(hObject,'String')) returns contents of edit_epsilon_f1 as a double


% --- Executes during object creation, after setting all properties.
function edit_epsilon_f1_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_epsilon_f1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit_sigma_f1_Callback(hObject, eventdata, handles)
% hObject    handle to edit_sigma_f1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_sigma_f1 as text
%        str2double(get(hObject,'String')) returns contents of edit_sigma_f1 as a double


% --- Executes during object creation, after setting all properties.
function edit_sigma_f1_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_sigma_f1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in popupmenu_repeat_att.
function popupmenu_repeat_att_Callback(hObject, eventdata, handles)
% hObject    handle to popupmenu_repeat_att (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns popupmenu_repeat_att contents as cell array
%        contents{get(hObject,'Value')} returns selected item from popupmenu_repeat_att


% --- Executes during object creation, after setting all properties.
function popupmenu_repeat_att_CreateFcn(hObject, eventdata, handles)
% hObject    handle to popupmenu_repeat_att (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in popupmenu_repeat_vel.
function popupmenu_repeat_vel_Callback(hObject, eventdata, handles)
% hObject    handle to popupmenu_repeat_vel (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns popupmenu_repeat_vel contents as cell array
%        contents{get(hObject,'Value')} returns selected item from popupmenu_repeat_vel


% --- Executes during object creation, after setting all properties.
function popupmenu_repeat_vel_CreateFcn(hObject, eventdata, handles)
% hObject    handle to popupmenu_repeat_vel (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in popupmenu_baseline_att.
function popupmenu_baseline_att_Callback(hObject, eventdata, handles)
% hObject    handle to popupmenu_baseline_att (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns popupmenu_baseline_att contents as cell array
%        contents{get(hObject,'Value')} returns selected item from popupmenu_baseline_att


% --- Executes during object creation, after setting all properties.
function popupmenu_baseline_att_CreateFcn(hObject, eventdata, handles)
% hObject    handle to popupmenu_baseline_att (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in popupmenu_baseline_vel.
function popupmenu_baseline_vel_Callback(hObject, eventdata, handles)
% hObject    handle to popupmenu_baseline_vel (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns popupmenu_baseline_vel contents as cell array
%        contents{get(hObject,'Value')} returns selected item from popupmenu_baseline_vel


% --- Executes during object creation, after setting all properties.
function popupmenu_baseline_vel_CreateFcn(hObject, eventdata, handles)
% hObject    handle to popupmenu_baseline_vel (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


function update_fig(handles)
h = getappdata(handles.fig_nanofluid, 'h');
str = get_str_locale();

gridx = 0.5*(h.panels(h.no_panel).grid.grx(1:end-1)+h.panels(h.no_panel).grid.grx(2:end));
gridz = 0.5*(h.panels(h.no_panel).grid.grz(1:end-1)+h.panels(h.no_panel).grid.grz(2:end));
n=length(gridz);
m=length(gridx);

imagesc(gridx,gridz,reshape(h.S,n,m),'Parent',handles.axes1)
set(handles.axes1,'DataAspectRatio',[1 1 1],'YDir','normal')
title(str.s308,'FontSize',14)
xlabel(str.s119)
ylabel(str.s120)
caxis(handles.axes1,[0 1])
colorbar('peer', handles.axes1)


% --------------------------------------------------------------------
function Save_menu_item_Callback(hObject, eventdata, handles)
% hObject    handle to Save_menu_item (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

h = getappdata(handles.fig_nanofluid, 'h');
if isfield(h.panels(h.no_panel),'saturation')
    n = length(h.panels(h.no_panel).saturation)+1;
else
    n = 1;
end

nameDefault = ['SaturationMap',num2str(n)];  
prompt = {'Saturation map name:                                                               '};
title = 'Save saturation results';
nblines = 1;
answer = inputdlg(prompt,title,nblines,{nameDefault});
name = [];
if ~isempty(answer)
    name=answer{1};
else
    return;
end

if isfield(h.panels(h.no_panel),'saturation')
	flag=0;
	for n=1:length(h.panels(h.no_panel).saturation)
		if strcmp(h.panels(h.no_panel).saturation(n).name, name)
			no_saturation = n;
			flag = 1;
			break;
		end
	end
	if flag==1
		str = get_str_locale();
		rep=questdlg([str.s187,' ',name,'?']);
		if ~strcmp(rep,'Yes')
			return
		end
	else
		no_saturation = 1+length(h.panels(h.no_panel).saturation);
	end
else
	no_saturation = 1;
end
h.panels(h.no_panel).saturation(no_saturation).name = name;
h.panels(h.no_panel).saturation(no_saturation).S = h.S;
h.panels(h.no_panel).saturation(no_saturation).ivel1 = h.ivel1;
h.panels(h.no_panel).saturation(no_saturation).ivel2 = h.ivel2;
h.panels(h.no_panel).saturation(no_saturation).iatt1 = h.iatt1;
h.panels(h.no_panel).saturation(no_saturation).iatt2 = h.iatt2;
h.panels(h.no_panel).saturation(no_saturation).epsilon_m = h.epsilon_m;
h.panels(h.no_panel).saturation(no_saturation).a = h.a;
h.panels(h.no_panel).saturation(no_saturation).m = h.m;
h.panels(h.no_panel).saturation(no_saturation).sigma_nf = h.sigma_nf;
h.panels(h.no_panel).saturation(no_saturation).epsilon_nf = h.epsilon_nf;
h.panels(h.no_panel).saturation(no_saturation).mu_nf = h.mu_nf;
h.panels(h.no_panel).saturation(no_saturation).epsilon_f1 = h.epsilon_f1;
h.panels(h.no_panel).saturation(no_saturation).sigma_f1 = h.sigma_f1;
h.panels(h.no_panel).saturation(no_saturation).omega = h.omega;
h.panels(h.no_panel).saturation(no_saturation).S0 = h.S0;

load(h.db_file,'panels')
panels(h.no_panel).saturation = h.panels(h.no_panel).saturation; %#ok<NASGU>
save(h.db_file,'panels','-append')
setappdata(handles.fig_nanofluid, 'h', h)

for n = 1:length(panels(h.no_panel).saturation)
    name_maps{n} = char( panels(h.no_panel).saturation(n).name );
end
set(handles.popupmenu_Smaps,'String',name_maps);
set(handles.popupmenu_Smaps,'Value',no_saturation);

% --- Executes on selection change in popupmenu_Smaps.
function popupmenu_Smaps_Callback(hObject, eventdata, handles)
% hObject    handle to popupmenu_Smaps (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

no = get(hObject,'Value');
h = getappdata(handles.fig_nanofluid, 'h');

sat = h.panels(h.no_panel).saturation(no);
    
set(handles.popupmenu_baseline_vel,'Value',sat.ivel1);
set(handles.popupmenu_baseline_att,'Value',sat.iatt1);
set(handles.popupmenu_repeat_vel,'Value',sat.ivel2);
set(handles.popupmenu_repeat_att,'Value',sat.iatt2);

set(handles.edit_epsilon_m,'String', num2str(sat.epsilon_m ) );
set(handles.edit_a,'String', num2str(sat.a) );
set(handles.edit_m,'String', num2str(sat.m) );
set(handles.edit_sigma_nf,'String', num2str(sat.sigma_nf) );
set(handles.edit_epsilon_nf,'String', num2str(sat.epsilon_nf) );
set(handles.edit_mu_nf,'String', num2str(sat.mu_nf) );
set(handles.edit_epsilon_f1,'String', num2str(sat.epsilon_f1) );
set(handles.edit_sigma_f1,'String', num2str(sat.sigma_f1) );
set(handles.edit_frequency,'String', num2str(sat.omega/(2*pi*1e6)) );
set(handles.edit_S0,'String', num2str(sat.S0) );

h.S = sat.S;
h.ivel1 = sat.ivel1;
h.ivel2 = sat.ivel2;
h.iatt1 = sat.iatt1;
h.iatt2 = sat.iatt2;

h.epsilon_m = sat.epsilon_m;
h.a = sat.a;
h.m = sat.m;
h.sigma_nf = sat.sigma_nf;
h.epsilon_nf = sat.epsilon_nf;
h.mu_nf = sat.mu_nf;
h.epsilon_f1 = sat.epsilon_f1;
h.sigma_f1 = sat.sigma_f1;
h.omega = sat.omega;
h.S0 = sat.S0;

setappdata(handles.fig_nanofluid, 'h', h)
update_fig(handles);


% --- Executes during object creation, after setting all properties.
function popupmenu_Smaps_CreateFcn(hObject, eventdata, handles)
% hObject    handle to popupmenu_Smaps (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --------------------------------------------------------------------
function Results_Menu_Callback(hObject, eventdata, handles)
% hObject    handle to Results_Menu (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --------------------------------------------------------------------
function Figure_menu_item_Callback(hObject, eventdata, handles)
% hObject    handle to Figure_menu_item (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

figure
h = getappdata(handles.fig_nanofluid, 'h');
str = get_str_locale();

gridx = 0.5*(h.panels(h.no_panel).grid.grx(1:end-1)+h.panels(h.no_panel).grid.grx(2:end));
gridz = 0.5*(h.panels(h.no_panel).grid.grz(1:end-1)+h.panels(h.no_panel).grid.grz(2:end));
n=length(gridz);
m=length(gridx);

imagesc(gridx,gridz,reshape(h.S,n,m))
set(gca,'DataAspectRatio',[1 1 1],'YDir','normal')
title(str.s308,'FontSize',14)
xlabel(str.s119)
ylabel(str.s120)
caxis([0 1])
colorbar
