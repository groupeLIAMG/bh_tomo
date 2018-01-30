function varargout = bh_tomo_tltomo2d(varargin)
% BH_TOMO_TLTOMO2D MATLAB code for bh_tomo_tltomo2d.fig
%      BH_TOMO_TLTOMO2D, by itself, creates a new BH_TOMO_TLTOMO2D or raises the existing
%      singleton*.
%
%      H = BH_TOMO_TLTOMO2D returns the handle to a new BH_TOMO_TLTOMO2D or the handle to
%      the existing singleton*.
%
%      BH_TOMO_TLTOMO2D('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in BH_TOMO_TLTOMO2D.M with the given input arguments.
%
%      BH_TOMO_TLTOMO2D('Property','Value',...) creates a new BH_TOMO_TLTOMO2D or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before bh_tomo_tltomo2d_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to bh_tomo_tltomo2d_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help bh_tomo_tltomo2d

% Last Modified by GUIDE v2.5 04-Jun-2015 21:23:53

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @bh_tomo_tltomo2d_OpeningFcn, ...
                   'gui_OutputFcn',  @bh_tomo_tltomo2d_OutputFcn, ...
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


% --- Executes just before bh_tomo_tltomo2d is made visible.
function bh_tomo_tltomo2d_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to bh_tomo_tltomo2d (see VARARGIN)

% Choose default command line output for bh_tomo_tltomo2d
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes bh_tomo_tltomo2d wait for user response (see UIRESUME)
% uiwait(handles.figure_tltomo2d);
h.db_file = get(handles.figure_tltomo2d,'UserData');
h.slo_cmin = 0.06;
h.slo_cmax = 0.12;
h.amp_cmin = 0;
h.amp_cmax = 3;
str = get_str_locale();

setappdata(handles.figure_tltomo2d, 'h', h)
setappdata(handles.figure_tltomo2d, 'str', str)

% --- Outputs from this function are returned to the command line.
function varargout = bh_tomo_tltomo2d_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


% --- Executes on selection change in listbox_baseline.
function listbox_baseline_Callback(hObject, eventdata, handles)
% hObject    handle to listbox_baseline (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns listbox_baseline contents as cell array
%        contents{get(hObject,'Value')} returns selected item from listbox_baseline


% --- Executes during object creation, after setting all properties.
function listbox_baseline_CreateFcn(hObject, eventdata, handles)
% hObject    handle to listbox_baseline (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: listbox controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in listbox_repeat.
function listbox_repeat_Callback(hObject, eventdata, handles)
% hObject    handle to listbox_repeat (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns listbox_repeat contents as cell array
%        contents{get(hObject,'Value')} returns selected item from listbox_repeat


% --- Executes during object creation, after setting all properties.
function listbox_repeat_CreateFcn(hObject, eventdata, handles)
% hObject    handle to listbox_repeat (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: listbox controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --------------------------------------------------------------------
function MenuFile_Callback(hObject, eventdata, handles)
% hObject    handle to MenuFile (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --------------------------------------------------------------------
function MenuQuit_Callback(hObject, eventdata, handles)
% hObject    handle to MenuQuit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
delete(handles.figure_tltomo2d)


% --------------------------------------------------------------------
function MenuChoosePanel_Callback(hObject, eventdata, handles)
% hObject    handle to MenuChoosePanel (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
h = getappdata(handles.figure_tltomo2d, 'h');
str = getappdata(handles.figure_tltomo2d, 'str');
[h.no_panel,h.db_file,hh] = choisirPanneau('UserData', h.db_file);
if ishandle(hh), delete(hh), end
if h.no_panel==0
	return
end
load(h.db_file,'panels')
h.panel = panels(h.no_panel);
if ~isfield(h.panel, 'grid')
	errordlg(str.s150, str.s45,'model')
	return
elseif isempty(h.panel.grid)
	errordlg(str.s150, str.s45,'model')
	return
end
set(handles.text_panel,'String',[str.s151,':   ', h.panel.name])

load(h.db_file,'mogs')
mog_date = cell(1,length(h.panel.mogs));
for n = 1:length(h.panel.mogs)
    mog_date{n} = [mogs(h.panel.mogs(n)).name, ',', ...
        mogs(h.panel.mogs(n)).date];
end
set(handles.listbox_baseline,'String',mog_date)
set(handles.listbox_repeat,'String',mog_date)
set(handles.listbox_baseline,'Value',1)
set(handles.listbox_repeat,'Value',1)

inv_name = cell(1,length(h.panel.inv_res));
for n=1:length(h.panel.inv_res)
    inv_name{n} = h.panel.inv_res(n).name;
end
set(handles.popupmenu_ref,'String',inv_name)
set(handles.popupmenu_di_rays,'String',inv_name)
set(handles.popupmenu_si_rays,'String',inv_name)

setappdata(handles.figure_tltomo2d, 'h', h)

% --------------------------------------------------------------------
function MenuSave_Callback(hObject, eventdata, handles)
% hObject    handle to MenuSave (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
h = getappdata(handles.figure_tltomo2d, 'h');
if isfield(h.panel,'inv_res')
    n = length(h.panel.inv_res)+1;
else
    n = 1;
end

ind = get(handles.popupmenu_type_data,'Value');
switch ind
    case 1
        dType = '-tt';
    case 2
        dType = '-amp';
    case 3
        dType = '-fce';
end

if get(handles.popupmenu_type_inversion,'Value')==1
    iType = '-diff_inv';
else
    iType = '-simult_inv';
end

nameDefault = ['tomo',num2str(n) iType dType];  
prompt = {'Inversion name:                                                               '};
title = 'Save inversion results';
nblines = 1;
answer = inputdlg(prompt,title,nblines,{nameDefault});
if ~isempty(answer)
    name=answer{1};
else
    return;
end

if isempty(name)
	str = getappdata(handles.figure_tltomo2d, 'str');
	uiwait(errordlg(str.s185,str.s45))
	return
end

if isfield(h.panel,'inv_res')
	flag=0;
	for n=1:length(h.panel.inv_res)
		if strcmp(h.panel.inv_res(n).name, name)
			no_inv_res = n;
			flag = 1;
			break;
		end
	end
	if flag==1
		str = getappdata(handles.figure_tltomo2d, 'str');
		rep=questdlg([str.s187,' ',name,'?']);
		if ~strcmp(rep,'Yes')
			return
		end
	else
		no_inv_res = 1+length(h.panel.inv_res);
	end
else
	no_inv_res = 1;
end
h.panel.inv_res(no_inv_res).name = name;
h.panel.inv_res(no_inv_res).tomo = h.tomo;
h.panel.inv_res(no_inv_res).param = h.param;

load(h.db_file,'panels')
panels(h.no_panel).inv_res = h.panel.inv_res; %#ok<NASGU>
save(h.db_file,'panels','-append')
setappdata(handles.figure_tltomo2d, 'h', h)

inv_name = cell(1,length(h.panel.inv_res));
for n=1:length(h.panel.inv_res)
    inv_name{n} = h.panel.inv_res(n).name;
end
set(handles.popupmenu_ref,'String',inv_name)
set(handles.popupmenu_di_rays,'String',inv_name)
set(handles.popupmenu_si_rays,'String',inv_name)


% --- Executes on selection change in popupmenu_type_data.
function popupmenu_type_data_Callback(hObject, eventdata, handles)
% hObject    handle to popupmenu_type_data (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns popupmenu_type_data contents as cell array
%        contents{get(hObject,'Value')} returns selected item from popupmenu_type_data

h = getappdata(handles.figure_tltomo2d, 'h');
if get(hObject,'Value')==1
    set(handles.popupmenu_di_rays,'Visible','off')
    set(handles.popupmenu_si_rays,'Visible','off')
    set(handles.text_di_rays,'Visible','off')
    set(handles.text_si_rays,'Visible','off')
    set(handles.edit_cmin,'String',num2str(h.slo_cmin))
    set(handles.edit_cmax,'String',num2str(h.slo_cmax))
else
    set(handles.popupmenu_di_rays,'Visible','on')
    set(handles.popupmenu_si_rays,'Visible','on')
    set(handles.text_di_rays,'Visible','on')
    set(handles.text_si_rays,'Visible','on')
    set(handles.edit_cmin,'String',num2str(h.amp_cmin))
    set(handles.edit_cmax,'String',num2str(h.amp_cmax))
end

% --- Executes during object creation, after setting all properties.
function popupmenu_type_data_CreateFcn(hObject, eventdata, handles)
% hObject    handle to popupmenu_type_data (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in popupmenu_type_inversion.
function popupmenu_type_inversion_Callback(hObject, eventdata, handles)
% hObject    handle to popupmenu_type_inversion (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns popupmenu_type_inversion contents as cell array
%        contents{get(hObject,'Value')} returns selected item from popupmenu_type_inversion
if get(hObject,'Value')==1
    set(handles.uipanel_diff_inv_params,'Visible','on')
    set(handles.uipanel_simult_inv_params,'Visible','off')
elseif get(hObject,'Value')==2
    set(handles.uipanel_diff_inv_params,'Visible','off')
    set(handles.uipanel_simult_inv_params,'Visible','on')    
end

% --- Executes during object creation, after setting all properties.
function popupmenu_type_inversion_CreateFcn(hObject, eventdata, handles)
% hObject    handle to popupmenu_type_inversion (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in popupmenu_ref.
function popupmenu_ref_Callback(hObject, eventdata, handles)
% hObject    handle to popupmenu_ref (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns popupmenu_ref contents as cell array
%        contents{get(hObject,'Value')} returns selected item from popupmenu_ref


% --- Executes during object creation, after setting all properties.
function popupmenu_ref_CreateFcn(hObject, eventdata, handles)
% hObject    handle to popupmenu_ref (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in pushbutton_inv.
function pushbutton_inv_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton_inv (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

h = getappdata(handles.figure_tltomo2d, 'h');
str = getappdata(handles.figure_tltomo2d, 'str');

set(handles.text_message_inv, 'String', 'Starting ...')

p.db_file = h.db_file;
p.saveInvData = 1;
if isfield(h.panel.grid.cont,'ind_reservoir') && ...
        get(handles.checkbox_wt,'Value')==1
    p.ind_reservoir = h.panel.grid.cont.ind_reservoir;
else
    p.ind_reservoir = false(length(h.panel.grid.grz)-1,...
        length(h.panel.grid.grx)-1);
end
p.weight_reservoir = str2double(get(handles.edit_wt,'String'));
p.max_it = str2double(get(handles.edit_numit,'String'));
p.ref_inv_no = get(handles.popupmenu_ref,'Value');

if get(handles.popupmenu_type_data,'Value')==1
    p.type_data = 'tt';
    p.tomo_amp = 0;
    
    if h.panel.inv_res(p.ref_inv_no).param.tomo_amp==1
        errordlg('Cannot use attenuation model as reference for traveltime inversion')
        set(handles.text_message_inv, 'String', '')
        return
    end
    
elseif get(handles.popupmenu_type_data,'Value')==2
    p.type_data = 'amp';
    p.tomo_amp = 1;

    if h.panel.inv_res(p.ref_inv_no).param.tomo_amp==0
        errordlg('Cannot use velocity model as reference for amplitude inversion')
        set(handles.text_message_inv, 'String', '')
        return
    end
    
else
    p.type_data = 'fce';
    p.tomo_amp = 1;

    if h.panel.inv_res(p.ref_inv_no).param.tomo_amp==0
        errordlg('Cannot use velocity model as reference for amplitude inversion')
        set(handles.text_message_inv, 'String', '')
        return
    end
    
end

if get(handles.popupmenu_type_data,'Value')==1
    clim = [h.slo_cmin h.slo_cmax];
else
    clim = [h.amp_cmin h.amp_cmax];
end
cla(handles.axes1);cla(handles.axes2)
cbh = findobj( 0, 'tag', 'Colorbar' );   %%%YH
for i = 1: length(cbh) 
    colorbar(cbh(i),'off') 
end

maps = get(handles.popupmenu_colormap,'String');
cmap = maps{get(handles.popupmenu_colormap,'Value')};

g_handles = {clim; cmap; handles.axes1; handles.axes2; ''};

if get(handles.popupmenu_type_inversion,'Value')==1
    
    % difference inversion
    p.beta = str2double(get(handles.edit_di_beta,'String'));
    p.lambda = str2double(get(handles.edit_di_lambda,'String'));
    p.mu = str2double(get(handles.edit_di_mu,'String'));
    p.eta = str2double(get(handles.edit_di_eta,'String'));
    p.damping = str2double(get(handles.edit_di_damping,'String'));

    p.mog_no = get(handles.listbox_repeat,'Value');
    p.L_tomo_no = get(handles.popupmenu_di_rays,'Value');

    M=length(h.panel.inv_res(p.ref_inv_no).tomo.z);
    N=length(h.panel.inv_res(p.ref_inv_no).tomo.x);
    if p.tomo_amp==0
        imagesc(h.panel.inv_res(p.ref_inv_no).tomo.x,...
            h.panel.inv_res(p.ref_inv_no).tomo.z,...
            1./reshape(h.panel.inv_res(p.ref_inv_no).tomo.s,M,N),...
            'Parent',handles.axes1)
        set(handles.axes1,'DataAspectRatio',[1 1 1],'YDir','normal')
        caxis(handles.axes1, clim)
        title(handles.axes1,'Baseline survey: velocity','FontSize',14)
    else
        imagesc(h.panel.inv_res(p.ref_inv_no).tomo.x,...
            h.panel.inv_res(p.ref_inv_no).tomo.z,...
            reshape(h.panel.inv_res(p.ref_inv_no).tomo.s,M,N),...
            'Parent',handles.axes1)
        set(handles.axes1,'DataAspectRatio',[1 1 1],'YDir','normal')
        caxis(handles.axes1, clim)
        title(handles.axes1,'Baseline survey: attenuation','FontSize',14)
    end
    xlabel(handles.axes1,str.s119)
    ylabel(handles.axes1,str.s120)
    if ~isempty(clim), caxis(handles.axes1, clim), end
    if get(handles.checkbox_colorbar,'Value')==1
        colorbar('peer',handles.axes1)
    end
    colormap(cmap)
    drawnow

    h.tomo = inv_difference2D(h.panel, p, handles.text_message_inv, g_handles);
else
    % simultaneous inversion
    p.alpha = str2double(get(handles.edit_si_alpha,'String'));
    p.beta = str2double(get(handles.edit_si_beta,'String'));
    p.lambda = str2double(get(handles.edit_si_lambda,'String'));
    p.mu = str2double(get(handles.edit_si_mu,'String'));
    p.eta = str2double(get(handles.edit_si_eta,'String'));
    p.damping = str2double(get(handles.edit_si_damping,'String'));

    p.mog_no0 = get(handles.listbox_baseline,'Value');
    p.mog_no1 = get(handles.listbox_repeat,'Value');
    p.L_tomo_no = get(handles.popupmenu_si_rays,'Value');
    
    if p.tomo_amp == 1 && ~isfield(h.panel.inv_res(p.L_tomo_no).tomo,'no_trace0')
        errordlg('Selected rays not obtained with simultaneous inversion')
        set(handles.text_message_inv, 'String', '')
        return
    end

    h.tomo = inv_simultaneous2D(h.panel, p, handles.text_message_inv, g_handles);
end

h.param = p;
set(handles.text_message_inv, 'String', '')
setappdata(handles.figure_tltomo2d, 'h', h)


% --- Executes on button press in checkbox_wt.
function checkbox_wt_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox_wt (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkbox_wt



function edit_wt_Callback(hObject, eventdata, handles)
% hObject    handle to edit_wt (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_wt as text
%        str2double(get(hObject,'String')) returns contents of edit_wt as a double


% --- Executes during object creation, after setting all properties.
function edit_wt_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_wt (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit_numit_Callback(hObject, eventdata, handles)
% hObject    handle to edit_numit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_numit as text
%        str2double(get(hObject,'String')) returns contents of edit_numit as a double


% --- Executes during object creation, after setting all properties.
function edit_numit_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_numit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit_di_beta_Callback(hObject, eventdata, handles)
% hObject    handle to edit_di_beta (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_di_beta as text
%        str2double(get(hObject,'String')) returns contents of edit_di_beta as a double


% --- Executes during object creation, after setting all properties.
function edit_di_beta_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_di_beta (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit_di_lambda_Callback(hObject, eventdata, handles)
% hObject    handle to edit_di_lambda (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_di_lambda as text
%        str2double(get(hObject,'String')) returns contents of edit_di_lambda as a double


% --- Executes during object creation, after setting all properties.
function edit_di_lambda_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_di_lambda (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit_di_mu_Callback(hObject, eventdata, handles)
% hObject    handle to edit_di_mu (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_di_mu as text
%        str2double(get(hObject,'String')) returns contents of edit_di_mu as a double


% --- Executes during object creation, after setting all properties.
function edit_di_mu_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_di_mu (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit_di_eta_Callback(hObject, eventdata, handles)
% hObject    handle to edit_di_eta (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_di_eta as text
%        str2double(get(hObject,'String')) returns contents of edit_di_eta as a double


% --- Executes during object creation, after setting all properties.
function edit_di_eta_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_di_eta (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit_di_damping_Callback(hObject, eventdata, handles)
% hObject    handle to edit_di_damping (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_di_damping as text
%        str2double(get(hObject,'String')) returns contents of edit_di_damping as a double


% --- Executes during object creation, after setting all properties.
function edit_di_damping_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_di_damping (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit_si_damping_Callback(hObject, eventdata, handles)
% hObject    handle to edit_si_damping (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_si_damping as text
%        str2double(get(hObject,'String')) returns contents of edit_si_damping as a double


% --- Executes during object creation, after setting all properties.
function edit_si_damping_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_si_damping (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit_si_eta_Callback(hObject, eventdata, handles)
% hObject    handle to edit_si_eta (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_si_eta as text
%        str2double(get(hObject,'String')) returns contents of edit_si_eta as a double


% --- Executes during object creation, after setting all properties.
function edit_si_eta_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_si_eta (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit_si_mu_Callback(hObject, eventdata, handles)
% hObject    handle to edit_si_mu (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_si_mu as text
%        str2double(get(hObject,'String')) returns contents of edit_si_mu as a double


% --- Executes during object creation, after setting all properties.
function edit_si_mu_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_si_mu (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit_si_lambda_Callback(hObject, eventdata, handles)
% hObject    handle to edit_si_lambda (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_si_lambda as text
%        str2double(get(hObject,'String')) returns contents of edit_si_lambda as a double


% --- Executes during object creation, after setting all properties.
function edit_si_lambda_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_si_lambda (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit_si_beta_Callback(hObject, eventdata, handles)
% hObject    handle to edit_si_beta (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_si_beta as text
%        str2double(get(hObject,'String')) returns contents of edit_si_beta as a double


% --- Executes during object creation, after setting all properties.
function edit_si_beta_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_si_beta (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit_si_alpha_Callback(hObject, eventdata, handles)
% hObject    handle to edit_si_alpha (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_si_alpha as text
%        str2double(get(hObject,'String')) returns contents of edit_si_alpha as a double


% --- Executes during object creation, after setting all properties.
function edit_si_alpha_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_si_alpha (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in popupmenu_di_rays.
function popupmenu_di_rays_Callback(hObject, eventdata, handles)
% hObject    handle to popupmenu_di_rays (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns popupmenu_di_rays contents as cell array
%        contents{get(hObject,'Value')} returns selected item from popupmenu_di_rays


% --- Executes during object creation, after setting all properties.
function popupmenu_di_rays_CreateFcn(hObject, eventdata, handles)
% hObject    handle to popupmenu_di_rays (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in popupmenu_si_rays.
function popupmenu_si_rays_Callback(hObject, eventdata, handles)
% hObject    handle to popupmenu_si_rays (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns popupmenu_si_rays contents as cell array
%        contents{get(hObject,'Value')} returns selected item from popupmenu_si_rays


% --- Executes during object creation, after setting all properties.
function popupmenu_si_rays_CreateFcn(hObject, eventdata, handles)
% hObject    handle to popupmenu_si_rays (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in popupmenu_si_rays1.
function popupmenu_si_rays1_Callback(hObject, eventdata, handles)
% hObject    handle to popupmenu_si_rays1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns popupmenu_si_rays1 contents as cell array
%        contents{get(hObject,'Value')} returns selected item from popupmenu_si_rays1


% --- Executes during object creation, after setting all properties.
function popupmenu_si_rays1_CreateFcn(hObject, eventdata, handles)
% hObject    handle to popupmenu_si_rays1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in checkbox_colorbar.
function checkbox_colorbar_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox_colorbar (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkbox_colorbar
updateFigures(handles)



function edit_cmin_Callback(hObject, eventdata, handles)
% hObject    handle to edit_cmin (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_cmin as text
%        str2double(get(hObject,'String')) returns contents of edit_cmin as a double
h = getappdata(handles.figure_tltomo2d, 'h');
if get(handles.popupmenu_type_data,'Value')==1
    h.slo_cmin = str2double(get(hObject,'String'));
else
    h.amp_cmin = str2double(get(hObject,'String'));
end
setappdata(handles.figure_tltomo2d, 'h', h)
updateFigures(handles)

% --- Executes during object creation, after setting all properties.
function edit_cmin_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_cmin (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit_cmax_Callback(hObject, eventdata, handles)
% hObject    handle to edit_cmax (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_cmax as text
%        str2double(get(hObject,'String')) returns contents of edit_cmax as a double
h = getappdata(handles.figure_tltomo2d, 'h');
if get(handles.popupmenu_type_data,'Value')==1
    h.slo_cmax = str2double(get(hObject,'String'));
else
    h.amp_cmax = str2double(get(hObject,'String'));
end
setappdata(handles.figure_tltomo2d, 'h', h)
updateFigures(handles)

% --- Executes during object creation, after setting all properties.
function edit_cmax_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_cmax (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in popupmenu_colormap.
function popupmenu_colormap_Callback(hObject, eventdata, handles)
% hObject    handle to popupmenu_colormap (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns popupmenu_colormap contents as cell array
%        contents{get(hObject,'Value')} returns selected item from popupmenu_colormap
updateFigures(handles)

% --- Executes during object creation, after setting all properties.
function popupmenu_colormap_CreateFcn(hObject, eventdata, handles)
% hObject    handle to popupmenu_colormap (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function updateFigures(handles)

maps = cellstr(get(handles.popupmenu_colormap,'String'));
cmap=maps{get(handles.popupmenu_colormap,'Value')};
colormap(handles.axes1, cmap)
colormap(handles.axes2, cmap)

clim=[str2double(get(handles.edit_cmin,'String')) ...
    str2double(get(handles.edit_cmax,'String'))];
caxis(handles.axes1, clim);
caxis(handles.axes2, clim);

if get(handles.checkbox_colorbar,'Value')==1
    axes(handles.axes1);
    colorbar
    axes(handles.axes2);
    colorbar
else
    axes(handles.axes1);
    colorbar off
    axes(handles.axes2);
    colorbar off
end


% --- Executes on button press in pushbutton_cla.
function pushbutton_cla_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton_cla (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
cla(handles.axes1,'reset')
cla(handles.axes2,'reset')
