function varargout = bh_tomo_interp(varargin)
% BH_TOMO_INTERP M-file for bh_tomo_interp.fig
%      BH_TOMO_INTERP, by itself, creates a new BH_TOMO_INTERP or raises the existing
%      singleton*.
%
%      H = BH_TOMO_INTERP returns the handle to a new BH_TOMO_INTERP or the handle to
%      the existing singleton*.
%
%      BH_TOMO_INTERP('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in BH_TOMO_INTERP.M with the given input arguments.
%
%      BH_TOMO_INTERP('Property','Value',...) creates a new BH_TOMO_INTERP or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before bh_tomo_interp_OpeningFunction gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to bh_tomo_interp_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help bh_tomo_interp

% Last Modified by GUIDE v2.5 17-Dec-2012 13:47:26

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @bh_tomo_interp_OpeningFcn, ...
                   'gui_OutputFcn',  @bh_tomo_interp_OutputFcn, ...
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


% --- Executes just before bh_tomo_interp is made visible.
function bh_tomo_interp_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to bh_tomo_interp (see VARARGIN)

% Choose default command line output for bh_tomo_interp
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes bh_tomo_interp wait for user response (see UIRESUME)
% uiwait(handles.fig_interp);
h.db_file = get(handles.fig_interp,'UserData');
h.panel = [];
h.cb = '';
h.str = get_str_locale();

setappdata(handles.fig_interp, 'h', h)

set_string_locale(handles)


% --- Outputs from this function are returned to the command line.
function varargout = bh_tomo_interp_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


function popupmenu_tomo_lent_Callback(hObject, eventdata, handles)
update_plots(handles)

function popupmenu_tomo_lent_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


function popupmenu_tomo_att_Callback(hObject, eventdata, handles)
update_plots(handles)

function popupmenu_tomo_att_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


function popupmenu_prop_Callback(hObject, eventdata, handles)
update_plots(handles)


function popupmenu_prop_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


function FileMenu_Callback(hObject, eventdata, handles)


function OpenMenuItem_Callback(hObject, eventdata, handles)
h = getappdata(handles.fig_interp, 'h');
[h.no_panel,h.db_file,hh] = choisirPanneau('UserData', h.db_file);
if ishandle(hh), delete(hh), end
if h.no_panel==0
	return
end
load(h.db_file,'panels')
h.panel = panels(h.no_panel);
if ~isfield(h.panel, 'grid')
    errordlg(h.str.s150, h.str.s45,'model')
    return
elseif isempty(h.panel.grid)
    errordlg(h.str.s150, h.str.s45,'model')
    return
end

h.slo_ind = [];
h.att_ind = [];
for n=1:length(h.panel.inv_res)
	if h.panel.inv_res(n).param.tomo_amp == 0
		h.slo_ind = [h.slo_ind n];
	else
		h.att_ind = [h.att_ind n];
	end
end
slo_name = cell(length(h.slo_ind),1);
att_name = cell(length(h.att_ind),1);
for n=1:length(h.slo_ind)
	slo_name{n} = h.panel.inv_res(h.slo_ind(n)).name;
end
for n=1:length(h.att_ind)
	att_name{n} = h.panel.inv_res(h.att_ind(n)).name;
end
att_name{length(h.att_ind)+1} = '-';
set(handles.popupmenu_tomo_lent,'String',slo_name)
if ~isempty(att_name)
	set(handles.popupmenu_tomo_att,'String',att_name)
end

setappdata(handles.fig_interp, 'h', h)
update_plots(handles)



function set_string_locale(handles)
h.str = getappdata(handles.fig_interp, 'h');


function QuitMenuItem_Callback(hObject, eventdata, handles)
delete(handles.fig_interp)



function update_plots(handles)
cla(handles.axes1,'reset')
h = getappdata(handles.fig_interp, 'h');
if ~isfield(h,'slo_ind')
	return
end
if ishandle(h.cb), delete(h.cb); end
switch get(handles.popupmenu_prop,'Value')
	case 1
		plot_vitesse(handles)
	case 2
		plot_lenteur(handles)
	case 3
		plot_attenuation(handles)
	case 4
		plot_teneurEnEau(handles)
	case 5
		plot_constante_dielectrique(handles)
	case 6
		plot_permittivite(handles)
	case 7
		plot_conductivite(handles)
	case 8
		plot_resistivite(handles)
    case 9
        plot_angle_perte(handles)
end
xlabel(h.str.s119,'FontSize',12)
ylabel(h.str.s120,'FontSize',12)
ca = caxis(handles.axes1);
set(handles.edit_cmin,'String',num2str(ca(1)))
set(handles.edit_cmax,'String',num2str(ca(2)))



function plot_vitesse(handles)
h = getappdata(handles.fig_interp, 'h');
[s,x,z] = get_lenteur(handles);
vit = 1./s;
axes(handles.axes1)
imagesc(x,z,vit)
set(handles.axes1,'DataAspectRatio',[1 1 1],'YDir','normal')
h.cb = colorbar('peer',handles.axes1);
setappdata(handles.fig_interp, 'h', h)



function plot_lenteur(handles)
h = getappdata(handles.fig_interp, 'h');
[s,x,z] = get_lenteur(handles);
axes(handles.axes1)
imagesc(x,z,s)
set(handles.axes1,'DataAspectRatio',[1 1 1],'YDir','normal')
h.cb = colorbar('peer',handles.axes1);
setappdata(handles.fig_interp, 'h', h)



function plot_attenuation(handles)
h = getappdata(handles.fig_interp, 'h');
[a,x,z] = get_attenuation(handles);
if ~isempty(a)
	axes(handles.axes1)
	imagesc(x,z,a)
	set(handles.axes1,'DataAspectRatio',[1 1 1],'YDir','normal')
	h.cb = colorbar('peer',handles.axes1);
	setappdata(handles.fig_interp, 'h', h)
end


function plot_teneurEnEau(handles)
h = getappdata(handles.fig_interp, 'h');
switch get(handles.popupmenu_petro,'Value')
	case 1
		[theta,x,z] = get_Topp(handles);
	case 2
		[theta,x,z] = get_CRIM(handles);
	case 3
		[theta,x,z] = get_HanaiBruggeman(handles);
end
axes(handles.axes1)
imagesc(x,z,theta)
set(handles.axes1,'DataAspectRatio',[1 1 1],'YDir','normal')
h.cb = colorbar('peer',handles.axes1);
setappdata(handles.fig_interp, 'h', h)


function plot_permittivite(handles)
h = getappdata(handles.fig_interp, 'h');
[k,x,z] = get_K_dielectrique(handles);
eps0 = 8.8541878e-12;
e=k*eps0;
axes(handles.axes1)
imagesc(x,z,e)
set(handles.axes1,'DataAspectRatio',[1 1 1],'YDir','normal')
h.cb = colorbar('peer',handles.axes1);
setappdata(handles.fig_interp, 'h', h)


function plot_constante_dielectrique(handles)
h = getappdata(handles.fig_interp, 'h');
[k,x,z] = get_K_dielectrique(handles);
axes(handles.axes1)
imagesc(x,z,k)
set(handles.axes1,'DataAspectRatio',[1 1 1],'YDir','normal')
h.cb = colorbar('peer',handles.axes1);
setappdata(handles.fig_interp, 'h', h)


function plot_resistivite(handles)
h = getappdata(handles.fig_interp, 'h');
[s,x,z] = get_conductivite(handles);
if ~isempty(s)
	s = 1./s;
	axes(handles.axes1)
	imagesc(x,z,s)
	set(handles.axes1,'DataAspectRatio',[1 1 1],'YDir','normal')
	h.cb = colorbar('peer',handles.axes1);
	setappdata(handles.fig_interp, 'h', h)
end


function plot_conductivite(handles)
h = getappdata(handles.fig_interp, 'h');
[s,x,z] = get_conductivite(handles);
if ~isempty(s)
	axes(handles.axes1)
	imagesc(x,z,s)
	set(handles.axes1,'DataAspectRatio',[1 1 1],'YDir','normal')
	h.cb = colorbar('peer',handles.axes1);
	setappdata(handles.fig_interp, 'h', h)
end


function popupmenu_type_tomo_Callback(hObject, eventdata, handles)
update_plots(handles)


function popupmenu_type_tomo_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


function [s,x,z] = get_lenteur(handles)
h = getappdata(handles.fig_interp, 'h');
no = h.slo_ind( get(handles.popupmenu_tomo_lent,'Value') );
x = h.panel.inv_res(no).tomo.x;
z = h.panel.inv_res(no).tomo.z;
if get(handles.popupmenu_type_tomo,'Value')==1
	s = reshape(h.panel.inv_res(no).tomo.s,length(z),length(x));
else
	s = reshape(h.panel.inv_res(no).tomo.sgr,length(z),length(x));
	s = s + h.panel.inv_res(no).tomo.lmoy;
end



function [a,x,z] = get_attenuation(handles)
a = [];
x = [];
z = [];

noms = get(handles.popupmenu_tomo_att,'String');

if strcmp(noms{get(handles.popupmenu_tomo_att,'Value')},'-')==0
	h = getappdata(handles.fig_interp, 'h');
	no = h.att_ind( get(handles.popupmenu_tomo_att,'Value') );
	x = h.panel.inv_res(no).tomo.x;
	z = h.panel.inv_res(no).tomo.z;
	if get(handles.popupmenu_type_tomo,'Value')==1
		a = reshape(h.panel.inv_res(no).tomo.s,length(z),length(x));
	else
		a = reshape(h.panel.inv_res(no).tomo.sgr,length(z),length(x));
		a = a + h.panel.inv_res(no).tomo.lmoy;
	end
end



function [k,x,z] = get_K_dielectrique(handles)
[s,x,z] = get_lenteur(handles);    % ns/m
a = get_attenuation(handles);      % Np/m
if isempty(a)
	k = 0.2298^2 * s.^2;
else
	f = str2double(get(handles.edit_freq,'String'))*1e6;
	s = s*1e-9;
	mu0 = 4*pi*1e-7;
    eps0 = 8.8541878e-12;
	w = 2*pi*f;
	k = (1/(mu0*eps0))*(s.^2 - (a./w).^2);
end


function [sig,x,z] = get_conductivite(handles)
[s,x,z] = get_lenteur(handles);    % ns/m
a = get_attenuation(handles);      % Np/m
if ~isempty(a)
	s = s*1e-9;
	mu0 = 4*pi*1e-7;
	sig = 2*a.*s./mu0;
else
	sig = [];
end


function popupmenu_petro_Callback(hObject, eventdata, handles)
if get(hObject,'Value')==1
	set(handles.text_k_eau,'Visible','off')
	set(handles.text_k_matrice,'Visible','off')
	set(handles.edit_k_eau,'Visible','off')
	set(handles.edit_k_matrice,'Visible','off')
else
	set(handles.text_k_eau,'Visible','on')
	set(handles.text_k_matrice,'Visible','on')
	set(handles.edit_k_eau,'Visible','on')
	set(handles.edit_k_matrice,'Visible','on')
end
update_plots(handles)


function popupmenu_petro_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


function edit_k_eau_Callback(hObject, eventdata, handles)
update_plots(handles)


function edit_k_eau_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


function edit_k_matrice_Callback(hObject, eventdata, handles)
update_plots(handles)


function edit_k_matrice_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function [theta,x,z] = get_Topp(handles)
[k,x,z] = get_K_dielectrique(handles);
theta = topp(k);


function [theta,x,z] = get_CRIM(handles)
[k,x,z] = get_K_dielectrique(handles);
k_w = str2double(get(handles.edit_k_eau,'String'));
k_m = str2double(get(handles.edit_k_matrice,'String'));
theta = (sqrt(k) - sqrt(k_m))./(sqrt(k_w) - sqrt(k_m));


function [theta,x,z] = get_HanaiBruggeman(handles)
[k,x,z] = get_K_dielectrique(handles);
k_w = str2double(get(handles.edit_k_eau,'String'));
k_m = str2double(get(handles.edit_k_matrice,'String'));
W = 1/3;
theta = 1-((k_w-k)./(k_w-k_m)).*(k_m./k).^W;


function edit_freq_Callback(hObject, eventdata, handles)
update_plots(handles)


function edit_freq_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


function plot_angle_perte(handles)
sig = get_conductivite(handles);
if ~isempty(sig)
    [k,x,z] = get_K_dielectrique(handles);
    f = str2double(get(handles.edit_freq,'String'))*1e6;
    w = 2*pi*f;
    eps0 = 8.8541878e-12;
    e=k*eps0;
    delta = sig./(w*e);
    axes(handles.axes1)
    imagesc(x,z,delta)
    set(handles.axes1,'DataAspectRatio',[1 1 1],'YDir','normal')
    h.cb = colorbar('peer',handles.axes1);
    setappdata(handles.fig_interp, 'h', h)
end

function edit_cmin_Callback(hObject, eventdata, handles)
ca = caxis(handles.axes1);
ca(1) = str2double(get(handles.edit_cmin,'String'));
caxis(handles.axes1, ca);


function edit_cmin_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


function edit_cmax_Callback(hObject, eventdata, handles)
ca = caxis(handles.axes1);
ca(2) = str2double(get(handles.edit_cmax,'String'));
caxis(handles.axes1, ca);


function edit_cmax_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
