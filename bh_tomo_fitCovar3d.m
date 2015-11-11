function varargout = bh_tomo_fitCovar3d(varargin)
% BH_TOMO_FITCOVAR3D M-file for bh_tomo_fitCovar3d.fig
%      BH_TOMO_FITCOVAR3D, by itself, creates a new BH_TOMO_FITCOVAR3D or raises the existing
%      singleton*.
%
%      H = BH_TOMO_FITCOVAR3D returns the handle to a new BH_TOMO_FITCOVAR3D or the handle to
%      the existing singleton*.
%
%      BH_TOMO_FITCOVAR3D('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in BH_TOMO_FITCOVAR3D.M with the given input arguments.
%
%      BH_TOMO_FITCOVAR3D('Property','Value',...) creates a new BH_TOMO_FITCOVAR3D or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before bh_tomo_fitCovar_OpeningFunction gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to bh_tomo_fitCovar3d_OpeningFcn via varargin.
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

% Edit the above text to modify the response to help bh_tomo_fitCovar3d

% Last Modified by GUIDE v2.5 14-Feb-2013 15:22:28

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @bh_tomo_fitCovar3d_OpeningFcn, ...
                   'gui_OutputFcn',  @bh_tomo_fitCovar3d_OutputFcn, ...
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


% --- Executes just before bh_tomo_fitCovar3d is made visible.
function bh_tomo_fitCovar3d_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to bh_tomo_fitCovar3d (see VARARGIN)

% Choose default command line output for bh_tomo_fitCovar3d
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

h.db_file = get(handles.fig_bh_FitCov,'UserData');
h.saved = true;
h.data = [];
h.xmin = [];
h.xmax = [];
h.dx = 0.4;
h.ymin = [];
h.ymax = [];
h.dy = 0.4;
h.zmin = [];
h.zmax = [];
h.dz = 0.4;
h.anglemin = 100;
h.Dir1=1;
h.Dir2=3;
h.L = [];
h.gridx = [];
h.gridy = [];
h.gridz = [];
h.iktt = [];
h.model = [2 4 4 4 0 0 0];
h.c = 1;
h.model_xi = [2 4 4 4 0 0 0];
h.c_xi = 1;
h.model_th = [2 50 50 0];
h.c_th = 0.001;
h.nugget_t = 2;
h.nugget_l = 0;
h.nugget_t_code = 0;
h.nugget_l_code = 0;
h.nugget_xi = 0;
h.nugget_th = 0;
h.nugget_xi_code = 0;
h.nugget_th_code = 0;
h.model_code=[1 0 0 0];
h.c_code=0;
h.model_xi_code=[1 0 0 0];
h.c_xi_code=0;
h.model_th_code=[1 0 0 0];
h.c_th_code=0;
h.afi=4;
h.lclas=50;  %1000;
h.c0 = [];
h.options=optimset('TolX',1e-12,'TolFun',1e-12,'MaxFunEvals',5,'Display','off');
h.aniso = 0;
h.str = get_str_locale();
set_String_locale(handles, h.str)

set(handles.dx,'String',num2str(h.dx));
set(handles.dy,'String',num2str(h.dy));
set(handles.dz,'String',num2str(h.dz));
set(handles.numItSimplex,'String','5');
set(handles.lclas,'String',num2str(h.lclas));
set(handles.afi,'String',num2str(1/h.afi));
hl=get(handles.axes1,'XLabel');
set(hl,'String', h.str.s108, 'FontSize',12);
hl=get(handles.axes1,'YLabel');
set(hl,'String', h.str.s109, 'FontSize',12);
hl=get(handles.axes2,'XLabel');
set(hl,'String', h.str.s110, 'FontSize',12);
hl=get(handles.axes2,'YLabel');
set(hl,'String', h.str.s111, 'FontSize',12);
set(handles.popup_model11,'Value',h.model(1,1)-1); % Nugget (effect (code 1) not included in list

setappdata(handles.fig_bh_FitCov, 'h_GUIdata', h);
updateCovAffichage(handles)

% UIWAIT makes bh_tomo_fitCovar3d wait for user response (see UIRESUME)
% uiwait(handles.fig_bh_FitCov);


% --- Outputs from this function are returned to the command line.
function varargout = bh_tomo_fitCovar3d_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


function OpenMenuItem_Callback(hObject, eventdata, handles)
h = getappdata(handles.fig_bh_FitCov, 'h_GUIdata');
h.data = [];
h.xmin = [];
h.xmax = [];
%h.dx = 0.4;
h.zmin = [];
h.zmax = [];
%h.dz = 0.4;
h.ymin = [];
h.ymax = [];
%h.dy = 0.4;
h.anglemin = 100;
h.Dir1=1;
h.Dir2=3;
h.L = [];
h.gridx = [];
h.gridy = [];
h.gridz = [];
h.iktt = [];
h.model = [2 4 4 4 0 0 0];
h.c = 1;
h.model_xi = [2 4 4 4 0 0 0];
h.c_xi = 1;
h.model_th = [2 50 50 0];
h.c_th = 0.001;
h.model_code=[1 0 0 0];
h.c_code = 0;
h.nugget_t = 2;
h.nugget_l = 0;
h.nugget_t_code = 0;
h.nugget_l_code = 0;
h.nugget_xi = 0;
h.nugget_th = 0;
h.nugget_xi_code = 0;
h.nugget_th_code = 0;
h.model_xi_code=[1 0 0 0];
h.c_xi_code = 0;
h.model_th_code=[1 0 0 0];
h.c_th_code = 0;
h.afi=4;
h.lclas=50;%1000;
h.c0 = [];
h.aniso = 0;

[h.no_model3d,h.db_file,hh] = choisirPanneau('UserData', h.db_file,'3d' );
if ishandle(hh), delete(hh), end
if h.no_model3d==0
	return
end
load(h.db_file,'model3d')
h.model3d = model3d(h.no_model3d);
if isempty(h.model3d.grid3d)
    errordlg(h.str.s150, s.str.s45,'modal')
    return
end
set(handles.NomData,'String',h.model3d.name);
load(h.db_file,'mogs')
mog_date = {};
%for i = 1:length(mogs)  BG: changed to load mogs of currently load pannel
for i = h.model3d.mogs
    mog_date{i} = [mogs(i).name, ',', mogs(i).date];
end
set(handles.listbox_mogs,'Value',1,'String',mog_date);
selected_mogs = get(handles.listbox_mogs,'Value');

h.xmin = h.model3d.grid3d.grx(1);
h.xmax = h.model3d.grid3d.grx(length(h.model3d.grid3d.grx));
h.ymin = h.model3d.grid3d.gry(1);
h.ymax = h.model3d.grid3d.gry(length(h.model3d.grid3d.gry));
h.zmin = h.model3d.grid3d.grz(1);
h.zmax = h.model3d.grid3d.grz(length(h.model3d.grid3d.grz));

h.dx = 2*(h.model3d.grid3d.grx(2)-h.model3d.grid3d.grx(1));
h.dy = 2*(h.model3d.grid3d.gry(2)-h.model3d.grid3d.gry(1));
h.dz = 2*(h.model3d.grid3d.grz(2)-h.model3d.grid3d.grz(1));

set(handles.xmin,'String',num2str(h.xmin));
set(handles.zmin,'String',num2str(h.zmin));
set(handles.xmax,'String',num2str(h.xmax));
set(handles.zmax,'String',num2str(h.zmax));
set(handles.dx,'String',num2str(h.dx));
set(handles.dz,'String',num2str(h.dz));
set(handles.ymin,'String',num2str(h.ymin));
set(handles.ymax,'String',num2str(h.ymax));
set(handles.dy,'String',num2str(h.dy));

if get(handles.popupmenu_type_data,'Value')==1
    if get(handles.checkbox_lim_vapp,'Value')==1
        vlim = str2double(get(handles.edit_lim_vapp,'String'));
        selected_mogs = 1:length(get(handles.listbox_mogs,'String'));
    else
        vlim = [];
    end
    [h.data,h.ind] = getPanneauData(h.model3d, h.db_file,'tt',selected_mogs, vlim);  
    if isfield(h.model3d,'tt_covar')
		if ~isempty( h.model3d.tt_covar )
			h.model = h.model3d.tt_covar.model;
			h.c = h.model3d.tt_covar.c;
			h.nugget_t = h.model3d.tt_covar.nugget_t;
			h.nugget_l = h.model3d.tt_covar.nugget_l;
            if isfield( h.model3d.tt_covar, 'aniso' )
                h.aniso = h.model3d.tt_covar.aniso;
                if h.aniso > 0 || isfield(h.model3d.tt_covar,'model_xi')
                    h.model_xi = h.model3d.tt_covar.model_xi;
                    h.c_xi = h.model3d.tt_covar.c_xi;
                    h.nugget_xi = h.model3d.tt_covar.nugget_xi;
                    if h.aniso > 1 || isfield(h.model3d.tt_covar,'model_th')
                        h.model_th = h.model3d.tt_covar.model_th;
                        h.c_th = h.model3d.tt_covar.c_th;
                        h.nugget_th = h.model3d.tt_covar.nugget_th;
                    end
                end
            end
			set(handles.checkbox_c0,'Value',h.model3d.tt_covar.use_c0)
		end
    end
else
	if get(handles.popupmenu_type_data,'Value')==2
		[h.data,h.ind] = getPanneauData(h.model3d, h.db_file,'amp',selected_mogs);
	elseif get(handles.popupmenu_type_data,'Value')==3
		[h.data,h.ind] = getPanneauData(h.model3d, h.db_file,'fce',selected_mogs); 
	elseif get(handles.popupmenu_type_data,'Value')==4
		[h.data,h.ind] = getPanneauData(h.model3d, h.db_file,'hyb',selected_mogs);
	end
    if isfield(h.model3d,'amp_covar')
		if ~isempty( h.model3d.amp_covar )
			h.model = h.model3d.amp_covar.model;
			h.c = h.model3d.amp_covar.c;
			h.nugget_t = h.model3d.amp_covar.nugget_t;
			h.nugget_l = h.model3d.amp_covar.nugget_l;
			set(handles.checkbox_c0,'Value',0)
		end
    end
end

if h.data(:,2)~=-1
    set(handles.checkbox_c0,'Enable','on')
end
texte{1} = '-';
if isfield(h.model3d,'inv_res')
	for n=1:length(h.model3d.inv_res)
		texte{n+1} = [char( h.model3d.inv_res(n).name ), ', ',char( h.model3d.inv_res(n).tomo.date)]; 
	end
end
set(handles.popupmenu_Ldc,'String',texte,'Value',1)  

set(handles.nombre_rais,'String',[num2str(size(h.data,1)),' ',lower(h.str.s64)])
setappdata(handles.fig_bh_FitCov, 'h_GUIdata', h);
updateCovAffichage(handles)
CalculCovExperimentale(handles)
updateFigures(handles)


function FileMenu_Callback(hObject, eventdata, handles)


function xmin_Callback(hObject, eventdata, handles)
val = str2double(get(hObject,'string'));
if isnan(val)
    h = getappdata(handles.fig_bh_FitCov, 'h_GUIdata');
    errordlg(h.str.s54, h.str.s45,'modal')
	return
end
h = getappdata(handles.fig_bh_FitCov, 'h_GUIdata');
h.xmin = val;
setappdata(handles.fig_bh_FitCov, 'h_GUIdata', h);
CalculRaisDroits(handles)
CalculCovExperimentale(handles)
updateFigures(handles)


function xmin_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


function xmax_Callback(hObject, eventdata, handles)
val = str2double(get(hObject,'string'));
if isnan(val)
    h = getappdata(handles.fig_bh_FitCov, 'h_GUIdata');
    errordlg(h.str.s54, h.str.s45,'modal')
	return
end
h = getappdata(handles.fig_bh_FitCov, 'h_GUIdata');
h.xmax = val;
setappdata(handles.fig_bh_FitCov, 'h_GUIdata', h);
CalculRaisDroits(handles)
CalculCovExperimentale(handles)
updateFigures(handles)


function xmax_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


function dx_Callback(hObject, eventdata, handles)
val = str2double(get(hObject,'string'));
if isnan(val)
    h = getappdata(handles.fig_bh_FitCov, 'h_GUIdata');
    errordlg(h.str.s54, h.str.s45,'modal')
	return
end
h = getappdata(handles.fig_bh_FitCov, 'h_GUIdata');
h.dx = val;
if ~isempty(h.xmax)
  tmp = h.xmin:val:h.xmax;
  if tmp(end)<h.xmax
	h.xmax = tmp(end)+val;
	set(handles.xmax,'string',num2str(h.xmax));
  end
end
setappdata(handles.fig_bh_FitCov, 'h_GUIdata', h);
CalculRaisDroits(handles)
CalculCovExperimentale(handles)
if ~isempty(h.data)
  updateFigures(handles)
end


function dx_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function ymin_Callback(hObject, eventdata, handles)
val = str2double(get(hObject,'string'));
if isnan(val)
    h = getappdata(handles.fig_bh_FitCov, 'h_GUIdata');
    errordlg(h.str.s54, h.str.s45,'modal')
	return
end
h = getappdata(handles.fig_bh_FitCov, 'h_GUIdata');
h.ymin = val;
setappdata(handles.fig_bh_FitCov, 'h_GUIdata', h);
CalculRaisDroits(handles)
CalculCovExperimentale(handles)
updateFigures(handles)

function ymin_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function ymax_Callback(hObject, eventdata, handles)
val = str2double(get(hObject,'string'));
if isnan(val)
    h = getappdata(handles.fig_bh_FitCov, 'h_GUIdata');
    errordlg(h.str.s54, h.str.s45,'modal')
	return
end
h = getappdata(handles.fig_bh_FitCov, 'h_GUIdata');
h.ymax = val;
setappdata(handles.fig_bh_FitCov, 'h_GUIdata', h);
CalculRaisDroits(handles)
CalculCovExperimentale(handles)
updateFigures(handles)
function ymax_CreateFcn(hObject, eventdata, handles)

if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


function dy_Callback(hObject, eventdata, handles)
val = str2double(get(hObject,'string'));
if isnan(val)
    h = getappdata(handles.fig_bh_FitCov, 'h_GUIdata');
    errordlg(h.str.s54, h.str.s45,'modal')
	return
end
h = getappdata(handles.fig_bh_FitCov, 'h_GUIdata');
h.dy = val;
if ~isempty(h.ymax)
  tmp = h.ymin:val:h.ymax;
  if tmp(end)<h.ymax
	h.ymax = tmp(end)+val;
	set(handles.ymax,'string',num2str(h.ymax));
  end
end
setappdata(handles.fig_bh_FitCov, 'h_GUIdata', h);
CalculRaisDroits(handles)
CalculCovExperimentale(handles)
if ~isempty(h.data)
  updateFigures(handles)
end

function dy_CreateFcn(hObject, eventdata, handles)

if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function zmin_Callback(hObject, eventdata, handles)
val = str2double(get(hObject,'string'));
if isnan(val)
    h = getappdata(handles.fig_bh_FitCov, 'h_GUIdata');
    errordlg(h.str.s54, h.str.s45,'modal')
	return
end
h = getappdata(handles.fig_bh_FitCov, 'h_GUIdata');
h.zmin = val;
setappdata(handles.fig_bh_FitCov, 'h_GUIdata', h);
CalculRaisDroits(handles)
CalculCovExperimentale(handles)
updateFigures(handles)


function zmin_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


function zmax_Callback(hObject, eventdata, handles)
val = str2double(get(hObject,'string'));
if isnan(val)
    h = getappdata(handles.fig_bh_FitCov, 'h_GUIdata');
    errordlg(h.str.s54, h.str.s45,'modal')
	return
end
h = getappdata(handles.fig_bh_FitCov, 'h_GUIdata');
h.zmax = val;
setappdata(handles.fig_bh_FitCov, 'h_GUIdata', h);
CalculRaisDroits(handles)
CalculCovExperimentale(handles)
updateFigures(handles)


function zmax_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


function dz_Callback(hObject, eventdata, handles)
val = str2double(get(hObject,'string'));
if isnan(val)
    h = getappdata(handles.fig_bh_FitCov, 'h_GUIdata');
    errordlg(h.str.s54, h.str.s45,'modal')
	return
end
h = getappdata(handles.fig_bh_FitCov, 'h_GUIdata');
h.dz = val;
if ~isempty(h.zmax)
  tmp = h.zmin:val:h.zmax;
  if tmp(end)<h.zmax
	h.zmax = tmp(end)+val;
	set(handles.zmax,'string',num2str(h.zmax));
  end
end
setappdata(handles.fig_bh_FitCov, 'h_GUIdata', h);
CalculRaisDroits(handles)
CalculCovExperimentale(handles)
if ~isempty(h.data)
  updateFigures(handles)
end


function dz_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


function CloseMenuItem_Callback(hObject, eventdata, handles)
h = getappdata(handles.fig_bh_FitCov, 'h_GUIdata');
if h.saved == false
	ButtonName=questdlg(h.str.s236);
	switch ButtonName,
		case 'Yes',
			SaveMenuItem_Callback(hObject, eventdata, handles)
		case 'No',
		case 'Cancel',
			return
	end % switch
end
delete(handles.fig_bh_FitCov)


function CalculCovExperimentale(handles)
h = getappdata(handles.fig_bh_FitCov, 'h_GUIdata');
if isempty(h.L)==1
	CalculRaisDroits(handles);
	h = getappdata(handles.fig_bh_FitCov, 'h_GUIdata');
	if isempty(h.L)==1, return, end
end
nt=length(h.L(:,1));

%Calcul de la lenteur moyenne
if h.aniso==0
    s0 = mean(h.data(:,1)./sum(h.L,2));
    mta = s0*sum(h.L,2);
else
    np = size(h.L,2)/2;
    l = sqrt(h.L(:,1:np).^2 + h.L(:,(np+1):end).^2);
    s0 = mean(h.data(:,1)./sum(l,2));
    mta = s0 * sum(l,2);
end
%vecteur des delta temps (t-moy(t))
dt = h.data(:,1)-mta;

%Matrice de covariance experimentale des dt: dt*dt'
%mise en vecteur
h.iktt=reshape(dt*dt',nt^2,1);

setappdata(handles.fig_bh_FitCov, 'h_GUIdata', h);


function CalculRaisDroits(handles)
h = getappdata(handles.fig_bh_FitCov, 'h_GUIdata');
if isempty(h.data)
    warndlg(h.str.s52, h.str.s53)
    return
end
grx = h.xmin:h.dx:h.xmax;
if (grx(end)<h.xmax)
  h.xmax = grx(end)+h.dx;
  grx = h.xmin:h.dx:h.xmax;
end
gry = h.ymin:h.dy:h.ymax;
if (gry(end)<h.ymax)
  h.xmax = grx(end)+h.dx;
  gry = h.ymin:h.dy:h.ymax;
end
grz = h.zmin:h.dz:h.zmax;
if (grz(end)<h.zmax)
  h.zmax = grz(end)+h.dz;
  grz = h.zmin:h.dz:h.zmax;
end
if h.aniso == 0
    [h.L,h.gridx,h.gridy, h.gridz] = Lsr3d(h.model3d.grid3d.Tx(h.ind,[1 2 3]),...
        h.model3d.grid3d.Rx(h.ind,[1 2 3]), grx, gry,grz);
else
    [h.L,h.gridx,h.gridy,h.gridz] = Lsr3da(h.model3d.grid3d.Tx(h.ind,[1 2 3]),...
        h.model3d.grid3d.Rx(h.ind,[1 2 3]), grx, gry,grz);
    h.L = abs(h.L);
end

n = length(h.gridz);
m = length(h.gridx);
l = length(h.gridy);
set(handles.text_nombre_cellules,'String',[num2str(n*m), ' ',h.str.s205])
% h.x: contient les coordonnees (z,x) des points de la grille
h.x = grille3d(min(h.gridz),h.gridz(2)-h.gridz(1),n,...
               min(h.gridy),h.gridy(2)-h.gridy(1),l,...
               min(h.gridx),h.gridx(2)-h.gridx(1),m);

setappdata(handles.fig_bh_FitCov, 'h_GUIdata', h);


function model14_Callback(hObject, eventdata, handles)
val = str2double(get(hObject,'string'));
if isnan(val)
    h = getappdata(handles.fig_bh_FitCov, 'h_GUIdata');
    errordlg(h.str.s54, h.str.s45,'modal')
	return
end
h = getappdata(handles.fig_bh_FitCov, 'h_GUIdata');
h.model(1,2) = val;
h.saved = false;
setappdata(handles.fig_bh_FitCov, 'h_GUIdata', h);
if get(handles.checkbox_calc_auto,'Value')==1, updateFigures(handles), end


function model14_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


function model12_Callback(hObject, eventdata, handles)
val = str2double(get(hObject,'string'));
if isnan(val)
    h = getappdata(handles.fig_bh_FitCov, 'h_GUIdata');
    errordlg(h.str.s54, h.str.s45,'modal')
	return
end
h = getappdata(handles.fig_bh_FitCov, 'h_GUIdata');
h.model(1,3) = val;
h.saved = false;
setappdata(handles.fig_bh_FitCov, 'h_GUIdata', h);
if get(handles.checkbox_calc_auto,'Value')==1, updateFigures(handles), end


function model12_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


function model15_Callback(hObject, eventdata, handles)
val = str2double(get(hObject,'string'));
if isnan(val)
    h = getappdata(handles.fig_bh_FitCov, 'h_GUIdata');
    errordlg(h.str.s54, h.str.s45,'modal')
	return
end
h = getappdata(handles.fig_bh_FitCov, 'h_GUIdata');
h.model(1,4) = val;
h.saved = false;
setappdata(handles.fig_bh_FitCov, 'h_GUIdata', h);
if get(handles.checkbox_calc_auto,'Value')==1, updateFigures(handles), end


function model15_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


function c1_Callback(hObject, eventdata, handles)
val = str2double(get(hObject,'string'));
if isnan(val)
    h = getappdata(handles.fig_bh_FitCov, 'h_GUIdata');
    errordlg(h.str.s54, h.str.s45,'modal')
	return
end
h = getappdata(handles.fig_bh_FitCov, 'h_GUIdata');
h.c(1,1) = val;
h.saved = false;
setappdata(handles.fig_bh_FitCov, 'h_GUIdata', h);
if get(handles.checkbox_calc_auto,'Value')==1, updateFigures(handles), end


function c1_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


function edit_pepite_t_Callback(hObject, eventdata, handles)
val = str2double(get(hObject,'string'));
if isnan(val)
    h = getappdata(handles.fig_bh_FitCov, 'h_GUIdata');
    errordlg(h.str.s54, h.str.s45,'modal')
	return
end
h = getappdata(handles.fig_bh_FitCov, 'h_GUIdata');
h.nugget_t = val;
h.saved = false;
setappdata(handles.fig_bh_FitCov, 'h_GUIdata', h);
if get(handles.checkbox_calc_auto,'Value')==1, updateFigures(handles), end


function edit_pepite_t_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


function lclas_Callback(hObject, eventdata, handles)
val = str2double(get(hObject,'string'));
if isnan(val)
    h = getappdata(handles.fig_bh_FitCov, 'h_GUIdata');
    errordlg(h.str.s54, h.str.s45,'modal')
	return
end
h = getappdata(handles.fig_bh_FitCov, 'h_GUIdata');
h.lclas = val;
setappdata(handles.fig_bh_FitCov, 'h_GUIdata', h);
updateFigures(handles)


function lclas_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


function updateFigures(handles)
h = getappdata(handles.fig_bh_FitCov, 'h_GUIdata');
hm = msgbox(h.str.s235); drawnow

covar.model = h.model;
covar.c = h.c;
covar.model_xi = h.model_xi;
covar.c_xi = h.c_xi;
covar.model_th = h.model_th;
covar.c_th = h.c_th;
covar.nugget_t = h.nugget_t;
covar.nugget_l = h.nugget_l;
covar.nugget_xi = h.nugget_xi;
covar.nugget_th = h.nugget_th;
covar.use_c0 = get(handles.checkbox_c0,'Value');
covar.aniso = h.aniso;

id{1} = h.model_code==0;
x0libre = h.model(id{1});
x0libre = x0libre(:);  % on s'assure qu'on a un vecteur colonne

id{2} = h.c_code==0;
x0libre = [x0libre; h.c(id{2})];

id{3} = h.nugget_t_code==0;
x0libre = [x0libre; h.nugget_t(id{3})];

id{4} = h.nugget_l_code==0;
x0libre = [x0libre; h.nugget_l(id{4})];

if covar.aniso > 0
    id{5} = h.model_xi_code==0;
    tmp = covar.model_xi(id{5});
    x0libre = [x0libre; tmp(:)];
    
    id{6} = h.c_xi_code==0;
    x0libre = [x0libre; covar.c_xi(id{6})];

    id{7} = h.nugget_xi_code==0;
    x0libre = [x0libre; h.nugget_xi(id{7})];
end
if covar.aniso > 1
    id{8} = h.model_th_code==0;
    tmp = covar.model_th(id{8});
    x0libre = [x0libre; tmp(:)];
    
    id{9} = h.c_th_code==0;
    x0libre = [x0libre; covar.c_th(id{9})];

    id{10} = h.nugget_th_code==0;
    x0libre = [x0libre; h.nugget_th(id{10})];
end


if h.aniso == 0
    modeliKss3d(x0libre,covar,id,...
        h.L,h.x,h.iktt,h.afi,h.lclas,h.c0,handles.axes1,...
        handles.axes2);
elseif h.aniso == 1
    np = size(h.L,2)/2;
    l = sqrt(h.L(:,1:np).^2 + h.L(:,(np+1):end).^2);
    s0 = mean(h.data(:,1)./sum(l,2))+zeros(np,1);
    xi0 = ones(np,1);
    J = calculJ(h.L, [s0; xi0]);
    modeliKss3d(x0libre,covar,id,...
        J,h.x,h.iktt,h.afi,h.lclas,h.c0,handles.axes1,...
        handles.axes2);
elseif h.aniso == 2
    np = size(h.L,2)/2;
    l = sqrt(h.L(:,1:np).^2 + h.L(:,(np+1):end).^2);
    s0 = mean(h.data(:,1)./sum(l,2))+zeros(np,1);
    xi0 = ones(np,1)+0.010;    % add 1/1000 so that J_th ~= 0
    theta0 = zeros(np,1)+0.0044;  % add a quarter of a degree so that J_th ~= 0
    J = calculJ2(h.L, [s0; xi0; theta0]);
    modeliKss3d(x0libre,covar,id,...
        J,h.x,h.iktt,h.afi,h.lclas,h.c0,handles.axes1,...
        handles.axes2);
end

legend(handles.axes2,h.str.s206,h.str.s207)
delete(hm)
setappdata(handles.fig_bh_FitCov, 'h_GUIdata', h)


function numItSimplex_Callback(hObject, eventdata, handles)
val = str2double(get(hObject,'string'));
if isnan(val)
    h = getappdata(handles.fig_bh_FitCov, 'h_GUIdata');
    errordlg(h.str.s54, h.str.s45,'modal')
	return
end
h = getappdata(handles.fig_bh_FitCov, 'h_GUIdata');
h.options=optimset('TolX',1e-12,'TolFun',1e-12,'MaxFunEvals',val,'Display','off');
setappdata(handles.fig_bh_FitCov, 'h_GUIdata', h);


function numItSimplex_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function load_selected_mogs(hObject, eventdata, handles)
selected_mogs = get(handles.listbox_mogs,'Value');
h = getappdata(handles.fig_bh_FitCov, 'h_GUIdata');

if get(handles.popupmenu_type_data,'Value')==1
		[h.data,h.ind] = getPanneauData(h.model3d, h.db_file,'tt',selected_mogs); %YH
elseif get(handles.popupmenu_type_data,'Value')==2
		[h.data,h.ind] = getPanneauData(h.model3d, h.db_file,'amp',selected_mogs); %YH
elseif get(handles.popupmenu_type_data,'Value')==3
		[h.data,h.ind] = getPanneauData(h.model3d, h.db_file,'fce',selected_mogs); %YH
elseif get(handles.popupmenu_type_data,'Value')==4
		[h.data,h.ind] = getPanneauData(h.model3d, h.db_file,'hyb',selected_mogs); %YH
end
setappdata(handles.fig_bh_FitCov, 'h_GUIdata', h);

function pushbutton_ajuster_Callback(hObject, eventdata, handles)
load_selected_mogs(hObject, eventdata, handles);
h = getappdata(handles.fig_bh_FitCov, 'h_GUIdata');
hm=msgbox(h.str.s113); drawnow
if size(h.model)~=size(h.model_code)
	% on avait seulement une structure
	h.model_code=[h.model_code; 1 0 0 0];
	h.c_code=[h.c_code; 0];
	h.model_xi_code=[h.model_xi_code; 1 0 0 0];
	h.c_xi_code=[h.c_xi_code; 0];
	h.model_th_code=[h.model_th_code; 1 0 0 0];
	h.c_th_code=[h.c_th_code; 0];
end

covar.model = h.model;
covar.c = h.c;
covar.model_xi = h.model_xi;
covar.c_xi = h.c_xi;
covar.model_th = h.model_th;
covar.c_th = h.c_th;
covar.nugget_t = h.nugget_t;
covar.nugget_l = h.nugget_l;
covar.nugget_xi = h.nugget_xi;
covar.nugget_th = h.nugget_th;
covar.use_c0 = get(handles.checkbox_c0,'Value');
covar.aniso = h.aniso;

covar_code.model = h.model_code;
covar_code.c = h.c_code;
covar_code.model_xi = h.model_xi_code;
covar_code.c_xi = h.c_xi_code;
covar_code.model_th = h.model_th_code;
covar_code.c_th = h.c_th_code;
covar_code.nugget_t = h.nugget_t_code;
covar_code.nugget_l = h.nugget_l_code;
covar_code.nugget_xi = h.nugget_xi_code;
covar_code.nugget_th = h.nugget_th_code;

if h.aniso == 0
[covar,dif]=ajuster_Cov3d(h.options,covar,...
	covar_code,h.L,h.gridx,h.gridz,h.x,h.iktt,...
	h.afi,h.lclas,h.c0,handles.axes1,handles.axes2);
elseif h.aniso == 1
    np = size(h.L,2)/2;
    l = sqrt(h.L(:,1:np).^2 + h.L(:,(np+1):end).^2);
    s0 = mean(h.data(:,1)./sum(l,2))+zeros(np,1);
    xi0 = ones(np,1);
    J = calculJ(h.L, [s0; xi0]);
    [covar,dif]=ajuster_Cov3d(h.options,covar,...
        covar_code,J,h.gridx,h.gridz,h.x,h.iktt,...
        h.afi,h.lclas,h.c0,handles.axes1,handles.axes2);
elseif h.aniso == 2
    np = size(h.L,2)/2;
    l = sqrt(h.L(:,1:np).^2 + h.L(:,(np+1):end).^2);
    s0 = mean(h.data(:,1)./sum(l,2))+zeros(np,1);
    xi0 = ones(np,1)+0.01;
    theta0 = zeros(np,1)+0.0087;
    J = calculJ2(h.L, [s0; xi0; theta0]);
    [covar,dif]=ajuster_Cov3d(h.options,covar,...
        covar_code,J,h.gridx,h.gridz,h.x,h.iktt,...
        h.afi,h.lclas,h.c0,handles.axes1,handles.axes2);
end

h.model = covar.model;
h.c = covar.c;
h.nugget_t = covar.nugget_t;
h.nugget_l = covar.nugget_l;
if h.aniso > 0
    h.model_xi = covar.model_xi;
    h.c_xi = covar.c_xi;
    h.nugget_xi = covar.nugget_xi;
end
if h.aniso > 1
    h.model_th = covar.model_th;
    h.c_th = covar.c_th;
    h.nugget_th = covar.nugget_th;
end
h.saved = false;
setappdata(handles.fig_bh_FitCov, 'h_GUIdata', h);
updateCovAffichage(handles)
delete(hm)


function updateCovAffichage(handles)
h = getappdata(handles.fig_bh_FitCov, 'h_GUIdata');
set(handles.popup_model11,'Value',h.model(1,1)-1) % Nugget (effect (code 1) not included in list
set(handles.model12,'String',num2str(h.model(1,2)));
set(handles.model13,'String',num2str(h.model(1,3)));
set(handles.model14,'String',num2str(h.model(1,4)));
set(handles.model15,'String',num2str(h.model(1,5)));
set(handles.model16,'String',num2str(h.model(1,6)));
set(handles.model17,'String',num2str(h.model(1,7)));
set(handles.c1,'String',num2str(h.c(1,1)));
set(handles.edit_pepite_t,'String',num2str(h.nugget_t));
set(handles.edit_pepite_l,'String',num2str(h.nugget_l));
if size(h.model,1)==2
	set(handles.popup_model21,'Value',h.model(2,1)-1) % Nugget (effect (code 1) not included in list
	set(handles.model22,'String',num2str(h.model(2,2)));
    set(handles.model23,'String',num2str(h.model(2,3)));
    set(handles.model24,'String',num2str(h.model(2,4)));
    set(handles.model25,'String',num2str(h.model(2,5)));
    set(handles.model26,'String',num2str(h.model(2,6)));
    set(handles.model27,'String',num2str(h.model(2,7)));
	set(handles.c2,'String',num2str(h.c(2,1)));
end
if h.aniso > 0
    
    set(handles.text_lenteur1, 'String', h.str.s283)
    set(handles.text_lenteur2, 'String', h.str.s283)
    set(handles.text_xi1, 'String', h.str.s284)
    set(handles.text_xi2, 'String', h.str.s284)
    set(handles.popup_model11x,'Visible','on')
    set(handles.model12x,'Visible','on')
    set(handles.model13x,'Visible','on')
    set(handles.model14x,'Visible','on')
    set(handles.c1x,'Visible','on')
    set(handles.popup_model21x,'Visible','on')
    set(handles.model22x,'Visible','on')
    set(handles.model23x,'Visible','on')
    set(handles.model24x,'Visible','on')
    set(handles.c2x,'Visible','on')
    set(handles.text_xi1,'Visible','on')
    set(handles.text_xi2,'Visible','on')
    set(handles.text_fixe1_x,'Visible','on')
    set(handles.text_fixe2_x,'Visible','on')
    set(handles.checkbox_mod12x,'Visible','on')
    set(handles.checkbox_mod13x,'Visible','on')
    set(handles.checkbox_mod14x,'Visible','on')
    set(handles.checkbox_c1x,'Visible','on')
    set(handles.checkbox_mod22x,'Visible','on')
    set(handles.checkbox_mod23x,'Visible','on')
    set(handles.checkbox_mod24x,'Visible','on')
    set(handles.checkbox_c2x,'Visible','on')
    set(handles.text_pepite_xi,'Visible','on')
    set(handles.edit_pepite_xi,'Visible','on')
    set(handles.checkbox_pepite_xi,'Visible','on')
    
    set(handles.popupmenu_Ldc,'Value',1);
    
    set(handles.popup_model11x,'Value',h.model_xi(1,1)-1) % Nugget (effect (code 1) not included in list
    set(handles.model12x,'String',num2str(h.model_xi(1,2)));
    set(handles.model13x,'String',num2str(h.model_xi(1,3)));
    set(handles.model14x,'String',num2str(h.model_xi(1,4)));
    set(handles.c1x,'String',num2str(h.c_xi(1,1)));
    set(handles.edit_pepite_xi,'String',num2str(h.nugget_xi));
    if size(h.model_xi,1)==2
        set(handles.popup_model21x,'Value',h.model_xi(2,1)-1) % Nugget (effect (code 1) not included in list
        set(handles.model22x,'String',num2str(h.model_xi(2,2)));
        set(handles.model23x,'String',num2str(h.model_xi(2,3)));
        set(handles.model24x,'String',num2str(h.model_xi(2,4)));
        set(handles.c2x,'String',num2str(h.c_xi(2,1)));
    end
    
    if h.aniso > 1
    
%        set(handles.text_theta1, 'String', h.str.s284)
%        set(handles.text_theta2, 'String', h.str.s284)
        set(handles.popup_model11t,'Visible','on')
        set(handles.model12t,'Visible','on')
        set(handles.model13t,'Visible','on')
        set(handles.model14t,'Visible','on')
        set(handles.c1t,'Visible','on')
        set(handles.popup_model21t,'Visible','on')
        set(handles.model22t,'Visible','on')
        set(handles.model23t,'Visible','on')
        set(handles.model24t,'Visible','on')
        set(handles.c2t,'Visible','on')
        set(handles.text_theta1,'Visible','on')
        set(handles.text_theta2,'Visible','on')
        set(handles.text_fixe1_t,'Visible','on')
        set(handles.text_fixe2_t,'Visible','on')
        set(handles.checkbox_mod12t,'Visible','on')
        set(handles.checkbox_mod13t,'Visible','on')
        set(handles.checkbox_mod14t,'Visible','on')
        set(handles.checkbox_c1t,'Visible','on')
        set(handles.checkbox_mod22t,'Visible','on')
        set(handles.checkbox_mod23t,'Visible','on')
        set(handles.checkbox_mod24t,'Visible','on')
        set(handles.checkbox_c2t,'Visible','on')
        set(handles.text_pepite_th,'Visible','on')
        set(handles.edit_pepite_th,'Visible','on')
        set(handles.checkbox_pepite_th,'Visible','on')

        set(handles.popup_model11t,'Value',h.model_th(1,1)-1) % Nugget (effect (code 1) not included in list
        set(handles.model12t,'String',num2str(h.model_th(1,2)));
        set(handles.model13t,'String',num2str(h.model_th(1,3)));
        set(handles.model14t,'String',num2str(h.model_th(1,4)));
        set(handles.c1t,'String',num2str(h.c_th(1,1)));
        set(handles.edit_pepite_th,'String',num2str(h.nugget_th));
        if size(h.model_th,1)==2
            set(handles.popup_model21t,'Value',h.model_th(2,1)-1) % Nugget (effect (code 1) not included in list
            set(handles.model22t,'String',num2str(h.model_th(2,2)));
            set(handles.model23t,'String',num2str(h.model_th(2,3)));
            set(handles.model24t,'String',num2str(h.model_th(2,4)));
            set(handles.c2t,'String',num2str(h.c_th(2,1)));
        end

    else
        set(handles.popup_model11t,'Visible','off')
        set(handles.model12t,'Visible','off')
        set(handles.model13t,'Visible','off')
        set(handles.model14t,'Visible','off')
        set(handles.c1t,'Visible','off')
        set(handles.popup_model21t,'Visible','off')
        set(handles.model22t,'Visible','off')
        set(handles.model23t,'Visible','off')
        set(handles.model24t,'Visible','off')
        set(handles.c2t,'Visible','off')
        set(handles.text_theta1,'Visible','off')
        set(handles.text_theta2,'Visible','off')
        set(handles.text_fixe1_t,'Visible','off')
        set(handles.text_fixe2_t,'Visible','off')
        set(handles.checkbox_mod12t,'Visible','off')
        set(handles.checkbox_mod13t,'Visible','off')
        set(handles.checkbox_mod14t,'Visible','off')
        set(handles.checkbox_c1t,'Visible','off')
        set(handles.checkbox_mod22t,'Visible','off')
        set(handles.checkbox_mod23t,'Visible','off')
        set(handles.checkbox_mod24t,'Visible','off')
        set(handles.checkbox_c2t,'Visible','off')
        set(handles.text_pepite_th,'Visible','off')
        set(handles.edit_pepite_th,'Visible','off')
        set(handles.checkbox_pepite_th,'Visible','off')
    end
    
else
%    set(handles.popup_model11x,'Visible','off')
%    set(handles.model12x,'Visible','off')
%    set(handles.model13x,'Visible','off')
%    set(handles.model14x,'Visible','off')
%    set(handles.c1x,'Visible','off')
%    set(handles.popup_model21x,'Visible','off')
%    set(handles.model22x,'Visible','off')
%    set(handles.model23x,'Visible','off')
%    set(handles.model24x,'Visible','off')
%    set(handles.c2x,'Visible','off')
%    set(handles.text_xi1,'Visible','off')
%    set(handles.text_xi2,'Visible','off')
%    set(handles.text_fixe1_x,'Visible','off')
%    set(handles.text_fixe1_t,'Visible','off')
%    set(handles.checkbox_mod12x,'Visible','off')
%    set(handles.checkbox_mod13x,'Visible','off')
%    set(handles.checkbox_mod14x,'Visible','off')
%    set(handles.checkbox_c1x,'Visible','off')
%    set(handles.checkbox_mod22x,'Visible','off')
%    set(handles.checkbox_mod23x,'Visible','off')
%    set(handles.checkbox_mod24x,'Visible','off')
%    set(handles.checkbox_c2x,'Visible','off')
    set(handles.text_pepite_xi,'Visible','off')
    set(handles.edit_pepite_xi,'Visible','off')
    set(handles.checkbox_pepite_xi,'Visible','off')

    set(handles.text_lenteur1, 'String', h.str.s95)
    set(handles.text_lenteur2, 'String', h.str.s95)
end


function checkbox_mod14_Callback(hObject, eventdata, handles)
h = getappdata(handles.fig_bh_FitCov, 'h_GUIdata');
h.model_code(1,2)=get(hObject,'Value');
setappdata(handles.fig_bh_FitCov, 'h_GUIdata', h);


function checkbox_mod12_Callback(hObject, eventdata, handles)
h = getappdata(handles.fig_bh_FitCov, 'h_GUIdata');
h.model_code(1,3)=get(hObject,'Value');
setappdata(handles.fig_bh_FitCov, 'h_GUIdata', h);


function checkbox_mod15_Callback(hObject, eventdata, handles)
h = getappdata(handles.fig_bh_FitCov, 'h_GUIdata');
h.model_code(1,4)=get(hObject,'Value');
setappdata(handles.fig_bh_FitCov, 'h_GUIdata', h);


function checkbox_c1_Callback(hObject, eventdata, handles)
h = getappdata(handles.fig_bh_FitCov, 'h_GUIdata');
h.c_code(1,1)=get(hObject,'Value');
setappdata(handles.fig_bh_FitCov, 'h_GUIdata', h);


function checkbox_pepite_t_Callback(hObject, eventdata, handles)
h = getappdata(handles.fig_bh_FitCov, 'h_GUIdata');
h.nugget_t_code=get(hObject,'Value');
setappdata(handles.fig_bh_FitCov, 'h_GUIdata', h);


function SaveMenuItem_Callback(hObject, eventdata, handles)
h = getappdata(handles.fig_bh_FitCov, 'h_GUIdata');

covar.model = h.model;
covar.c = h.c;
covar.nugget_t = h.nugget_t;
covar.nugget_l = h.nugget_l;
covar.use_c0 = get(handles.checkbox_c0,'Value');
covar.aniso = h.aniso;

load(h.db_file,'model3d')
if get(handles.popupmenu_type_data,'Value')==1
    covar.model_xi = h.model_xi;
    covar.c_xi = h.c_xi;
    covar.model_th = h.model_th;
    covar.c_th = h.c_th;
    covar.nugget_xi = h.nugget_xi;
    covar.nugget_th = h.nugget_th;
    
    h.model.tt_covar = covar;
    model3d(h.no_model3d).tt_covar = covar; %#ok<NASGU>
else
    h.model3d.amp_covar = covar;
    model3d(h.no_model3d).amp_covar = covar; %#ok<NASGU>
end
save(h.db_file,'model3d','-append')
h.saved = true;
setappdata(handles.fig_bh_FitCov, 'h_GUIdata', h);


function afi_Callback(hObject, eventdata, handles)
val = str2double(get(hObject,'string'));
if isnan(val)
    h = getappdata(handles.fig_bh_FitCov, 'h_GUIdata');
    errordlg(h.str.s54, h.str.s45,'modal')
	return
end
h = getappdata(handles.fig_bh_FitCov, 'h_GUIdata');
h.afi = 1/val;
setappdata(handles.fig_bh_FitCov, 'h_GUIdata', h);
CalculCovExperimentale(handles)
updateFigures(handles)


function afi_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


function HelpMenu_Callback(hObject, eventdata, handles)


function lclasHelpMenuItem_Callback(hObject, eventdata, handles)
h = getappdata(handles.fig_bh_FitCov, 'h_GUIdata');
helpdlg(h.str.s114, h.str.s102);


function typeModelHelpMenuItem_Callback(hObject, eventdata, handles)
h = getappdata(handles.fig_bh_FitCov, 'h_GUIdata');
helpdlg(h.str.s116, h.str.s117);


function afiHelpMenuItem_Callback(hObject, eventdata, handles)
h = getappdata(handles.fig_bh_FitCov, 'h_GUIdata');
helpdlg(h.str.s115, h.str.s103);


 function checkbox_c0_Callback(hObject, eventdata, handles)
h = getappdata(handles.fig_bh_FitCov, 'h_GUIdata');
if get(hObject,'Value') == 1
    h.c0 = h.data(:,2).^2;
else
    h.c0 = [];
end
setappdata(handles.fig_bh_FitCov, 'h_GUIdata', h);
updateFigures(handles)


function rais_droits_Callback(hObject, eventdata, handles)


function rais_courbes_Callback(hObject, eventdata, handles)


function GetRaisCourbes(handles)
[fichier, rep] = uigetfile('*.mat','Fichier contenant la matrice L');
if fichier==0
    set(handles.rais_droits,'Value',1);
    return
end
tmp=load([rep,fichier]);
h = getappdata(handles.fig_bh_FitCov, 'h_GUIdata');
no_trace = h.data(:,3);
nn = 0;
ind = zeros(size(tmp.no_trace));
for n=1:length(no_trace)
  ii=find(no_trace(n)==tmp.no_trace);
  if ~isempty(ii)
	nn = nn+1;
	ind(nn) = ii;
  end
end
%ind=ind(find(ind));
ind=ind(ind~=0);
h.L = sparse(tmp.L(ind,:));

h.gridx = tmp.gridx;
h.gridy = tmp.gridy;
h.gridz = tmp.gridz;
h.xmin = h.gridx(1);
h.xmax = h.gridx(length(h.gridx));
h.dx = h.gridx(2)-h.gridx(1);
h.ymin = h.gridy(1);
h.ymax = h.gridy(length(h.gridy));
h.dy = h.gridy(2)-h.gridy(1);
h.zmin = h.gridz(1);
h.zmax = h.gridz(length(h.gridz));
h.dz = h.gridz(2)-h.gridz(1);
n=length(h.gridz);
m=length(h.gridx);
l=length(h.gridy);
h.x=grille3d(min(h.gridz),h.gridz(2)-h.gridz(1),n,...
             min(h.gridy),h.gridy(2)-h.gridy(1),l,...
             min(h.gridx),h.gridx(2)-h.gridx(1),m);
set(handles.xmin,'String',num2str(h.xmin));
set(handles.xmax,'String',num2str(h.xmax));
set(handles.dx,'String',num2str(h.dx));
set(handles.ymin,'String',num2str(h.ymin));
set(handles.ymax,'String',num2str(h.ymax));
set(handles.dy,'String',num2str(h.dy));
set(handles.zmin,'String',num2str(h.zmin));
set(handles.zmax,'String',num2str(h.zmax));
set(handles.dz,'String',num2str(h.dz));
set(handles.text_nombre_cellules,'String',[num2str(n*m), ' ',h.str.s205])
setappdata(handles.fig_bh_FitCov, 'h_GUIdata', h);


function Boutons_Rais_SelectionChangeFcn(hObject, eventdata, handles)
switch get(hObject,'Tag')   % Get Tag of selected object
    case 'rais_droits'
        CalculRaisDroits(handles);
    case 'rais_courbes'
        GetRaisCourbes(handles);
end
h = getappdata(handles.fig_bh_FitCov, 'h_GUIdata');
if isempty(h.data), return, end
CalculCovExperimentale(handles)
updateFigures(handles)


function popup_model11_Callback(hObject, eventdata, handles)
val = get(hObject,'Value');
h = getappdata(handles.fig_bh_FitCov, 'h_GUIdata');
h.model(1,1) = val + 1; % Nugget (effect (code 1) not included in list
h.saved = false;
setappdata(handles.fig_bh_FitCov, 'h_GUIdata', h);
if get(handles.checkbox_calc_auto,'Value')==1, updateFigures(handles), end


function popup_model11_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


function set_String_locale(handles, str)
set(handles.uipanel1,                'Title',  str.s23)
set(handles.uipanel2,                'Title',  str.s92)
set(handles.panel_covar2d,                'Title',  str.s93)
set(handles.uipanel5,                'Title',  str.s104)
set(handles.checkbox_c0,             'String', str.s107)
set(handles.popup_model11,           'String', str.s94)
set(handles.popup_model21,           'String', str.s94)
%set(handles.popup_model11x,          'String', str.s94)
%set(handles.popup_model21x,          'String', str.s94)
%set(handles.popup_model11t,          'String', str.s94)
%set(handles.popup_model21t,          'String', str.s94)
set(handles.text_lenteur1,           'String', str.s95)
set(handles.text_lenteur2,           'String', str.s95)
set(handles.text_fixe1_s,            'String', str.s96)
set(handles.text_fixe2_s,            'String', str.s96)
%set(handles.text_fixe1_x,            'String', str.s96)
%set(handles.text_fixe2_x,            'String', str.s96)
%set(handles.text_fixe1_t,            'String', str.s96)
%set(handles.text_fixe2_t,            'String', str.s96)
set(handles.model14l,                'String', str.s97)
set(handles.model14l,                'String', str.s98)
set(handles.model15l,                'String', str.s99)
set(handles.c1l,                     'String', str.s100)

set(handles.model12l,                'String', str.s97)
set(handles.model13l,                'String', str.s977)
set(handles.model14l,                'String', str.s98)
set(handles.model15l,                'String', str.s99)
set(handles.model16l,                'String', str.s991)
set(handles.model17l,                'String', str.s992)
set(handles.c2l,                     'String', str.s100)
set(handles.model22l,                'String', str.s97)
set(handles.model23l,                'String', str.s977)
set(handles.model24l,                'String', str.s98)
set(handles.model25l,                'String', str.s99)
set(handles.model26l,                'String', str.s991)
set(handles.model27l,                'String', str.s992)
set(handles.c2l,                     'String', str.s100)
set(handles.text_pepite_t,           'String', str.s101)
set(handles.text_pepite_l,           'String', str.s112)
set(handles.c4l,                     'String', str.s102)
set(handles.text23,                  'String', str.s103)
set(handles.text20,                  'String', str.s105)
set(handles.pushbutton_ajuster,      'String', str.s106)
set(handles.FileMenu,                'Label',  str.s25)
set(handles.OpenMenuItem,            'Label',  str.s118)
set(handles.SaveMenuItem,            'Label',  str.s29)
set(handles.CloseMenuItem,           'Label',  str.s31)
set(handles.HelpMenu,                'Label',  str.s32)
set(handles.typeModelHelpMenuItem,   'Label',  str.s93)
set(handles.lclasHelpMenuItem,       'Label',  str.s102)
set(handles.afiHelpMenuItem,         'Label',  str.s103)
set(handles.EditMenuItem,            'Label',  str.s122)
set(handles.RemoveES2MenuItem,       'Label',  str.s57)
set(handles.popupmenu_type_data,     'String', str.s156)
set(handles.uipanel_rais,            'Title',  str.s153)
set(handles.text5,                   'String', lower(str.s162))
set(handles.uipanel_param_struct1,   'Title',  str.s164)
set(handles.uipanel_param_struct2,   'Title',  str.s164)
set(handles.checkbox_lim_vapp,       'String', str.s286)
%set(handles.checkbox_aniso,          'String', str.s287)
%set(handles.checkbox_aniso_rot,      'String', str.s288)
set(handles.text_pepite_xi,          'String', str.s289)
set(handles.text_pepite_th,          'String', str.s290)
set(handles.checkbox_calc_auto,      'String', str.s291)
set(handles.pushbutton_calcul_covar, 'String', str.s292)


function popupmenu_type_data_Callback(hObject, eventdata, handles)
h = getappdata(handles.fig_bh_FitCov, 'h_GUIdata');
if ~isfield(h, 'model3d')
    return
end
selected_mogs = get(handles.listbox_mogs,'Value');
if get(hObject,'Value')==1
    [h.data,h.ind] = getPanneauData(h.model3d, h.db_file,'tt',selected_mogs);
    if isfield(h.model3d,'tt_covar')
		if ~isempty( h.model3d.tt_covar )
			h.model = h.model3d.tt_covar.model;
			h.c = h.model3d.tt_covar.c;
			h.nugget_t = h.model3d.tt_covar.nugget_t;
			h.nugget_l = h.model3d.tt_covar.nugget_l;
			set(handles.checkbox_c0,'Value',h.model3d.tt_covar.use_c0)
            if isfield( h.model3d.tt_covar, 'aniso' )
                h.aniso = h.model3d.tt_covar.aniso;
                if h.aniso > 0
                    h.model_xi = h.model3d.tt_covar.model_xi;
                    h.c_xi = h.model3d.tt_covar.c_xi;
                    h.nugget_xi = h.model3d.tt_covar.nugget_xi;
                    if h.aniso > 1
                        h.model_th = h.model3d.tt_covar.model_th;
                        h.c_th = h.model3d.tt_covar.c_th;
                        h.nugget_th = h.model3d.tt_covar.nugget_th;
                    end
                end
            end
		end
    end
    
else
    
    if h.aniso > 0
        warndlg('Elliptic anisotropy not available for amplitude data');
        set(hObject,'Value',1)
        return
    end
    
	if get(hObject,'Value')==2
		[h.data,h.ind] = getPanneauData(h.model3d, h.db_file,'amp',selected_mogs);
	elseif get(hObject,'Value')==3
		[h.data,h.ind] = getPanneauData(h.model3d, h.db_file,'fce',selected_mogs);
	elseif get(hObject,'Value')==4
		[h.data,h.ind] = getPanneauData(h.model3d, h.db_file,'hyb',selected_mogs);
	end
    if isfield(h.model3d,'amp_covar')
		if ~isempty( h.model3d.amp_covar )
			h.model = h.model3d.amp_covar.model;
			h.c = h.model3d.amp_covar.c;
			h.nugget_t = h.model3d.amp_covar.nugget_t;
			h.nugget_l = h.model3d.amp_covar.nugget_l;
			set(handles.checkbox_c0,'Value',0)
		end
    end
end
h.L = [];
setappdata(handles.fig_bh_FitCov, 'h_GUIdata', h);
if isempty(h.data)
	warndlg(h.str.s131,h.str.s45)
% 	if get(hObject,'Value')==1
% 		set(hObject,'Value',2)
% 	else
% 		set(hObject,'Value',1)
% 	end
	return
end

set(handles.nombre_rais,'String',[num2str(size(h.data,1)),' ',...
					lower(h.str.s64)])
updateCovAffichage(handles)
CalculCovExperimentale(handles)
updateFigures(handles)


function popupmenu_type_data_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


function popupmenu_Ldc_Callback(hObject, eventdata, handles)
h = getappdata(handles.fig_bh_FitCov, 'h_GUIdata');
if isempty(h.data), return, end
no_Ldc = get(handles.popupmenu_Ldc,'Value');
if no_Ldc>1
%    ind = zeros(1,length(h.model3d.inv_res(no_Ldc-1).tomo.no_trace));
    ind = [];
    for n=1:length(h.data(:,3)) %#ok<ALIGN>
%        ind(n) = findnear(h.model3d.inv_res(no_Ldc-1).tomo.no_trace(n), h.data(:,3));
        ii = find( h.model3d.inv_res(no_Ldc-1).tomo.no_trace(:)==h.data(n,3) );
% 		if isempty(ii)
%             uiwait(warndlg([h.str.s186{1}, num2str(h.data(n,3)), h.str.s186{2}]))
% 		else
    if ~isempty(ii)
			   ind = [ind ii];
		end
	end
	if isempty(ind)
		uiwait(warndlg(h.str.s188{1}))
		return
	end
    h.L = h.model3d.inv_res(no_Ldc-1).tomo.L(ind,:);
    h.gridx = h.model3d.inv_res(no_Ldc-1).tomo.x;
    h.gridy = h.model3d.inv_res(no_Ldc-1).tomo.y;
    h.gridz = h.model3d.inv_res(no_Ldc-1).tomo.z;
    h.xmin = h.gridx(1);
    h.xmax = h.gridx(length(h.gridx));
    h.dx = h.gridx(2)-h.gridx(1);
    h.ymin = h.gridy(1);
    h.ymax = h.gridy(length(h.gridy));
    h.dy = h.gridy(2)-h.gridy(1);
    h.zmin = h.gridz(1);
    h.zmax = h.gridz(length(h.gridz));
    h.dz = h.gridz(2)-h.gridz(1);
    n=length(h.gridz);
    m=length(h.gridx);
    l=length(h.gridy);
    h.x=grille3d(min(h.gridz),h.gridz(2)-h.gridz(1),n,...
                 min(h.gridy),h.gridy(2)-h.gridy(1),l,...
                 min(h.gridx),h.gridx(2)-h.gridx(1),m);
    set(handles.xmin,'String',num2str(h.xmin));
    set(handles.xmax,'String',num2str(h.xmax));
    set(handles.dx,'String',num2str(h.dx));
    set(handles.ymin,'String',num2str(h.ymin));
    set(handles.ymax,'String',num2str(h.ymax));
    set(handles.dy,'String',num2str(h.dy));
    set(handles.zmin,'String',num2str(h.zmin));
    set(handles.zmax,'String',num2str(h.zmax));
    set(handles.dz,'String',num2str(h.dz));
    set(handles.text_nombre_cellules,'String',[num2str(n*m), ' ',h.str.s205])
    setappdata(handles.fig_bh_FitCov, 'h_GUIdata', h);
else
    CalculRaisDroits(handles);
end
CalculCovExperimentale(handles)
updateFigures(handles)


function popupmenu_Ldc_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


function popupmenu_no_structure_Callback(hObject, eventdata, handles)
no = get(hObject,'Value');
if no==1
	set(handles.uipanel_param_struct2,'Visible','off')
	set(handles.uipanel_param_struct1,'Visible','on')
else
	h=getappdata(handles.fig_bh_FitCov, 'h_GUIdata');
	if size(h.model,1)==1
		h.model = [h.model; 1 0 0 0];
		h.c = [h.c; 1];
		h.model_code=[h.model_code; 1 0 0 0];
		h.c_code=[h.c_code; 0];
		h.model_xi = [h.model_xi; 1 0 0 0];
		h.c_xi = [h.c_xi; 1];
		h.model_xi_code=[h.model_xi_code; 1 0 0 0];
		h.c_xi_code=[h.c_xi_code; 0];
		setappdata(handles.fig_bh_FitCov, 'h_GUIdata',h)
		updateCovAffichage(handles)
	end
	if size(h.model,1)>size(h.model_code,1)
		h.model_code=[h.model_code; 1 0 0 0];
		h.c_code=[h.c_code; 0];
		h.model_xi_code=[h.model_xi_code; 1 0 0 0];
		h.c_xi_code=[h.c_xi_code; 0];
		setappdata(handles.fig_bh_FitCov, 'h_GUIdata',h)
		updateCovAffichage(handles)
	end
	set(handles.uipanel_param_struct1,'Visible','off')
	set(handles.uipanel_param_struct2,'Visible','on')
end


function popupmenu_no_structure_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


function model22_Callback(hObject, eventdata, handles)
val = str2double(get(hObject,'string'));
if isnan(val)
    h = getappdata(handles.fig_bh_FitCov, 'h_GUIdata');
    errordlg(h.str.s54, h.str.s45,'modal')
	return
end
h = getappdata(handles.fig_bh_FitCov, 'h_GUIdata');
h.model(2,2) = val;
h.saved = false;
setappdata(handles.fig_bh_FitCov, 'h_GUIdata', h);
if get(handles.checkbox_calc_auto,'Value')==1, updateFigures(handles), end


function model22_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


function model23_Callback(hObject, eventdata, handles)
val = str2double(get(hObject,'string'));
if isnan(val)
    h = getappdata(handles.fig_bh_FitCov, 'h_GUIdata');
    errordlg(h.str.s54, h.str.s45,'modal')
	return
end
h = getappdata(handles.fig_bh_FitCov, 'h_GUIdata');
h.model(2,3) = val;
h.saved = false;
setappdata(handles.fig_bh_FitCov, 'h_GUIdata', h);
if get(handles.checkbox_calc_auto,'Value')==1, updateFigures(handles), end


function model23_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


function model24_Callback(hObject, eventdata, handles)
val = str2double(get(hObject,'string'));
if isnan(val)
    h = getappdata(handles.fig_bh_FitCov, 'h_GUIdata');
    errordlg(h.str.s54, h.str.s45,'modal')
	return
end
h = getappdata(handles.fig_bh_FitCov, 'h_GUIdata');
h.model(2,4) = val;
h.saved = false;
setappdata(handles.fig_bh_FitCov, 'h_GUIdata', h);
if get(handles.checkbox_calc_auto,'Value')==1, updateFigures(handles), end


function model24_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


function c2_Callback(hObject, eventdata, handles)
val = str2double(get(hObject,'string'));
if isnan(val)
    h = getappdata(handles.fig_bh_FitCov, 'h_GUIdata');
    errordlg(h.str.s54, h.str.s45,'modal')
	return
end
h = getappdata(handles.fig_bh_FitCov, 'h_GUIdata');
h.c(2,1) = val;
h.saved = false;
setappdata(handles.fig_bh_FitCov, 'h_GUIdata', h);
if get(handles.checkbox_calc_auto,'Value')==1, updateFigures(handles), end


function c2_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


function checkbox_mod22_Callback(hObject, eventdata, handles)
h = getappdata(handles.fig_bh_FitCov, 'h_GUIdata');
h.model_code(2,2)=get(hObject,'Value');
setappdata(handles.fig_bh_FitCov, 'h_GUIdata', h);


function checkbox_mod23_Callback(hObject, eventdata, handles)
h = getappdata(handles.fig_bh_FitCov, 'h_GUIdata');
h.model_code(2,3)=get(hObject,'Value');
setappdata(handles.fig_bh_FitCov, 'h_GUIdata', h);


function checkbox_mod24_Callback(hObject, eventdata, handles)
h = getappdata(handles.fig_bh_FitCov, 'h_GUIdata');
h.model_code(2,4)=get(hObject,'Value');
setappdata(handles.fig_bh_FitCov, 'h_GUIdata', h);


function checkbox_c2_Callback(hObject, eventdata, handles)
h = getappdata(handles.fig_bh_FitCov, 'h_GUIdata');
h.c_code(2,1)=get(hObject,'Value');
setappdata(handles.fig_bh_FitCov, 'h_GUIdata', h);


function popup_model21_Callback(hObject, eventdata, handles)
val = get(hObject,'Value');
h = getappdata(handles.fig_bh_FitCov, 'h_GUIdata');
h.model(2,1) = val + 1; % Nugget (effect (code 1) not included in list
h.saved = false;
setappdata(handles.fig_bh_FitCov, 'h_GUIdata', h);
if get(handles.checkbox_calc_auto,'Value')==1, updateFigures(handles), end


function popup_model21_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


function edit_pepite_l_Callback(hObject, eventdata, handles)
val = str2double(get(hObject,'string'));
if isnan(val)
    h = getappdata(handles.fig_bh_FitCov, 'h_GUIdata');
    errordlg(h.str.s54, h.str.s45,'modal')
	return
end
h = getappdata(handles.fig_bh_FitCov, 'h_GUIdata');
h.nugget_l = val;
h.saved = false;
setappdata(handles.fig_bh_FitCov, 'h_GUIdata', h);
if get(handles.checkbox_calc_auto,'Value')==1, updateFigures(handles), end


function edit_pepite_l_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


function checkbox_pepite_l_Callback(hObject, eventdata, handles)
h = getappdata(handles.fig_bh_FitCov, 'h_GUIdata');
h.nugget_l_code=get(hObject,'Value');
h.saved = false;
setappdata(handles.fig_bh_FitCov, 'h_GUIdata', h);


function EditMenuItem_Callback(hObject, eventdata, handles)


function RemoveES2MenuItem_Callback(hObject, eventdata, handles)
h = getappdata(handles.fig_bh_FitCov, 'h_GUIdata');
h.model = h.model(1,:);
h.c = h.c(1);
h.model_code = h.model_code(1,:);
h.c_code = h.c_code(1);

h.model_xi = h.model_xi(1,:);
h.c_xi = h.c_xi(1);
h.model_xi_code = h.model_xi_code(1,:);
h.c_xi_code = h.c_xi_code(1);

h.model_th = h.model_th(1,:);
h.c_th = h.c_th(1);
h.model_th_code = h.model_th_code(1,:);
h.c_th_code = h.c_th_code(1);

set(handles.popupmenu_no_structure,'Value',1)
set(handles.uipanel_param_struct1,'Visible','on')
set(handles.uipanel_param_struct2,'Visible','off')

h.saved = false;
setappdata(handles.fig_bh_FitCov, 'h_GUIdata', h);
updateFigures(handles)


function checkbox_lim_vapp_Callback(hObject, eventdata, handles)
update_vlim(handles)


function edit_lim_vapp_Callback(hObject, eventdata, handles)
update_vlim(handles)


function edit_lim_vapp_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


function update_vlim(handles)
h = getappdata(handles.fig_bh_FitCov, 'h_GUIdata');
selected_mogs=get(handles.listbox_mogs,'Value');
if get(handles.popupmenu_type_data,'Value')==1
    if get(handles.checkbox_lim_vapp,'Value')==1
        vlim = str2double(get(handles.edit_lim_vapp,'String'));
        selected_mogs = 1:length(get(handles.listbox_mogs,'String'));
    else
        vlim = [];
    end
    [h.data,h.ind] = getPanneauData(h.model3d, h.db_file,'tt',selected_mogs, vlim);
end

no_Ldc = get(handles.popupmenu_Ldc,'Value');
if no_Ldc>1 % rais courbes
    ind = [];
    for n=1:length(h.data(:,3)) %#ok<ALIGN>
%        ind(n) = findnear(h.model3d.inv_res(no_Ldc-1).tomo.no_trace(n), h.data(:,3));
        ii = find( h.model3d.inv_res(no_Ldc-1).tomo.no_trace(:)==h.data(n,3) );
% 		if isempty(ii)
%             uiwait(warndlg([h.str.s186{1}, num2str(h.data(n,3)), h.str.s186{2}]))
% 		else
        if ~isempty(ii)
			ind = [ind ii];
		end
	end
	if isempty(ind)
		uiwait(warndlg(h.str.s188{1}))
		return
	end
    h.L = h.model3d.inv_res(no_Ldc-1).tomo.L(ind,:);
else % rais droits
	h.L = [];
end

set(handles.nombre_rais,'String',[num2str(size(h.data,1)),' ',lower(h.str.s64)])
setappdata(handles.fig_bh_FitCov, 'h_GUIdata', h);
CalculCovExperimentale(handles)
updateFigures(handles)


function checkbox_aniso_Callback(hObject, eventdata, handles)
h = getappdata(handles.fig_bh_FitCov, 'h_GUIdata');
h.aniso = get(hObject,'Value');
setappdata(handles.fig_bh_FitCov, 'h_GUIdata', h);
if h.aniso==1
    if get(handles.popupmenu_type_data,'Value')~=1
        warndlg('Elliptic anisotropy not available for amplitude data');
        set(hObject,'Value',0)
        h.aniso = 0;
        setappdata(handles.fig_bh_FitCov, 'h_GUIdata', h);
        return
    end
end
updateCovAffichage(handles)
CalculRaisDroits(handles)
CalculCovExperimentale(handles)

if ~isempty(h.data)
  updateFigures(handles)
end


function model13x_Callback(hObject, eventdata, handles)
val = str2double(get(hObject,'string'));
if isnan(val)
    h = getappdata(handles.fig_bh_FitCov, 'h_GUIdata');
    errordlg(h.str.s54, h.str.s45,'modal')
	return
end
h = getappdata(handles.fig_bh_FitCov, 'h_GUIdata');
h.model_xi(1,3) = val;
h.saved = false;
setappdata(handles.fig_bh_FitCov, 'h_GUIdata', h);
if get(handles.checkbox_calc_auto,'Value')==1, updateFigures(handles), end


function model13x_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


function model12x_Callback(hObject, eventdata, handles)
val = str2double(get(hObject,'string'));
if isnan(val)
    h = getappdata(handles.fig_bh_FitCov, 'h_GUIdata');
    errordlg(h.str.s54, h.str.s45,'modal')
	return
end
h = getappdata(handles.fig_bh_FitCov, 'h_GUIdata');
h.model_xi(1,2) = val;
h.saved = false;
setappdata(handles.fig_bh_FitCov, 'h_GUIdata', h);
if get(handles.checkbox_calc_auto,'Value')==1, updateFigures(handles), end


function model12x_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


function model14x_Callback(hObject, eventdata, handles)
val = str2double(get(hObject,'string'));
if isnan(val)
    h = getappdata(handles.fig_bh_FitCov, 'h_GUIdata');
    errordlg(h.str.s54, h.str.s45,'modal')
	return
end
h = getappdata(handles.fig_bh_FitCov, 'h_GUIdata');
h.model_xi(1,4) = val;
h.saved = false;
setappdata(handles.fig_bh_FitCov, 'h_GUIdata', h);
if get(handles.checkbox_calc_auto,'Value')==1, updateFigures(handles), end


function model14x_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


function c1x_Callback(hObject, eventdata, handles)
val = str2double(get(hObject,'string'));
if isnan(val)
    h = getappdata(handles.fig_bh_FitCov, 'h_GUIdata');
    errordlg(h.str.s54, h.str.s45,'modal')
	return
end
h = getappdata(handles.fig_bh_FitCov, 'h_GUIdata');
h.c_xi(1,1) = val;
h.saved = false;
setappdata(handles.fig_bh_FitCov, 'h_GUIdata', h);
if get(handles.checkbox_calc_auto,'Value')==1, updateFigures(handles), end


function c1x_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


function model23x_Callback(hObject, eventdata, handles)
val = str2double(get(hObject,'string'));
if isnan(val)
    h = getappdata(handles.fig_bh_FitCov, 'h_GUIdata');
    errordlg(h.str.s54, h.str.s45,'modal')
	return
end
h = getappdata(handles.fig_bh_FitCov, 'h_GUIdata');
h.model_xi(2,3) = val;
h.saved = false;
setappdata(handles.fig_bh_FitCov, 'h_GUIdata', h);
if get(handles.checkbox_calc_auto,'Value')==1, updateFigures(handles), end


function model23x_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


function model22x_Callback(hObject, eventdata, handles)
val = str2double(get(hObject,'string'));
if isnan(val)
    h = getappdata(handles.fig_bh_FitCov, 'h_GUIdata');
    errordlg(h.str.s54, h.str.s45,'modal')
	return
end
h = getappdata(handles.fig_bh_FitCov, 'h_GUIdata');
h.model_xi(2,2) = val;
h.saved = false;
setappdata(handles.fig_bh_FitCov, 'h_GUIdata', h);
if get(handles.checkbox_calc_auto,'Value')==1, updateFigures(handles), end


function model22x_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


function model24x_Callback(hObject, eventdata, handles)
val = str2double(get(hObject,'string'));
if isnan(val)
    h = getappdata(handles.fig_bh_FitCov, 'h_GUIdata');
    errordlg(h.str.s54, h.str.s45,'modal')
	return
end
h = getappdata(handles.fig_bh_FitCov, 'h_GUIdata');
h.model_xi(2,4) = val;
h.saved = false;
setappdata(handles.fig_bh_FitCov, 'h_GUIdata', h);
if get(handles.checkbox_calc_auto,'Value')==1, updateFigures(handles), end


function model24x_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


function c2x_Callback(hObject, eventdata, handles)
val = str2double(get(hObject,'string'));
if isnan(val)
    h = getappdata(handles.fig_bh_FitCov, 'h_GUIdata');
    errordlg(h.str.s54, h.str.s45,'modal')
	return
end
h = getappdata(handles.fig_bh_FitCov, 'h_GUIdata');
h.c_xi(2,1) = val;
h.saved = false;
setappdata(handles.fig_bh_FitCov, 'h_GUIdata', h);
if get(handles.checkbox_calc_auto,'Value')==1, updateFigures(handles), end


function c2x_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


function model13t_Callback(hObject, eventdata, handles)
val = str2double(get(hObject,'string'));
if isnan(val)
    h = getappdata(handles.fig_bh_FitCov, 'h_GUIdata');
    errordlg(h.str.s54, h.str.s45,'modal')
	return
end
h = getappdata(handles.fig_bh_FitCov, 'h_GUIdata');
h.model_th(1,3) = val;
h.saved = false;
setappdata(handles.fig_bh_FitCov, 'h_GUIdata', h);
if get(handles.checkbox_calc_auto,'Value')==1, updateFigures(handles), end


function model13t_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


function model12t_Callback(hObject, eventdata, handles)
val = str2double(get(hObject,'string'));
if isnan(val)
    h = getappdata(handles.fig_bh_FitCov, 'h_GUIdata');
    errordlg(h.str.s54, h.str.s45,'modal')
	return
end
h = getappdata(handles.fig_bh_FitCov, 'h_GUIdata');
h.model_th(1,2) = val;
h.saved = false;
setappdata(handles.fig_bh_FitCov, 'h_GUIdata', h);
if get(handles.checkbox_calc_auto,'Value')==1, updateFigures(handles), end


function model12t_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


function model14t_Callback(hObject, eventdata, handles)
val = str2double(get(hObject,'string'));
if isnan(val)
    h = getappdata(handles.fig_bh_FitCov, 'h_GUIdata');
    errordlg(h.str.s54, h.str.s45,'modal')
	return
end
h = getappdata(handles.fig_bh_FitCov, 'h_GUIdata');
h.model_th(1,4) = val;
h.saved = false;
setappdata(handles.fig_bh_FitCov, 'h_GUIdata', h);
if get(handles.checkbox_calc_auto,'Value')==1, updateFigures(handles), end


function model14t_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


function c1t_Callback(hObject, eventdata, handles)
val = str2double(get(hObject,'string'));
if isnan(val)
    h = getappdata(handles.fig_bh_FitCov, 'h_GUIdata');
    errordlg(h.str.s54, h.str.s45,'modal')
	return
end
h = getappdata(handles.fig_bh_FitCov, 'h_GUIdata');
h.c_th(1,1) = val;
h.saved = false;
setappdata(handles.fig_bh_FitCov, 'h_GUIdata', h);
if get(handles.checkbox_calc_auto,'Value')==1, updateFigures(handles), end


function c1t_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


function model23t_Callback(hObject, eventdata, handles)
val = str2double(get(hObject,'string'));
if isnan(val)
    h = getappdata(handles.fig_bh_FitCov, 'h_GUIdata');
    errordlg(h.str.s54, h.str.s45,'modal')
	return
end
h = getappdata(handles.fig_bh_FitCov, 'h_GUIdata');
h.model_th(2,3) = val;
h.saved = false;
setappdata(handles.fig_bh_FitCov, 'h_GUIdata', h);
if get(handles.checkbox_calc_auto,'Value')==1, updateFigures(handles), end


function model23t_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


function model22t_Callback(hObject, eventdata, handles)
val = str2double(get(hObject,'string'));
if isnan(val)
    h = getappdata(handles.fig_bh_FitCov, 'h_GUIdata');
    errordlg(h.str.s54, h.str.s45,'modal')
	return
end
h = getappdata(handles.fig_bh_FitCov, 'h_GUIdata');
h.model_th(2,2) = val;
h.saved = false;
setappdata(handles.fig_bh_FitCov, 'h_GUIdata', h);
if get(handles.checkbox_calc_auto,'Value')==1, updateFigures(handles), end


function model22t_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


function model24t_Callback(hObject, eventdata, handles)
val = str2double(get(hObject,'string'));
if isnan(val)
    h = getappdata(handles.fig_bh_FitCov, 'h_GUIdata');
    errordlg(h.str.s54, h.str.s45,'modal')
	return
end
h = getappdata(handles.fig_bh_FitCov, 'h_GUIdata');
h.model_th(2,4) = val;
h.saved = false;
setappdata(handles.fig_bh_FitCov, 'h_GUIdata', h);
if get(handles.checkbox_calc_auto,'Value')==1, updateFigures(handles), end


function model24t_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


function c2t_Callback(hObject, eventdata, handles)
val = str2double(get(hObject,'string'));
if isnan(val)
    h = getappdata(handles.fig_bh_FitCov, 'h_GUIdata');
    errordlg(h.str.s54, h.str.s45,'modal')
	return
end
h = getappdata(handles.fig_bh_FitCov, 'h_GUIdata');
h.c_th(2,1) = val;
h.saved = false;
setappdata(handles.fig_bh_FitCov, 'h_GUIdata', h);
if get(handles.checkbox_calc_auto,'Value')==1, updateFigures(handles), end


function c2t_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


function popup_model11x_Callback(hObject, eventdata, handles)
val = get(hObject,'Value');
h = getappdata(handles.fig_bh_FitCov, 'h_GUIdata');
h.model_xi(1,1) = val + 1; % Nugget (effect (code 1) not included in list
h.saved = false;
setappdata(handles.fig_bh_FitCov, 'h_GUIdata', h);
if get(handles.checkbox_calc_auto,'Value')==1, updateFigures(handles), end


function popup_model11x_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


function popup_model21x_Callback(hObject, eventdata, handles)
val = get(hObject,'Value');
h = getappdata(handles.fig_bh_FitCov, 'h_GUIdata');
h.model_xi(2,1) = val + 1; % Nugget (effect (code 1) not included in list
h.saved = false;
setappdata(handles.fig_bh_FitCov, 'h_GUIdata', h);
if get(handles.checkbox_calc_auto,'Value')==1, updateFigures(handles), end


function popup_model21x_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


function popup_model11t_Callback(hObject, eventdata, handles)
val = get(hObject,'Value');
h = getappdata(handles.fig_bh_FitCov, 'h_GUIdata');
h.model_th(1,1) = val + 1; % Nugget (effect (code 1) not included in list
h.saved = false;
setappdata(handles.fig_bh_FitCov, 'h_GUIdata', h);
if get(handles.checkbox_calc_auto,'Value')==1, updateFigures(handles), end


function popup_model11t_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


function popup_model21t_Callback(hObject, eventdata, handles)
val = get(hObject,'Value');
h = getappdata(handles.fig_bh_FitCov, 'h_GUIdata');
h.model_th(2,1) = val + 1; % Nugget (effect (code 1) not included in list
h.saved = false;
setappdata(handles.fig_bh_FitCov, 'h_GUIdata', h);
if get(handles.checkbox_calc_auto,'Value')==1, updateFigures(handles), end


function popup_model21t_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


function checkbox_mod13x_Callback(hObject, eventdata, handles)
h = getappdata(handles.fig_bh_FitCov, 'h_GUIdata');
h.model_xi_code(1,3)=get(hObject,'Value');
setappdata(handles.fig_bh_FitCov, 'h_GUIdata', h);


function checkbox_mod12x_Callback(hObject, eventdata, handles)
h = getappdata(handles.fig_bh_FitCov, 'h_GUIdata');
h.model_xi_code(1,2)=get(hObject,'Value');
setappdata(handles.fig_bh_FitCov, 'h_GUIdata', h);


function checkbox_mod14x_Callback(hObject, eventdata, handles)
h = getappdata(handles.fig_bh_FitCov, 'h_GUIdata');
h.model_xi_code(1,4)=get(hObject,'Value');
setappdata(handles.fig_bh_FitCov, 'h_GUIdata', h);


function checkbox_c1x_Callback(hObject, eventdata, handles)
h = getappdata(handles.fig_bh_FitCov, 'h_GUIdata');
h.c_xi_code(1,1)=get(hObject,'Value');
setappdata(handles.fig_bh_FitCov, 'h_GUIdata', h);


function checkbox_mod23x_Callback(hObject, eventdata, handles)
h = getappdata(handles.fig_bh_FitCov, 'h_GUIdata');
h.model_xi_code(2,3)=get(hObject,'Value');
setappdata(handles.fig_bh_FitCov, 'h_GUIdata', h);


function checkbox_mod22x_Callback(hObject, eventdata, handles)
h = getappdata(handles.fig_bh_FitCov, 'h_GUIdata');
h.model_xi_code(2,2)=get(hObject,'Value');
setappdata(handles.fig_bh_FitCov, 'h_GUIdata', h);


function checkbox_mod24x_Callback(hObject, eventdata, handles)
h = getappdata(handles.fig_bh_FitCov, 'h_GUIdata');
h.model_xi_code(2,4)=get(hObject,'Value');
setappdata(handles.fig_bh_FitCov, 'h_GUIdata', h);


function checkbox_c2x_Callback(hObject, eventdata, handles)
h = getappdata(handles.fig_bh_FitCov, 'h_GUIdata');
h.c_xi_code(2,1)=get(hObject,'Value');
setappdata(handles.fig_bh_FitCov, 'h_GUIdata', h);


function checkbox_mod13t_Callback(hObject, eventdata, handles)
h = getappdata(handles.fig_bh_FitCov, 'h_GUIdata');
h.model_th_code(1,3)=get(hObject,'Value');
setappdata(handles.fig_bh_FitCov, 'h_GUIdata', h);


function checkbox_mod12t_Callback(hObject, eventdata, handles)
h = getappdata(handles.fig_bh_FitCov, 'h_GUIdata');
h.model_th_code(1,2)=get(hObject,'Value');
setappdata(handles.fig_bh_FitCov, 'h_GUIdata', h);


function checkbox_mod14t_Callback(hObject, eventdata, handles)
h = getappdata(handles.fig_bh_FitCov, 'h_GUIdata');
h.model_th_code(1,4)=get(hObject,'Value');
setappdata(handles.fig_bh_FitCov, 'h_GUIdata', h);


function checkbox_c1t_Callback(hObject, eventdata, handles)
h = getappdata(handles.fig_bh_FitCov, 'h_GUIdata');
h.c_th_code(1,1)=get(hObject,'Value');
setappdata(handles.fig_bh_FitCov, 'h_GUIdata', h);


function checkbox_mod23t_Callback(hObject, eventdata, handles)
h = getappdata(handles.fig_bh_FitCov, 'h_GUIdata');
h.model_th_code(2,3)=get(hObject,'Value');
setappdata(handles.fig_bh_FitCov, 'h_GUIdata', h);


function checkbox_mod22t_Callback(hObject, eventdata, handles)
h = getappdata(handles.fig_bh_FitCov, 'h_GUIdata');
h.model_th_code(2,2)=get(hObject,'Value');
setappdata(handles.fig_bh_FitCov, 'h_GUIdata', h);


function checkbox_mod24t_Callback(hObject, eventdata, handles)
h = getappdata(handles.fig_bh_FitCov, 'h_GUIdata');
h.model_th_code(2,4)=get(hObject,'Value');
setappdata(handles.fig_bh_FitCov, 'h_GUIdata', h);


function checkbox_c2t_Callback(hObject, eventdata, handles)
h = getappdata(handles.fig_bh_FitCov, 'h_GUIdata');
h.c_th_code(2,1)=get(hObject,'Value');
setappdata(handles.fig_bh_FitCov, 'h_GUIdata', h);


function checkbox_calc_auto_Callback(hObject, eventdata, handles)


function pushbutton_calcul_covar_Callback(hObject, eventdata, handles)
load_selected_mogs(hObject, eventdata, handles);
updateFigures(handles)


function checkbox_aniso_rot_Callback(hObject, eventdata, handles)
h = getappdata(handles.fig_bh_FitCov, 'h_GUIdata');
h.aniso = get(hObject,'Value');
if h.aniso == 1
    h.aniso = h.aniso+1;
    set(handles.checkbox_aniso,'Value',1)
else
    h.aniso = get(handles.checkbox_aniso,'Value');
end
setappdata(handles.fig_bh_FitCov, 'h_GUIdata', h);
if h.aniso > 0
    if get(handles.popupmenu_type_data,'Value')~=1
        warndlg('Elliptic anisotropy not available for amplitude data');
        set(hObject,'Value',0)
        h.aniso = 0;
        setappdata(handles.fig_bh_FitCov, 'h_GUIdata', h);
        return
    end
end
updateCovAffichage(handles)
CalculRaisDroits(handles)
CalculCovExperimentale(handles)

if ~isempty(h.data)
  updateFigures(handles)
end


function popupmenu11_Callback(hObject, eventdata, handles)


function popupmenu11_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


function popupmenu12_Callback(hObject, eventdata, handles)


function popupmenu12_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit_pepite_th_Callback(hObject, eventdata, handles)
val = str2double(get(hObject,'string'));
if isnan(val)
    h = getappdata(handles.fig_bh_FitCov, 'h_GUIdata');
    errordlg(h.str.s54, h.str.s45,'modal')
	return
end
h = getappdata(handles.fig_bh_FitCov, 'h_GUIdata');
h.nugget_th = val;
h.saved = false;
setappdata(handles.fig_bh_FitCov, 'h_GUIdata', h);
if get(handles.checkbox_calc_auto,'Value')==1, updateFigures(handles), end


function edit_pepite_th_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


function checkbox_pepite_th_Callback(hObject, eventdata, handles)
h = getappdata(handles.fig_bh_FitCov, 'h_GUIdata');
h.nugget_th_code=get(hObject,'Value');
h.saved = false;
setappdata(handles.fig_bh_FitCov, 'h_GUIdata', h);



function edit_pepite_xi_Callback(hObject, eventdata, handles)
val = str2double(get(hObject,'string'));
if isnan(val)
    h = getappdata(handles.fig_bh_FitCov, 'h_GUIdata');
    errordlg(h.str.s54, h.str.s45,'modal')
	return
end
h = getappdata(handles.fig_bh_FitCov, 'h_GUIdata');
h.nugget_xi = val;
h.saved = false;
setappdata(handles.fig_bh_FitCov, 'h_GUIdata', h);
if get(handles.checkbox_calc_auto,'Value')==1, updateFigures(handles), end


function edit_pepite_xi_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


function checkbox_pepite_xi_Callback(hObject, eventdata, handles)
h.nugget_xi_code=get(hObject,'Value');
h.saved = false;
setappdata(handles.fig_bh_FitCov, 'h_GUIdata', h);


function pushbutton_show_vapp_stats_Callback(hObject, eventdata, handles)
h = getappdata(handles.fig_bh_FitCov, 'h_GUIdata');

if isempty(h.L)==1
	CalculRaisDroits(handles);
	h = getappdata(handles.fig_bh_FitCov, 'h_GUIdata');
	if isempty(h.L)==1, return, end
end
nt=length(h.L(:,1));

%Calcul de la lenteur apparente

s = h.data(:,1)./sum(h.L,2);
s0 = mean( s );
vs = var( s );


Tx = h.model3d.grid3d.Tx(h.ind,[1 3]);
Rx = h.model3d.grid3d.Rx(h.ind,[1 3]);
hyp = sqrt( sum((Tx-Rx).^2, 2) );
dz = Tx(:,2)-Rx(:,2);
theta = 180/pi*asin(dz./hyp);

figure
subplot(221)
hist(s,30)
title(['App. slowness, mean = ',num2str(s0),', var = ',num2str(vs)])
subplot(222)
plot(hyp,1./s,'+')
subplot(223)
plot(theta,1./s,'+')
subplot(224)
Tx = unique(Tx,'rows');
Rx = unique(Rx,'rows');
plot(-Tx(:,1),Tx(:,2),'o',-Rx(:,1),Rx(:,2),'+')
%plot(h.iktt,'+')


function fig_bh_FitCov_CloseRequestFcn(hObject, eventdata, handles)
CloseMenuItem_Callback(hObject, eventdata, handles);



function listbox_mogs_Callback(hObject, eventdata, handles)

function listbox_mogs_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function listbox2_Callback(hObject, eventdata, handles)

function listbox2_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function pushbutton4_Callback(hObject, eventdata, handles)

function popupmenu13_Callback(hObject, eventdata, handles)

function popupmenu13_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function checkbox51_Callback(hObject, eventdata, handles)

function checkbox50_Callback(hObject, eventdata, handles)

function edit60_Callback(hObject, eventdata, handles)

function edit60_CreateFcn(hObject, eventdata, handles)

if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function checkbox49_Callback(hObject, eventdata, handles)

function checkbox48_Callback(hObject, eventdata, handles)

function popupmenu14_Callback(hObject, eventdata, handles)

function popupmenu14_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function pushbutton5_Callback(hObject, eventdata, handles)

function checkbox54_Callback(hObject, eventdata, handles)

function checkbox53_Callback(hObject, eventdata, handles)

function edit62_Callback(hObject, eventdata, handles)

function edit62_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function checkbox52_Callback(hObject, eventdata, handles)

function edit61_Callback(hObject, eventdata, handles)

function edit61_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function edit77_Callback(hObject, eventdata, handles)

function edit77_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function edit78_Callback(hObject, eventdata, handles)

function edit78_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit79_Callback(hObject, eventdata, handles)

function edit79_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit80_Callback(hObject, eventdata, handles)

function edit80_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit81_Callback(hObject, eventdata, handles)

function edit81_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit82_Callback(hObject, eventdata, handles)

function edit82_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit83_Callback(hObject, eventdata, handles)

function edit83_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit84_Callback(hObject, eventdata, handles)

function edit84_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit85_Callback(hObject, eventdata, handles)

function edit85_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit86_Callback(hObject, eventdata, handles)

function edit86_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit87_Callback(hObject, eventdata, handles)

function edit87_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit88_Callback(hObject, eventdata, handles)

function edit88_CreateFcn(hObject, eventdata, handles)

if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function checkbox69_Callback(hObject, eventdata, handles)

function checkbox70_Callback(hObject, eventdata, handles)

function checkbox71_Callback(hObject, eventdata, handles)

function checkbox72_Callback(hObject, eventdata, handles)

function checkbox73_Callback(hObject, eventdata, handles)

function checkbox74_Callback(hObject, eventdata, handles)

function checkbox75_Callback(hObject, eventdata, handles)

function checkbox76_Callback(hObject, eventdata, handles)

function checkbox77_Callback(hObject, eventdata, handles)

function checkbox78_Callback(hObject, eventdata, handles)

function checkbox79_Callback(hObject, eventdata, handles)

function checkbox80_Callback(hObject, eventdata, handles)

function popupmenu18_Callback(hObject, eventdata, handles)

function popupmenu18_CreateFcn(hObject, eventdata, handles)

if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function popupmenu19_Callback(hObject, eventdata, handles)

function popupmenu19_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function popupmenu20_Callback(hObject, eventdata, handles)

function popupmenu20_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit65_Callback(hObject, eventdata, handles)

function edit65_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function edit66_Callback(hObject, eventdata, handles)

function edit66_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function edit67_Callback(hObject, eventdata, handles)

function edit67_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit68_Callback(hObject, eventdata, handles)

function edit68_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function checkbox57_Callback(hObject, eventdata, handles)

function checkbox58_Callback(hObject, eventdata, handles)

function checkbox59_Callback(hObject, eventdata, handles)

function checkbox60_Callback(hObject, eventdata, handles)

function popupmenu16_Callback(hObject, eventdata, handles)

function popupmenu16_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function edit92_Callback(hObject, eventdata, handles)

function edit92_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function checkbox81_Callback(hObject, eventdata, handles)

function model16_Callback(hObject, eventdata, handles)

function model16_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function model17_Callback(hObject, eventdata, handles)

function model17_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function checkbox_mod16_Callback(hObject, eventdata, handles)

function checkbox_mod17_Callback(hObject, eventdata, handles)

function checkbox_mod27_Callback(hObject, eventdata, handles)

function checkbox_mod26_Callback(hObject, eventdata, handles)

function model27_Callback(hObject, eventdata, handles)

function model27_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function model26_Callback(hObject, eventdata, handles)

function model26_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function checkbox86_Callback(hObject, eventdata, handles)

function edit97_Callback(hObject, eventdata, handles)

function edit97_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function popupmenu21_Callback(hObject, eventdata, handles)

function popupmenu21_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function checkbox87_Callback(hObject, eventdata, handles)

function checkbox_mod25_Callback(hObject, eventdata, handles)

function checkbox89_Callback(hObject, eventdata, handles)

function checkbox90_Callback(hObject, eventdata, handles)

function edit98_Callback(hObject, eventdata, handles)

function edit98_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function model25_Callback(hObject, eventdata, handles)

function model25_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function edit100_Callback(hObject, eventdata, handles)

function edit100_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function edit101_Callback(hObject, eventdata, handles)

function edit101_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function model13_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
