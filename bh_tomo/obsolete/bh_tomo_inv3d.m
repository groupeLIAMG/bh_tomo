function varargout = bh_tomo_inv3d(varargin)
% BH_TOMO_INV3D M-file for bh_tomo_inv3d.fig
%      BH_TOMO_INV3D, by itself, creates a new BH_TOMO_INV3D or raises the existing
%      singleton*.
%
%      H = BH_TOMO_INV3D returns the handle to a new BH_TOMO_INV3D or the handle to
%      the existing singleton*.
%
%      BH_TOMO_INV3D('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in BH_TOMO_INV3D.M with the given input arguments.
%
%      BH_TOMO_INV3D('Property','Value',...) creates a new BH_TOMO_INV3D or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before bh_tomo_inv_OpeningFunction gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to bh_tomo_inv3d_OpeningFcn via varargin.
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

% Edit the above text to modify the response to help bh_tomo_inv3d

% Last Modified by GUIDE v2.5 22-Feb-2016 19:53:21

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
	'gui_Singleton',  gui_Singleton, ...
	'gui_OpeningFcn', @bh_tomo_inv3d_OpeningFcn, ...
	'gui_OutputFcn',  @bh_tomo_inv3d_OutputFcn, ...
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


% --- Executes just before bh_tomo_inv3d is made visible.
function bh_tomo_inv3d_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to bh_tomo_inv3d (see VARARGIN)

% Choose default command line output for bh_tomo_inv3d
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes bh_tomo_inv3d wait for user response (see UIRESUME)
% uiwait(handles.figure1);
h.db_file = get(handles.fig_bh_inv,'UserData');
h.tt_b_coul_min = 0.06;
h.tt_b_coul_max = 0.12;
h.amp_b_coul_min = 0;
h.amp_b_coul_max = 3;
h.model_initial = [];
str = get_str_locale();

setappdata(handles.fig_bh_inv, 'h', h)
setappdata(handles.fig_bh_inv, 'str', str)

set_string_locale(handles)

% --- Outputs from this function are returned to the command line.
function varargout = bh_tomo_inv3d_OutputFcn(hObject, eventdata, handles)
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


function popupmenu_type_inv_Callback(hObject, eventdata, handles)
h = getappdata(handles.fig_bh_inv, 'h');
val = get(hObject,'Value');
if val==1 % cokri
	set(handles.uipanel_lsqr,'Visible','off')
	set(handles.uipanel_cokri,'Visible','on')
	if ~isfield(h, 'model3d'), return, end
	if get(handles.popupmenu_type_data,'Value')==1 % tt
		if isfield(h.model3d,'tt_covar')
			update_covar_info(handles, h.model3d.tt_covar)
		else
			str = getappdata(handles.fig_bh_inv, 'str');
			uiwait(warndlg(str.s184))
			update_covar_info(handles, [])
			return
		end
	else
		if isfield(h.model3d,'amp_covar')
			update_covar_info(handles, h.model3d.amp_covar)
		else
			str = getappdata(handles.fig_bh_inv, 'str');
			uiwait(warndlg(str.s184))
			update_covar_info(handles, [])
			return
		end
	end
elseif val==2 % lsqr
	set(handles.uipanel_cokri,'Visible','off')
	set(handles.uipanel_lsqr,'Visible','on')
end


function popupmenu_type_inv_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
	set(hObject,'BackgroundColor','white');
end


function edit_cokri_n_iter_Callback(hObject, eventdata, handles)


function edit_cokri_n_iter_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
	set(hObject,'BackgroundColor','white');
end


function edit_lsqr_n_iter_Callback(hObject, eventdata, handles)


function edit_lsqr_n_iter_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
	set(hObject,'BackgroundColor','white');
end


function edit_lsqr_alpha_Callback(hObject, eventdata, handles)


function edit_lsqr_alpha_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
	set(hObject,'BackgroundColor','white');
end


function edit_lsqr_correlmin_Callback(hObject, eventdata, handles)


function edit_lsqr_correlmin_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
	set(hObject,'BackgroundColor','white');
end


function checkbox_cokri_simu_Callback(hObject, eventdata, handles)


function edit_cokri_n_simu_Callback(hObject, eventdata, handles)
val = str2double(get(hObject,'String'));
if isnan(val)
	str = getappdata(handles.fig_bh_inv, 'str');
	errordlg(str.s54, str.s45,'modal')
	return
end
set(hObject,'String', num2str(2^nextpow2(val)))

function edit_cokri_n_simu_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
	set(hObject,'BackgroundColor','white');
end


function checkbox_use_c0_Callback(hObject, eventdata, handles)


function FileMenu_Callback(hObject, eventdata, handles)


function OpenMenuItem_Callback(hObject, eventdata, handles)
h = getappdata(handles.fig_bh_inv, 'h');
str = getappdata(handles.fig_bh_inv, 'str');
[h.no_model3d,h.db_file,hh] = choisirPanneau('UserData', h.db_file,'3d');
if ishandle(hh), delete(hh), end
if h.no_model3d==0
	return
end
load(h.db_file,'model3d')
h.model3d = model3d(h.no_model3d);
if ~isfield(h.model3d, 'grid3d')
	errordlg(str.s150, str.s45,'model3d')
	return
elseif isempty(h.model3d.grid3d)
	errordlg(str.s150, str.s45,'model3d')
	return
end
set(handles.text_nom_panneau,'String',[str.s151,':   ', h.model3d.name])

load(h.db_file,'mogs')
mog_date = {};
for i = 1:length(h.model3d.mogs)
    mog_date{i} = [mogs(h.model3d.mogs(i)).name, ',', mogs(h.model3d.mogs(i)).date];
end
set(handles.listbox_mogs,'String',mog_date)

set(handles.text_xmin,'String',num2str(h.model3d.grid3d.grx(1)))
set(handles.text_xmax,'String',num2str(h.model3d.grid3d.grx(length(h.model3d.grid3d.grx))))
set(handles.text_dx,'String',num2str(h.model3d.grid3d.grx(2)-h.model3d.grid3d.grx(1)))
set(handles.text_ymin,'String',num2str(h.model3d.grid3d.gry(1)))
set(handles.text_ymax,'String',num2str(h.model3d.grid3d.gry(length(h.model3d.grid3d.gry))))
set(handles.text_dy,'String',num2str(h.model3d.grid3d.gry(2)-h.model3d.grid3d.gry(1)))
set(handles.text_zmin,'String',num2str(h.model3d.grid3d.grz(1)))
set(handles.text_zmax,'String',num2str(h.model3d.grid3d.grz(length(h.model3d.grid3d.grz))))
set(handles.text_dz,'String',num2str(h.model3d.grid3d.grz(2)-h.model3d.grid3d.grz(1)))
if get(handles.popupmenu_type_inv,'Value')==1
	if get(handles.popupmenu_type_data,'Value')==1 && isfield(h.model3d,'tt_covar')
		update_covar_info(handles, h.model3d.tt_covar)
	elseif get(handles.popupmenu_type_data,'Value')==2 && isfield(h.model3d,'amp_covar')
		update_covar_info(handles, h.model3d.amp_covar)
	end
end
setappdata(handles.fig_bh_inv, 'h', h)
update_Ldc(handles)


function SaveMenuItem_Callback(hObject, eventdata, handles)
h = getappdata(handles.fig_bh_inv, 'h');
if isfield(h.model3d,'inv_res')
    n = length(h.model3d.inv_res)+1;
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
    case 4
        dType = '-hyb';
end
nameDefault = ['tomo',num2str(n) dType];
prompt = {'Inversion name:                                                               '};
title = 'Save inversion results';
nblines = 1;
answer = inputdlg(prompt,title,nblines,{nameDefault});
name = [];
if ~isempty(answer)
    name=answer{1};
end
%name = get(handles.edit_nom_inv,'String');
if isempty(name)
	str = getappdata(handles.fig_bh_inv, 'str');
	uiwait(errordlg(str.s185,str.s45))
	return
end

if isfield(h.model3d,'inv_res')
	flag=0;
	for n=1:length(h.model3d.inv_res)
		if strcmp(h.model3d.inv_res(n).name, name)
			no_inv_res = n;
			flag = 1;
			break;
		end
	end
	if flag==1
		str = getappdata(handles.fig_bh_inv, 'str');
		rep=questdlg([str.s187,' ',name,'?']);
		if ~strcmp(rep,'Yes')
			return
		end
	else
		no_inv_res = 1+length(h.model3d.inv_res);
	end
else
	no_inv_res = 1;
end
h.model3d.inv_res(no_inv_res).name = name;
h.model3d.inv_res(no_inv_res).tomo = h.tomo;
h.model3d.inv_res(no_inv_res).param = h.param;
if get(handles.popupmenu_type_inv,'Value')==1
	h.model3d.inv_res(no_inv_res).covar = h.covar;
end
load(h.db_file,'model3d')
model3d(h.no_model3d).inv_res = h.model3d.inv_res; %#ok<NASGU>
save(h.db_file,'model3d','-append')
setappdata(handles.fig_bh_inv, 'h', h)

update_Ldc(handles)


function CloseMenuItem_Callback(hObject, eventdata, handles)
delete(handles.fig_bh_inv)


function popupmenu_type_data_Callback(hObject, eventdata, handles)
h = getappdata(handles.fig_bh_inv, 'h');
if get(handles.popupmenu_type_inv,'Value')==1
	if get(hObject,'Value')==1 % tt
        set(handles.checkbox_use_Ldc,'Value',0);
		if isfield(h.model3d,'tt_covar')
			update_covar_info(handles, h.model3d.tt_covar)
		else
			str = getappdata(handles.fig_bh_inv,'str');
			uiwait(warndlg(str.s184))
			update_covar_info(handles, [])
			return
        end
        set(handles.edit_n_iter_courbe,'String','2')
		set(handles.edit_b_coul_min,'String',num2str(h.tt_b_coul_min))
		set(handles.edit_b_coul_max,'String',num2str(h.tt_b_coul_max))
    else
        set(handles.checkbox_use_Ldc,'Value',1);
		if isfield(h.model3d,'amp_covar')
			update_covar_info(handles, h.model3d.amp_covar)
		else
			str = getappdata(handles.fig_bh_inv,'str');
			uiwait(warndlg(str.s184))
			update_covar_info(handles, [])
			return
		end
		set(handles.edit_n_iter_courbe,'String','0')
		set(handles.edit_b_coul_min,'String',num2str(h.amp_b_coul_min))
		set(handles.edit_b_coul_max,'String',num2str(h.amp_b_coul_max))
	end
end


function popupmenu_type_data_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
	set(hObject,'BackgroundColor','white');
end


function edit_n_iter_courbe_Callback(hObject, eventdata, handles)


function edit_n_iter_courbe_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
	set(hObject,'BackgroundColor','white');
end


function pushbutton_inv_Callback(hObject, eventdata, handles)
if get(handles.popupmenu_type_data,'Value')~=1 && get(handles.checkbox_use_Ldc,'Value')~=1
    msgbox('Please select a previous inversion and click checkbox ''Use rays''.','Message','modal');
    return;
end
if get(handles.popupmenu_type_inv,'Value')==1 && ...
		isempty( get(handles.edit_model14,'String') )
	str = getappdata(handles.fig_bh_inv, 'str');
	uiwait(warndlg(str.s184))
	return
end

cla(handles.axes1)
cla(handles.axes2)
cbh = findobj( 0, 'tag', 'Colorbar' );   
for i = 1: length(cbh) 
    colorbar(cbh(i),'off') 
end

h = getappdata(handles.fig_bh_inv, 'h');  

selected_mogs = get(handles.listbox_mogs,'Value'); 
str_mogs = get(handles.listbox_mogs,'String'); 
h.param.model_initial = h.model_initial;
clim = [];
if get(handles.popupmenu_type_data,'Value')==1
	h.param.tomo_amp = 0;    
	[data,ind] = getPanneauData(h.model3d, h.db_file, 'tt',selected_mogs);
	data = [h.model3d.grid3d.Tx(ind,:) h.model3d.grid3d.Rx(ind,:) data ];%...
 %       h.model3d.grid3d.TxCosDir(ind,:) h.model3d.grid3d.RxCosDir(ind,:)];
	if get(handles.checkbox_b_coul,'Value')==1
		clim = [h.tt_b_coul_min h.tt_b_coul_max];
	end
else
	h.param.tomo_amp = 1;
	if get(handles.popupmenu_type_data,'Value')==2
		[data,ind] = getPanneauData(h.model3d, h.db_file,'amp',selected_mogs);
	elseif get(handles.popupmenu_type_data,'Value')==3
		[data,ind] = getPanneauData(h.model3d, h.db_file,'fce',selected_mogs);
	elseif get(handles.popupmenu_type_data,'Value')==4
		[data,ind] = getPanneauData(h.model3d, h.db_file,'hyb',selected_mogs);
	end
	data = [h.model3d.grid3d.Tx(ind,:) h.model3d.grid3d.Rx(ind,:) data];
	if get(handles.checkbox_b_coul,'Value')==1
		clim = [h.amp_b_coul_min h.amp_b_coul_max];
	end
end

h.param.dv_max = 0.01*str2double( get(handles.edit_dv_max, 'String') );
h.param.nbreitrd = str2double(get(handles.edit_n_iter_droit,'String'));
h.param.nbreitrc = str2double(get(handles.edit_n_iter_courbe,'String'));
h.param.nbreitra = str2double(get(handles.edit_n_iter_antenne,'String'));
h.param.it1_rays_horiz = 0;
h.param.plot_xi = 0;
h.param.radar = 1;
%
% correction longueur des antennes
%
% if h.param.nbreitra > 0
% 	% on prend le type d'antenne sur le 1er mog
% 	load(h.db_file,'mogs')  
%     
% 	if ~isempty( regexpi(mogs(h.model3d.mogs(1)).data.antennas, 'ramac') )
% 		type = ['Ramac ', num2str(mogs(h.model3d.mogs(1)).data.rnomfreq, '%d')];
% 		h.param.cla.corr = corrAntL(type);
% 	end
% 	h.param.cla.diam = getPanneauData(h.model3d, h.db_file,'ant');
% 	clear mogs
% end

h.param.pixnonvisite = 'non';
h.param.saveInvData = 1;
Ldc = [];
no_Ldc = get(handles.popupmenu_Ldc,'Value');
if get(handles.checkbox_use_Ldc,'Value')==1
%     if isempty(h.model3d.inv_res)
%     if get(handles.popupmenu_type_data,'Value') ~=1
%         if 
	
	%    ind = zeros(1,length(h.model3d.inv_res(no_Ldc).tomo.no_trace));
	ind = [];
	for n=1:size(data,1)%length(h.model3d.inv_res(no_Ldc).tomo.no_trace)
		ii = find( h.model3d.inv_res(no_Ldc).tomo.no_trace==data(n,9) );
		if isempty(ii)
			str = getappdata(handles.fig_bh_inv, 'str');
			uiwait(warndlg([str.s186{1}, num2str(data(n,9)), str.s186{2}]))
			%Ldc = [];
			return
		else
			ind = [ind ii];
		end
		%            ind(n) = findnear(h.model3d.inv_res(no_Ldc).tomo.no_trace(n), data(:,9));
	end
	Ldc = h.model3d.inv_res(no_Ldc).tomo.L(ind,:);
end

h.param.use_cont = true;
if get(handles.checkbox_use_cont,'Value')==0, h.param.use_cont = false; end

maps = get(handles.popupmenu_colormap,'String');
cmap = maps{get(handles.popupmenu_colormap,'Value')};

gh = {clim; cmap; handles.axes1; ''; ''};
set(handles.axes1,'Visible','on','Position',[0.45 0.081 0.5 0.71])
set(handles.axes2,'Visible','off')

if get(handles.popupmenu_type_inv,'Value')==1
	% cokrigeage
	h.param.ival = 2; % 2 = cokri
	h.param.mode = 1; % 1 = cokrigreage, 2 = simulations
	h.param.nbresimu = 0;
	if get(handles.checkbox_cokri_simu,'Value')==1
		h.param.mode = 2;
		h.param.nbresimu = str2double(get(handles.edit_cokri_n_simu,'String'));
	end
	%	h.covar.model3d = ones(1,4);
	%    h.covar.c = zeros(1,1);
	%	h.covar.nugget_t = 1;
	%	h.covar.nugget_l = 0;
	h.covar.use_c0 = get(handles.checkbox_use_c0, 'Value');
	h.covar.model(1,1) = get(handles.popupmenu_model11,'Value')+1;
	h.covar.model(1,2) = str2double(get(handles.edit_model12,'String'));
	h.covar.model(1,3) = str2double(get(handles.edit_model13,'String'));
	h.covar.model(1,4) = str2double(get(handles.edit_model14,'String'));
    h.covar.model(1,5) = str2double(get(handles.edit_model15,'String'));
	h.covar.model(1,6) = str2double(get(handles.edit_model16,'String'));
	h.covar.model(1,7) = str2double(get(handles.edit_model17,'String'));    
	h.covar.c(1,1) = str2double(get(handles.edit_c1,'String'));
	h.covar.nugget_t = str2double(get(handles.edit_pepite_t_val,'String'));
	h.covar.nugget_l = str2double(get(handles.edit_pepite_l_val,'String'));
	if ~isempty(get(handles.edit_model22,'String'))
		h.covar.model(2,1) = get(handles.popupmenu_model21,'Value')+1;
		h.covar.model(2,2) = str2double(get(handles.edit_model22,'String'));
		h.covar.model(2,3) = str2double(get(handles.edit_model23,'String'));
		h.covar.model(2,4) = str2double(get(handles.edit_model24,'String'));
        h.covar.model(2,5) = str2double(get(handles.edit_model25,'String'));
		h.covar.model(2,6) = str2double(get(handles.edit_model26,'String'));
		h.covar.model(2,7) = str2double(get(handles.edit_model27,'String'));
		h.covar.c(2,1) = str2double(get(handles.edit_c2,'String'));
    end

	if h.param.mode==2 
		gh = {clim; cmap; handles.axes1; handles.axes2; ''};
		set(handles.axes1,'Visible','on','Position',[0.425 0.081 0.19 0.71])
		set(handles.axes2,'Visible','on','Position',[0.75 0.081 0.19 0.71])
    else
        gh = {clim; cmap; handles.axes1; ''; ''};
        set(handles.axes1,'Visible','on','Position',[0.425 0.081 0.5 0.71])
		set(handles.axes2,'Visible','off')
    end
    h.tomo = inv_cokri3D(h.param, data, h.covar, h.model3d.grid3d, Ldc, ...
            handles.text_message_inv, gh);

elseif get(handles.popupmenu_type_inv,'Value')==2
	% LSQR
	h.param.nbreiter = str2double(get(handles.edit_lsqr_n_iter,'String'));
	h.param.alpha = str2double(get(handles.edit_lsqr_alpha,'String'));
	h.param.tol = str2double(get(handles.edit_lsqr_tol,'String'));
	h.param.gradmin = str2double(get(handles.edit_lsqr_gradmin,'String'));
	h.param.correlmin = str2double(get(handles.edit_lsqr_correlmin,'String'));
    set(handles.axes1,'Visible','on','Position',[0.425 0.081 0.5 0.71])
    set(handles.axes2,'Visible','off');

	h.tomo = inv_lsqr3D(h.param, data, h.model3d.grid3d, Ldc, ...
		handles.text_message_inv, gh);
   
end

[mog mdate] = strtok(str_mogs{selected_mogs(1)},',');
h.tomo.date = mdate(2:end);
for n=2:length(selected_mogs)
    [mog mdate] = strtok(str_mogs{selected_mogs(n)},',');
    d = mdate(2:end);
    if ( datenum(d) < datenum(h.tomo.date) )
        h.tomo.date = d;
    end
end
setappdata(handles.fig_bh_inv, 'h', h)
edit_b_coul_min_Callback(hObject, eventdata, handles);
edit_b_coul_max_Callback(hObject, eventdata, handles);
if get(handles.popupmenu_type_data,'Value') == 2
    h.tomo.rays = h.model3d.inv_res(no_Ldc).tomo.rays;
end
set(handles.text_message_inv, 'String', '')
%setappdata(handles.fig_bh_inv, 'h', h)

function update_covar_info(handles, covar)
if isempty( covar )
	set(handles.edit_model12,'String','')
	set(handles.edit_model13,'String','')
	set(handles.edit_model14,'String','')
    set(handles.edit_model15,'String','')
	set(handles.edit_model16,'String','')
	set(handles.edit_model17,'String','')
	set(handles.edit_c1,'String','')
	set(handles.edit_pepite_t_val,'String','')
	set(handles.edit_pepite_l_val,'String','')
    
    set(handles.text_lenteur,'Visible','off')

%    set(handles.checkbox_covar_aniso,'Value',0)
%    set(handles.checkbox_aniso_rot,'Value',0)
    set(handles.edit_pepite_xi_val,'Visible','off')
    set(handles.edit_pepite_th_val,'Visible','off')
    set(handles.text_pepite_xi,'Visible','off')
    set(handles.text_pepite_th,'Visible','off')
else
	set(handles.popupmenu_model11,'Value',covar.model(1,1)-1)
	set(handles.edit_model12,'String',num2str(covar.model(1,2)))
	set(handles.edit_model13,'String',num2str(covar.model(1,3)))
	set(handles.edit_model14,'String',num2str(covar.model(1,4)))
    set(handles.edit_model15,'String',num2str(covar.model(1,5)))
	set(handles.edit_model16,'String',num2str(covar.model(1,6)))
	set(handles.edit_model17,'String',num2str(covar.model(1,7)))
	set(handles.edit_c1,'String',num2str(covar.c(1,1)))
	set(handles.edit_pepite_t_val,'String',num2str(covar.nugget_t))
	set(handles.edit_pepite_l_val,'String',num2str(covar.nugget_l))
	set(handles.checkbox_use_c0, 'Value', covar.use_c0)
	if size(covar.model,1)==2
		set(handles.popupmenu_model21,'Value',covar.model(2,1)-1)
		set(handles.edit_model22,'String',num2str(covar.model(2,2)))
		set(handles.edit_model23,'String',num2str(covar.model(2,3)))
		set(handles.edit_model24,'String',num2str(covar.model(2,4)))
        set(handles.edit_model25,'String',num2str(covar.model(2,5)))
		set(handles.edit_model26,'String',num2str(covar.model(2,6)))
		set(handles.edit_model27,'String',num2str(covar.model(2,7)))
		set(handles.edit_c2,'String',num2str(covar.c(2,1)))
	else
		set(handles.edit_model22,'String','')
		set(handles.edit_model23,'String','')
		set(handles.edit_model24,'String','')
        set(handles.edit_model25,'String','')
		set(handles.edit_model26,'String','')
		set(handles.edit_model27,'String','')
		set(handles.edit_c2,'String','')
	end
end

function edit_lsqr_gradmin_CreateFcn(hObject, eventdata, handles)

function edit_lsqr_tol_CreateFcn(hObject, eventdata, handles)

function edit_nom_inv_Callback(hObject, eventdata, handles)


function edit_nom_inv_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
	set(hObject,'BackgroundColor','white');
end


function popupmenu_Ldc_Callback(hObject, eventdata, handles)


function popupmenu_Ldc_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
	set(hObject,'BackgroundColor','white');
end


function edit_n_iter_courbe_cokri_Callback(hObject, eventdata, handles)


function edit_n_iter_courbe_cokri_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
	set(hObject,'BackgroundColor','white');
end


function text_message_inv_CreateFcn(hObject, eventdata, handles)


function update_Ldc(handles)
h = getappdata(handles.fig_bh_inv, 'h');
texte = {};
if isfield(h.model3d,'inv_res')
	for n=1:length(h.model3d.inv_res)
		texte{n} = [char( h.model3d.inv_res(n).name ), ', ',char( h.model3d.inv_res(n).tomo.date)]; 
	end
end
if ~isempty( texte )
	set(handles.popupmenu_Ldc,'String',texte,'Value',1)  
else
	set(handles.popupmenu_Ldc,'String','--','Value',1) 
end


function ResMenu_Callback(hObject, eventdata, handles)


function RaisMenuItem_Callback(hObject, eventdata, handles)
h = getappdata(handles.fig_bh_inv, 'h');
selected_mogs = get(handles.listbox_mogs,'Value');
if get(handles.checkbox_use_Ldc,'Value')==1
	if get(handles.popupmenu_type_data,'Value')==1
		h.param.tomo_amp = 0;
		[data] = getPanneauData(h.model3d, h.db_file, 'tt',selected_mogs);
	else
		if get(handles.popupmenu_type_data,'Value')==2
			[data] = getPanneauData(h.model3d, h.db_file,'amp',selected_mogs);
		elseif get(handles.popupmenu_type_data,'Value')==3
			[data] = getPanneauData(h.model3d, h.db_file,'fce',selected_mogs);
		elseif get(handles.popupmenu_type_data,'Value')==4
			[data] = getPanneauData(h.model3d, h.db_file,'hyb',selected_mogs);
		end
	end
	no_Ldc = get(handles.popupmenu_Ldc,'Value');
	ind = [];
	for n=1:size(data,1)
		ii = find( h.model3d.inv_res(no_Ldc).tomo.no_trace==data(n,3) );
		if isempty(ii)
			str = getappdata(handles.fig_bh_inv, 'str');
			uiwait(warndlg([str.s186{1}, num2str(data(n,9)), str.s186{2}]))
			%Ldc = [];
			return
		else
			ind = [ind ii];
		end
	end
	for n=1:length(ind)
		t.rays{n} = h.model3d.inv_res(no_Ldc).tomo.rays{ind(n)};
    end
    t.invData = h.model3d.inv_res(no_Ldc).tomo.invData;
elseif isfield(h, 'tomo')
	if ~isfield(h.tomo,'rays')
		return
	end
	t = h.tomo;
else
    msgbox('Please make inversion or load previous inversion.','Message','modal'); 
    return;
end
figure
%ha=plotRais2D(t,true);
for i = 1:length(t.rays)
    ha=plot3(t.rays{i}(:,1),t.rays{i}(:,2),t.rays{i}(:,3));
    hold on;
end
grid on;

set(gca,'XLim',[h.model3d.grid3d.grx(1) h.model3d.grid3d.grx(length(h.model3d.grid3d.grx))])
set(gca,'YLim',[h.model3d.grid3d.gry(1) h.model3d.grid3d.gry(length(h.model3d.grid3d.gry))])
set(gca,'ZLim',[h.model3d.grid3d.grz(1) h.model3d.grid3d.grz(length(h.model3d.grid3d.grz))])
str = getappdata(handles.fig_bh_inv, 'str');
title(str.s179,'FontSize',14)
xlabel(['X ',str.s119],'FontSize',12)
ylabel(['Y ',str.s119],'FontSize',12)
zlabel(str.s120,'FontSize',12)


function TomoMenuItem_Callback(hObject, eventdata, handles)
h = getappdata(handles.fig_bh_inv, 'h');
str = getappdata(handles.fig_bh_inv, 'str');
if isfield(h, 'tomo')
    
    ssplot = reshape((h.tomo.s),length(h.tomo.x),length(h.tomo.y),length(h.tomo.z));
    [nx,ny,nz] = size(ssplot);
    ix = h.tomo.x(round(nx/2));
    iy = h.tomo.y(round(ny/2));
    iz = h.tomo.z(round(nz/2));
    figure;
	if h.param.tomo_amp==0
		slice(h.tomo.x,h.tomo.y,h.tomo.z,1./ssplot,ix,iy,iz);
		titre = [str.s121,' [m/ns]'];
		clim = [h.tt_b_coul_min h.tt_b_coul_max];
	else
		slice(h.tomo.x,h.tomo.y,h.tomo.z,ssplot,ix,iy,iz);
		titre = [str.s177,' [Np/m]'];
		clim = [h.amp_b_coul_min h.amp_b_coul_max];
	end
	
	set(gca,'DataAspectRatio',[1 1 1])
	axis tight
	title(titre,'FontSize',14)
	xlabel(['X ',str.s119],'FontSize',12)
	ylabel(['Y ',str.s119],'FontSize',12)
    zlabel(str.s120,'FontSize',12)
	if get(handles.checkbox_b_coul,'Value')==1, caxis(clim), end
	colorbar
else
    msgbox('Please make inversion or load previous inversion.','Message','modal'); 
    return;
end

function SimuMenuItem_Callback(hObject, eventdata, handles)
h = getappdata(handles.fig_bh_inv, 'h');
str = getappdata(handles.fig_bh_inv, 'str');
if isfield(h, 'tomo')
	if h.param.mode ~=2
		return
    end
	if h.param.tomo_amp==0
		data = reshape(1./(h.tomo.lmoy+h.tomo.sgr),length(h.tomo.x),length(h.tomo.y), length(h.tomo.z));
		ec_type = reshape(std((1./(h.tomo.lmoy+h.tomo.simu)),[],2),length(h.tomo.x),length(h.tomo.y), length(h.tomo.z));
		titre = [str.s121,' [m/ns]'];
		clim = [h.tt_b_coul_min h.tt_b_coul_max];
	else
		data = reshape(h.tomo.lmoy+h.tomo.sgr,length(h.tomo.x),length(h.tomo.y), length(h.tomo.z));
		ec_type = reshape(std((h.tomo.lmoy+h.tomo.simu),[],2),length(h.tomo.x),length(h.tomo.y), length(h.tomo.z));
		titre = [str.s177,' [Np/m]'];
		clim = [h.amp_b_coul_min h.amp_b_coul_max];
	end
	figure
	set(gcf,'Position',[250 200 900 600])
    [nx,ny,nz] = size(data);
    ix = h.tomo.x(round(nx/2));
    iy = h.tomo.y(round(ny/2));
    iz = h.tomo.z(round(nz/2));
	subplot(121)
    slice(h.tomo.x,h.tomo.y,h.tomo.z,data,ix,iy,iz);
	set(gca,'DataAspectRatio',[1 1 1])
	axis tight
	title(titre,'FontSize',12)
	xlabel(['X ',str.s119],'FontSize',10)
    ylabel(['Y ',str.s119],'FontSize',10)
	zlabel(str.s120,'FontSize',10)
	if get(handles.checkbox_b_coul,'Value')==1, caxis(clim), end
	colorbar
	subplot(122)
    slice(h.tomo.x,h.tomo.y,h.tomo.z,ec_type,ix,iy,iz);
	set(gca,'DataAspectRatio',[1 1 1])
	axis tight
	title(str.s43,'FontSize',12)
	xlabel(['X ',str.s119],'FontSize',10)
    ylabel(['Y ',str.s119],'FontSize',10)
	zlabel(str.s120,'FontSize',10)
    if get(handles.checkbox_b_coul,'Value')==1, caxis(clim), end
	colorbar
else
    msgbox('Please make inversion or load previous inversion.','Message','modal'); 
    return;
end


function edit_n_iter_courbe_lsqr_Callback(hObject, eventdata, handles)


function edit_n_iter_courbe_lsqr_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
	set(hObject,'BackgroundColor','white');
end

function set_string_locale(handles)
str = getappdata(handles.fig_bh_inv, 'str');

set(handles.uipanel_info,         'Title',  str.s23)
set(handles.popupmenu_type_data,  'String', str.s156)
set(handles.text_nom_panneau,     'String', [str.s151, ': '])
set(handles.uipanel_n_iter,       'Title', str.s157)
set(handles.popupmenu_type_inv,   'String', str.s161)
set(handles.pushbutton_inv,       'String', upper(str.s33))
set(handles.uipanel_grille,       'Title',  str.s92)
set(handles.text_pas,             'String', lower(str.s162))
%set(handles.text_nom_inv,         'String', str.s163)
set(handles.uipanel_lsqr,         'Title',  [str.s164,' - ',str.s159])
set(handles.uipanel_cokri,        'Title',  [str.s164,' - ',str.s158])
set(handles.text_lsqr_tol,        'String', str.s165)
set(handles.text_lsqr_n_iter,     'String', str.s105)
set(handles.checkbox_cokri_simu,  'String', str.s172)
set(handles.text_cokri_n_simu,    'String', str.s168)
set(handles.uipanel_mod1_covar,   'Title',  str.s93)
set(handles.uipanel_mod2_covar,   'Title',  str.s93)
set(handles.popupmenu_model11,    'String', str.s94)
set(handles.popupmenu_model21,    'String', str.s94)

set(handles.text_model12,         'String', str.s97)
set(handles.text_model13,         'String', str.s977)
set(handles.text_model14,     	  'String', str.s98)
set(handles.text_model15,     	  'String', str.s99)
set(handles.text_model16,     	  'String', str.s991)
set(handles.text_model17,     	  'String', str.s992)
set(handles.text_c1,              'String', str.s100)
set(handles.text_model22,         'String', str.s97)
set(handles.text_model23,         'String', str.s977)
set(handles.text_model24,     	  'String', str.s98)
set(handles.text_model25,     	  'String', str.s99)
set(handles.text_model26,     	  'String', str.s991)
set(handles.text_model27,     	  'String', str.s992)
set(handles.text_c2,              'String', str.s100)

set(handles.text_pepite_t,        'String', str.s101)
set(handles.text_pepite_l,        'String', str.s112)
set(handles.checkbox_use_c0,      'String', str.s107)

set(handles.FileMenu,             'Label',  str.s25)
set(handles.OpenMenuItem,         'Label',  str.s169)
set(handles.SaveMenuItem,         'Label',  str.s29)
set(handles.CloseMenuItem,        'Label',  str.s31)
set(handles.ResMenu,              'Label',  str.s170)
set(handles.TomoMenuItem,         'Label',  str.s171)
set(handles.SimuMenuItem,         'Label',  str.s172)
set(handles.RaisMenuItem,         'Label',  str.s64)
set(handles.DataMenu,             'Label',  str.s23)
set(handles.PlotTTMenuItem,       'Label',  str.s198)
set(handles.PlotAmpMenuItem,      'Label',  str.s199)
set(handles.ProcTTMenuItem,       'Label',  str.s209)
set(handles.ProcAmpMenuItem,      'Label',  str.s210)
set(handles.ExportMenuItem,       'Label',  str.s204)
set(handles.EditGridMenuItem,     'Label',  str.s237)
set(handles.ModeleInitialMenuItem,'Label',  str.s28)
set(handles.uipanel_prev_inv,     'Title',  str.s173)
set(handles.checkbox_use_Ldc,     'String', str.s174)
set(handles.pushbutton_visu_prev, 'String', str.s175)
set(handles.pushbutton_load_prev, 'String', str.s176)
set(handles.text_dv_max,          'String', str.s197)
set(handles.checkbox_use_cont,    'String', str.s208)
set(handles.checkbox_b_coul,      'String', str.s218)
set(handles.text_rais_droits,     'String', str.s166)
set(handles.text_rais_courbes,    'String', str.s153)


load rotate
set(handles.togglebutton_rotate1,'CData',cdata)
set(handles.togglebutton_rotate2,'CData',cdata)
load zoom
set(handles.togglebutton_zoom1,'CData',zoomCData)
set(handles.togglebutton_zoom2,'CData',zoomCData)
clear cdata zoomCData



function pushbutton_visu_prev_Callback(hObject, eventdata, handles)
if isempty( get(handles.popupmenu_Ldc,'String') )
	return
elseif strcmp( get(handles.popupmenu_Ldc,'String'), '--')==1
	return
end
no = get(handles.popupmenu_Ldc,'Value');
h = getappdata(handles.fig_bh_inv, 'h');
str = getappdata(handles.fig_bh_inv, 'str');
tomo = h.model3d.inv_res(no).tomo;
param = h.model3d.inv_res(no).param;
ssplot = reshape((tomo.s),length(tomo.x),length(tomo.y),length(tomo.z));
[nx,ny,nz] = size(ssplot);
ix = round(nx/2);
iy = round(ny/2);
iz = round(nz/2);
figure;
if sum(sum(sum(ssplot))) == 0
	slice(ssplot,iy,ix,iz);
else
    slice(1./ssplot,iy,ix,iz);
end
if param.tomo_amp==0
	titre = [str.s121,' [m/ns]'];
else
	titre = [str.s177,' [Np/m]'];
end

set(gca,'DataAspectRatio',[1 1 1],'YDir','normal')
axis tight
title(titre,'FontSize',14)
xlabel(str.s119,'FontSize',12)
ylabel(str.s119,'FontSize',12)
zlabel(str.s120,'FontSize',8)
colorbar;
suptitle([h.model3d.name,' - ',h.model3d.inv_res(no).name],'Interpreter','none')

function pushbutton_load_prev_Callback(hObject, eventdata, handles)
if isempty( get(handles.popupmenu_Ldc,'String') )
	return
elseif strcmp( get(handles.popupmenu_Ldc,'String'), '--')==1
	return
end
no = get(handles.popupmenu_Ldc,'Value');
h = getappdata(handles.fig_bh_inv, 'h');
h.tomo = h.model3d.inv_res(no).tomo;
h.param = h.model3d.inv_res(no).param;
update_covar_info(handles, h.model3d.inv_res(no).covar);
setappdata(handles.fig_bh_inv, 'h', h)


function edit_dv_max_Callback(hObject, eventdata, handles)


function edit_dv_max_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
	set(hObject,'BackgroundColor','white');
end


function DataMenu_Callback(hObject, eventdata, handles)


function PlotTTMenuItem_Callback(hObject, eventdata, handles)
h = getappdata(handles.fig_bh_inv, 'h');
str = getappdata(handles.fig_bh_inv, 'str');
selected_mogs = get(handles.listbox_mogs,'Value');
[data,ind] = getPanneauData(h.model3d, h.db_file, 'tt',selected_mogs);
if isempty(data)
	uiwait(warndlg(str.s131))
	return
end
data = [h.model3d.grid3d.Tx(ind,:) h.model3d.grid3d.Rx(ind,:) data];
load(h.db_file,'mogs')
mog = mogs(h.model3d.mogs(1));

hyp = sqrt( sum((data(:,1:3)-data(:,4:6)).^2, 2) );
dz = data(:,6)-data(:,3);
theta = 180/pi*asin(dz./hyp);

figure;

subplot(221)
plot(hyp, data(:,7), 'o')
xlabel(str.s35)
ylabel([str.s22,' [',mog.data.tunits,']'])

subplot(222)
plot(theta, hyp./data(:,7), 'o')
xlabel(str.s36)
ylabel([str.s37, '[',mog.data.cunits,'/',mog.data.tunits,']'])
title(str.s39)

subplot(223)
plot(hyp,data(:,8),'o')
xlabel(str.s35)
ylabel(str.s43)

subplot(224)
plot(theta,data(:,8),'o')
xlabel(str.s36)
ylabel(str.s43)



function PlotAmpMenuItem_Callback(hObject, eventdata, handles)
h = getappdata(handles.fig_bh_inv, 'h');
str = getappdata(handles.fig_bh_inv, 'str');
selected_mogs = get(handles.listbox_mogs,'Value');
data = [];
if get(handles.popupmenu_type_data,'Value')==2
	[data,ind] = getPanneauData(h.model3d, h.db_file,'amp',selected_mogs);
elseif get(handles.popupmenu_type_data,'Value')==3
	[data,ind] = getPanneauData(h.model3d, h.db_file,'fce',selected_mogs);
elseif get(handles.popupmenu_type_data,'Value')==4
	[data,ind] = getPanneauData(h.model3d, h.db_file,'hyb',selected_mogs);
end
if isempty(data)
	uiwait(warndlg(str.s131))
	return
end
data = [h.model3d.grid3d.Tx(ind,:) h.model3d.grid3d.Rx(ind,:) data];

hyp = sqrt( sum((data(:,1:3)-data(:,4:6)).^2, 2) );
dz = data(:,6)-data(:,3);
theta = 180/pi*asin(dz./hyp);

figure;

subplot(221)
plot(hyp, data(:,7), 'o')
xlabel(str.s35)
ylabel(str.s200)

subplot(222)
plot(theta, data(:,7)./hyp, 'o')
xlabel(str.s36)
ylabel(str.s201)

subplot(223)
plot(hyp,data(:,8),'o')
xlabel(str.s35)
ylabel(str.s43)

subplot(224)
plot(theta,data(:,8),'o')
xlabel(str.s36)
ylabel(str.s43)


function ExportMenuItem_Callback(hObject, eventdata, handles)
h = getappdata(handles.fig_bh_inv, 'h');
if isfield(h, 'tomo')
	[file, rep] = uiputfile('*.mat','Exporter rsultats');
	if isequal(file,0)
		return
	end
	tomo = h.tomo; %#ok<NASGU>
	param = h.param; %#ok<NASGU>
	save([rep,file],'tomo','param')
end


function edit_b_coul_min_Callback(hObject, eventdata, handles)
val = str2double(get(handles.edit_b_coul_min,'String'));
h = getappdata(handles.fig_bh_inv, 'h');
if get(handles.popupmenu_type_data,'Value')==1
	h.tt_b_coul_min = val;
	if get(handles.checkbox_b_coul,'Value')==1
		caxis(handles.axes1,[h.tt_b_coul_min h.tt_b_coul_max])
        caxis(handles.axes2,[h.tt_b_coul_min h.tt_b_coul_max])
	end
else
	h.amp_b_coul_min = val;
	if get(handles.checkbox_b_coul,'Value')==1
		caxis(handles.axes1,[h.amp_b_coul_min h.amp_b_coul_max])
		caxis(handles.axes2,[h.amp_b_coul_min h.amp_b_coul_max])
	end
end
setappdata(handles.fig_bh_inv, 'h',h);


function edit_b_coul_min_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
	set(hObject,'BackgroundColor','white');
end


function edit_b_coul_max_Callback(hObject, eventdata, handles)
val = str2double(get(handles.edit_b_coul_max,'String'));
h = getappdata(handles.fig_bh_inv, 'h');
if get(handles.popupmenu_type_data,'Value')==1
	h.tt_b_coul_max = val;
	if get(handles.checkbox_b_coul,'Value')==1
		caxis(handles.axes1,[h.tt_b_coul_min h.tt_b_coul_max])
        caxis(handles.axes2,[h.tt_b_coul_min h.tt_b_coul_max])
	end
else
	h.amp_b_coul_max = val;
	if get(handles.checkbox_b_coul,'Value')==1
		caxis(handles.axes1,[h.amp_b_coul_min h.amp_b_coul_max])
		caxis(handles.axes2,[h.amp_b_coul_min h.amp_b_coul_max])
	end
end
setappdata(handles.fig_bh_inv, 'h',h);

function edit_b_coul_max_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
	set(hObject,'BackgroundColor','white');
end


function checkbox_use_Ldc_Callback(hObject, eventdata, handles)
if isempty( get(handles.popupmenu_Ldc,'String') )
	set(hObject,'Value',0)
elseif strcmp( get(handles.popupmenu_Ldc,'String'), '--')==1
	set(hObject,'Value',0)
end


function checkbox_use_cont_Callback(hObject, eventdata, handles)


function ProcTTMenuItem_Callback(hObject, eventdata, handles)
h = getappdata(handles.fig_bh_inv, 'h');
uiwait(bh_tomo_tt(h.db_file))
str = getappdata(handles.fig_bh_inv, 'str');
uiwait(warndlg(str.s211))

function ProcAmpMenuItem_Callback(hObject, eventdata, handles)
h = getappdata(handles.fig_bh_inv, 'h');
uiwait(bh_tomo_amp(h.db_file))
str = getappdata(handles.fig_bh_inv, 'str');
uiwait(warndlg(str.s211))


function EditGridMenuItem_Callback(hObject, eventdata, handles)
h = getappdata(handles.fig_bh_inv, 'h');
str = getappdata(handles.fig_bh_inv, 'str');
if ~isfield(h, 'model3d')
	errordlg(str.s150, str.s45,'model')
	return
elseif ~isfield(h.model3d, 'grid3d')
	errordlg(str.s150, str.s45,'model')
	return
elseif isempty(h.model3d.grid3d)
	errordlg(str.s150, str.s45,'model')
	return
end

if strcmp( h.model3d.type, 'normal' )
	tmp = prepare_grille_data(handles);

	[hh,grid3d] = bh_tomo_grid3d(tmp,h.model3d.grid3d);
	if ishandle(hh), delete(hh), end
	if ~isempty(grid3d.grx)
		h.model3d.grid3d = grid3d;
    end
end
setappdata(handles.fig_bh_inv, 'h',h);

function tmp = prepare_grille_data(handles)
h = getappdata(handles.fig_bh_inv, 'h');
load(h.db_file,'boreholes')
load(h.db_file,'mogs')

for n=1:length(h.model3d.boreholes)
	toto(n) = boreholes(h.model3d.boreholes(n));
end
tmp.boreholes = toto;
tmp.Tx = [];
tmp.Rx = [];
tmp.in = [];
for n=1:length(h.model3d.mogs)
    tmp.in = [tmp.in; mogs(h.model3d.mogs(n)).in'];
	tmp.Tx = [tmp.Tx; [mogs(h.model3d.mogs(n)).data.Tx_x' ...
		mogs(h.model3d.mogs(n)).data.Tx_y' ...
		mogs(h.model3d.mogs(n)).data.Tx_z'] ];
	tmp.Rx = [tmp.Rx; [mogs(h.model3d.mogs(n)).data.Rx_x' ...
		mogs(h.model3d.mogs(n)).data.Rx_y' ...
		mogs(h.model3d.mogs(n)).data.Rx_z'] ];
end
tmp.in = logical(tmp.in);

function InversionMenuItem_Callback(hObject, eventdata, handles)


function ModeleInitialMenuItem_Callback(hObject, eventdata, handles)
h = getappdata(handles.fig_bh_inv, 'h');
[hh,h.model_initial] = bh_tomo_initial3d(h.model3d);
if ishandle(hh), delete(hh), end
setappdata(handles.fig_bh_inv, 'h',h);


function edit_n_iter_droit_Callback(hObject, eventdata, handles)


function edit_n_iter_droit_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
	set(hObject,'BackgroundColor','white');
end


function popupmenu_no_struct_Callback(hObject, eventdata, handles)
no = get(hObject,'Value');
if no==1
	set(handles.uipanel_mod2_covar,'Visible','off')
	set(handles.uipanel_mod1_covar,'Visible','on')
else
	if ~isempty(get(handles.edit_model22,'String'))
		set(handles.uipanel_mod1_covar,'Visible','off')
		set(handles.uipanel_mod2_covar,'Visible','on')
	else
		set(hObject,'Value',1)
	end
end

function popupmenu_no_struct_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
	set(hObject,'BackgroundColor','white');
end


function popupmenu_model21_Callback(hObject, eventdata, handles)


function popupmenu_model21_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
	set(hObject,'BackgroundColor','white');
end


function edit_model22_Callback(hObject, eventdata, handles)


function edit_model22_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
	set(hObject,'BackgroundColor','white');
end


function edit_model23_Callback(hObject, eventdata, handles)


function edit_model23_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
	set(hObject,'BackgroundColor','white');
end


function edit_model24_Callback(hObject, eventdata, handles)


function edit_model24_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
	set(hObject,'BackgroundColor','white');
end


function edit_c2_Callback(hObject, eventdata, handles)


function edit_c2_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
	set(hObject,'BackgroundColor','white');
end


function edit_pepite_l_val_Callback(hObject, eventdata, handles)


function edit_pepite_l_val_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
	set(hObject,'BackgroundColor','white');
end


function fig_bh_inv_CreateFcn(hObject, eventdata, handles)


function fig_bh_inv_DeleteFcn(hObject, eventdata, handles)


function popupmenu_colormap_Callback(hObject, eventdata, handles)
maps = get(hObject,'String');
axes(handles.axes1)
eval(['colormap(',maps{get(hObject,'Value')},')'])
if get(handles.popupmenu_type_data,'Value')==2
	axes(handles.axes2)
	eval(['colormap(',maps{get(hObject,'Value')},')'])
end

function popupmenu_colormap_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
	set(hObject,'BackgroundColor','white');
end


function edit_n_iter_antenne_Callback(hObject, eventdata, handles)


function edit_n_iter_antenne_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


function ResidusMenuItem_Callback(hObject, eventdata, handles)
h = getappdata(handles.fig_bh_inv, 'h');
str = getappdata(handles.fig_bh_inv, 'str');
selected_mogs = get(handles.listbox_mogs,'Value');
if get(handles.popupmenu_type_data,'Value')==1
    [data,ind] = getPanneauData(h.model3d, h.db_file, 'tt',selected_mogs);
else
    [data,ind] = getPanneauData(h.model3d, h.db_file, 'amp',selected_mogs);
end
if isempty(data)
	uiwait(warndlg(str.s131))
	return
end
data = [h.model3d.grid3d.Tx(ind,:) h.model3d.grid3d.Rx(ind,:) data];

hyp = sqrt( sum((data(:,1:3)-data(:,4:6)).^2, 2) );
dz = data(:,6)-data(:,3);
theta = 180/pi*asin(dz./hyp);

if isfield(h, 'tomo')
	if ~isfield(h.tomo,'rays')
		return
	end
	tomo = h.tomo;
else
    msgbox('Please make inversion or load previous inversion.','Message','modal'); 
    return;
end
nIt = length(tomo.invData);
rms = zeros(nIt,1);
for n=1:nIt
    rms(n) = rmsv(tomo.invData(n).res);
end

ns = 2;
if isfield(tomo,'corr_tt_Lant')
if ~isempty(tomo.corr_tt_Lant)
    ns = 3;
end
end
figure
subplot(ns,2,1)
plot(1:nIt, rms,'o')
ylabel('||res||')
xlabel(str.s280)

res = tomo.invData(nIt).res;
subplot(ns,2,2)
plot(theta, res,'o')
xlabel(str.s36)
ylabel(str.s279)

vres = var(res);
h1=subplot(ns,2,3);
hist(res,30)
ylabel(str.s279)
title(['\sigma^2 = ', num2str(vres)])
if get(handles.popupmenu_type_data,'Value')==1
    [depth,ind] = getPanneauData(h.model3d, h.db_file, 'depth',selected_mogs, [], 'tt');
else
    [depth,ind] = getPanneauData(h.model3d, h.db_file, 'depth',selected_mogs, [], 'amp');
end
dTx = sort(unique(depth(:,1)));
dRx = sort(unique(depth(:,2)));
imdata = nan(length(dTx),length(dRx));
for i=1:length(dTx)
    for j=1:length(dRx)
        ind = dTx(i)==depth(:,1) & dRx(j)==depth(:,2);
        if sum(ind)==1
            imdata(i,j) = res(ind);
        end
    end
end

p = [1 0 0;1 1 1;0 0 1];
p = interp1((-1:1)',p,(-1:0.02:1)');

z=imdata; 
z(isnan(imdata))=0; 
z(~isnan(imdata))=1;

set(get(h1,'Children'),'FaceColor',[0 0 1])
if isfield(tomo,'corr_tt_Lant')
    if ~isempty(tomo.corr_tt_Lant)
        subplot(ns,2,5)
        plot(theta, tomo.corr_tt_Lant, 'o')
        xlabel(str.s36)
        ylabel(str.s281)
    
        corr_rel = tomo.corr_tt_Lant./data(:,7);
        subplot(ns,2,6)
        plot(theta, corr_rel, 'o')
        xlabel(str.s36)
        ylabel([str.s281,' [%]'])
    end
end


function edit_model13x_Callback(hObject, eventdata, handles)


function edit_model13x_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit_model12x_Callback(hObject, eventdata, handles)


function edit_model12x_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit_model14x_Callback(hObject, eventdata, handles)


function edit_model14x_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit_c1x_Callback(hObject, eventdata, handles)


function edit_c1x_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit_model23x_Callback(hObject, eventdata, handles)


function edit_model23x_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit_model22x_Callback(hObject, eventdata, handles)


function edit_model22x_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit_model24x_Callback(hObject, eventdata, handles)


function edit_model24x_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit_c2x_Callback(hObject, eventdata, handles)


function edit_c2x_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


function popupmenu_model11x_Callback(hObject, eventdata, handles)


function popupmenu_model11x_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


function popupmenu_model21x_Callback(hObject, eventdata, handles)


function popupmenu_model21x_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function checkbox_show_simu_Callback(hObject, eventdata, handles)


function popupmenu_model11t_Callback(hObject, eventdata, handles)


function popupmenu_model11t_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit_model13t_Callback(hObject, eventdata, handles)


function edit_model13t_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit_model12t_Callback(hObject, eventdata, handles)


function edit_model12t_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit_model14t_Callback(hObject, eventdata, handles)


function edit_model14t_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit_c1t_Callback(hObject, eventdata, handles)


function edit_c1t_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit_model23t_Callback(hObject, eventdata, handles)


function edit_model23t_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit_model22t_Callback(hObject, eventdata, handles)


function edit_model22t_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit_model24t_Callback(hObject, eventdata, handles)


function edit_model24t_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit_c2t_Callback(hObject, eventdata, handles)


function edit_c2t_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


function popupmenu_model21t_Callback(hObject, eventdata, handles)


function popupmenu_model21t_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit_pepite_xi_val_Callback(hObject, eventdata, handles)


function edit_pepite_xi_val_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function edit_pepite_th_val_Callback(hObject, eventdata, handles)


function edit_pepite_th_val_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


function listbox_mogs_Callback(hObject, eventdata, handles)

function listbox_mogs_CreateFcn(hObject, eventdata, handles)

if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


function checkbox_b_coul_Callback(hObject, eventdata, handles)
cbar = get(handles.checkbox_b_coul,'Value'); 
h = getappdata(handles.fig_bh_inv,'h');
if cbar
    if h.param.mode==2 
		axes(handles.axes1); 
        colorbar
        axes(handles.axes2); 
        colorbar
    else
        axes(handles.axes1); 
        colorbar
    end
else
    if h.param.mode==2 
		axes(handles.axes1); 
        colorbar off
        axes(handles.axes2); 
        colorbar off
    else
        axes(handles.axes1); 
        colorbar off
    end
end



function edit_model13_Callback(hObject, eventdata, handles)

function edit_model13_CreateFcn(hObject, eventdata, handles)

if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit_model16_Callback(hObject, eventdata, handles)

function edit_model16_CreateFcn(hObject, eventdata, handles)

if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit_model17_Callback(hObject, eventdata, handles)

function edit_model17_CreateFcn(hObject, eventdata, handles)

if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit_model25_Callback(hObject, eventdata, handles)

function edit_model25_CreateFcn(hObject, eventdata, handles)

if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit_model26_Callback(hObject, eventdata, handles)

function edit_model26_CreateFcn(hObject, eventdata, handles)

if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit_model27_Callback(hObject, eventdata, handles)

function edit_model27_CreateFcn(hObject, eventdata, handles)

if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


function pushbutton_delete_inv_Callback(hObject, eventdata, handles)
handles = guidata(hObject);
h = getappdata(handles.fig_bh_inv, 'h');
ind = get(handles.popupmenu_Ldc,'Value');
h.model3d.inv_res(ind) = [];
setappdata(handles.fig_bh_inv, 'h',h);
load(h.db_file,'model3d')
model3d(h.no_model3d) = h.model3d;
save(h.db_file,'model3d','-append');
update_Ldc(handles);


% --- Executes on button press in togglebutton_zoom1.
function togglebutton_zoom1_Callback(hObject, eventdata, handles)
% hObject    handle to togglebutton_zoom1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of togglebutton_zoom1
button_state = get(hObject,'Value');
if button_state == get(hObject,'Max')
    % toggle button is pressed
	zoom(handles.axes1,'on')
elseif button_state == get(hObject,'Min')
    % toggle button is not pressed
	zoom(handles.axes1,'off')
end


% --- Executes on button press in togglebutton_rotate1.
function togglebutton_rotate1_Callback(hObject, eventdata, handles)
% hObject    handle to togglebutton_rotate1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of togglebutton_rotate1
button_state = get(hObject,'Value');
if button_state == get(hObject,'Max')
    % toggle button is pressed
	rotate3d(handles.axes1,'on')
elseif button_state == get(hObject,'Min')
    % toggle button is not pressed
	rotate3d(handles.axes1,'off')
end


% --- Executes on button press in togglebutton_zoom2.
function togglebutton_zoom2_Callback(hObject, eventdata, handles)
% hObject    handle to togglebutton_zoom2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of togglebutton_zoom2
button_state = get(hObject,'Value');
if button_state == get(hObject,'Max')
    % toggle button is pressed
	zoom(handles.axes2,'on')
elseif button_state == get(hObject,'Min')
    % toggle button is not pressed
	zoom(handles.axes2,'off')
end


% --- Executes on button press in togglebutton_rotate2.
function togglebutton_rotate2_Callback(hObject, eventdata, handles)
% hObject    handle to togglebutton_rotate2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of togglebutton_rotate2
button_state = get(hObject,'Value');
if button_state == get(hObject,'Max')
    % toggle button is pressed
	rotate3d(handles.axes2,'on')
elseif button_state == get(hObject,'Min')
    % toggle button is not pressed
	rotate3d(handles.axes2,'off')
end
