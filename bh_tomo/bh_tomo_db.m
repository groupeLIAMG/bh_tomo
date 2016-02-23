function varargout = bh_tomo_db(varargin)
% BH_TOMO_DB M-file for bh_tomo_db.fig
%      BH_TOMO_DB, by itself, creates a new BH_TOMO_DB or raises the existing
%      singleton*.
%
%      H = BH_TOMO_DB returns the handle to a new BH_TOMO_DB or the handle to
%      the existing singleton*.
%
%      BH_TOMO_DB('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in BH_TOMO_DB.M with the given input arguments.
%
%      BH_TOMO_DB('Property','Value',...) creates a new BH_TOMO_DB or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before bh_tomo_db_OpeningFunction gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to bh_tomo_db_OpeningFcn via varargin.
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

% Edit the above text to modify the response to help bh_tomo_db

% Last Modified by GUIDE v2.5 20-Dec-2012 15:09:04

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @bh_tomo_db_OpeningFcn, ...
                   'gui_OutputFcn',  @bh_tomo_db_OutputFcn, ...
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


% --- Executes just before bh_tomo_db is made visible.
function bh_tomo_db_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to bh_tomo_db (see VARARGIN)

% Choose default command line output for bh_tomo_db
handles.output = [];
handles.second_output = hObject;

h.db_file = get(handles.fig_bh_tomo_db,'UserData');
boreholes = [];
mogs = [];
names_mog = {};
air = [];
panels = [];
par_decime = [];
auto_pick = [];
model3d = [];
h.data_rep = '.';
h.saved = true;

if ~isempty(h.db_file)
	load(h.db_file)
	handles.output = h.db_file;
	
	for no=1:length(boreholes)
		if ~isfield(boreholes(no), 'fdata')
			boreholes(no).fdata = [boreholes(no).X boreholes(no).Y boreholes(no).Z; ...
				boreholes(no).Xmax boreholes(no).Ymax boreholes(no).Zmax]; %#ok<AGROW>
		end
	end
end

% Update handles structure
guidata(hObject, handles);

setappdata(handles.fig_bh_tomo_db,'h',h)

setappdata(handles.fig_bh_tomo_db,'boreholes',boreholes)
setappdata(handles.fig_bh_tomo_db,'mogs',mogs)
setappdata(handles.fig_bh_tomo_db,'names_mog',names_mog)
setappdata(handles.fig_bh_tomo_db,'air',air)
setappdata(handles.fig_bh_tomo_db,'panels',panels)
setappdata(handles.fig_bh_tomo_db,'par_decime',par_decime)
setappdata(handles.fig_bh_tomo_db,'auto_pick',auto_pick)
setappdata(handles.fig_bh_tomo_db,'model3d',model3d)
str = get_str_locale();
setappdata(handles.fig_bh_tomo_db,'str',str)
set_String_locale(handles, str)

update_tout(hObject, eventdata,handles)  %YH

% UIWAIT makes bh_tomo_db wait for user response (see UIRESUME)
% uiwait(handles.fig_bh_tomo_db);



% --- Outputs from this function are returned to the command line.
function varargout = bh_tomo_db_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;
varargout{2} = handles.second_output;


function pushbutton_ajouteForage_Callback(hObject, eventdata, handles)
loc_str = getappdata(handles.fig_bh_tomo_db,'str');
name = inputdlg(loc_str.s145);
if isempty(name)
	return;
end
boreholes = getappdata(handles.fig_bh_tomo_db,'boreholes');
no = length(boreholes)+1;
boreholes(no).name = char( name );
boreholes(no).X = 0;
boreholes(no).Y = 0;
boreholes(no).Z = 0;
boreholes(no).Xmax = 0;
boreholes(no).Ymax = 0;
boreholes(no).Zmax = 0;
boreholes(no).Z_surf = 0;
boreholes(no).scont = [];
boreholes(no).acont = [];
boreholes(no).diam = 0;
boreholes(no).Z_water = NaN;

set(handles.edit_forage_X,'String',num2str( boreholes(no).X ));
set(handles.edit_forage_Y,'String',num2str( boreholes(no).Y ));
set(handles.edit_forage_Z,'String',num2str( boreholes(no).Z ));
set(handles.edit_forage_Xmax,'String',num2str( boreholes(no).Xmax ));
set(handles.edit_forage_Ymax,'String',num2str( boreholes(no).Ymax ));
set(handles.edit_forage_Zmax,'String',num2str( boreholes(no).Zmax ));
set(handles.edit_forage_surface,'String',num2str( boreholes(no).Z_surf ));
set(handles.edit_forage_eau,'String',num2str( boreholes(no).Z_water ));
set(handles.edit_forage_diam,'String',num2str( boreholes(no).diam ));
str = get(handles.listbox_forages,'String');
str{no} = char( boreholes(no).name );
set(handles.listbox_forages,'String',str)
set(handles.listbox_forages,'Value',no)
set(handles.popupmenu_TxCDM,'String',str)
set(handles.popupmenu_RxCDM,'String',str)
h = getappdata(handles.fig_bh_tomo_db,'h');
h.saved = false;
setappdata(handles.fig_bh_tomo_db,'h',h)
setappdata(handles.fig_bh_tomo_db,'boreholes',boreholes);


function listbox_forages_Callback(hObject, eventdata, handles) %#ok<*INUSL>
boreholes = getappdata(handles.fig_bh_tomo_db,'boreholes');
%no = get(hObject,'Value');
no = get(handles.listbox_forages,'Value');  %YH
set(handles.edit_forage_X,'String',num2str( boreholes(no).X ));
set(handles.edit_forage_Y,'String',num2str( boreholes(no).Y ));
set(handles.edit_forage_Z,'String',num2str( boreholes(no).Z ));
set(handles.edit_forage_Xmax,'String',num2str( boreholes(no).Xmax ));
set(handles.edit_forage_Ymax,'String',num2str( boreholes(no).Ymax ));
set(handles.edit_forage_Zmax,'String',num2str( boreholes(no).Zmax ));
set(handles.edit_forage_surface,'String',num2str( boreholes(no).Z_surf ));
set(handles.edit_forage_eau,'String',num2str( boreholes(no).Z_water ));
set(handles.edit_forage_diam,'String',num2str( boreholes(no).diam ));

function listbox_forages_CreateFcn(hObject, eventdata, handles) %#ok<*INUSD>
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


function pushbutton_supprimeForage_Callback(hObject, eventdata, handles)
boreholes = getappdata(handles.fig_bh_tomo_db,'boreholes');
no = get(handles.listbox_forages,'Value');
ind=1:length(boreholes);
ind = ind~=no;

boreholes = boreholes(ind);
str2 = get(handles.listbox_forages,'String');
n1=1;
str = cell(0);
for n=1:length(ind)
	if ind(n)==1
		str{n1} = str2{n};
		n1 = n1+1;
	end
end
setappdata(handles.fig_bh_tomo_db,'boreholes',boreholes);
no = no-1;
if no==0
	no=1;
end
if isempty(boreholes)
	set(handles.listbox_forages,'String','')
	set(handles.popupmenu_TxCDM,'String','')
	set(handles.popupmenu_RxCDM,'String','')
	set(handles.edit_forage_X,'String','');
	set(handles.edit_forage_Y,'String','');
	set(handles.edit_forage_Z,'String','');
	set(handles.edit_forage_Xmax,'String','');
	set(handles.edit_forage_Ymax,'String','');
	set(handles.edit_forage_Zmax,'String','');
	set(handles.edit_forage_surface,'String','');
	set(handles.edit_forage_eau,'String','');
	set(handles.edit_forage_diam,'String','');
else
	if length(boreholes)==1
		set(handles.listbox_forages,'String',str{1})
		set(handles.popupmenu_TxCDM,'String',str{1})
		set(handles.popupmenu_RxCDM,'String',str{1})
	else
		set(handles.listbox_forages,'Value',no,'String',str)
		set(handles.popupmenu_TxCDM,'Value',1,'String',str)
		set(handles.popupmenu_RxCDM,'Value',1,'String',str)
	end
	set(handles.edit_forage_X,'String',num2str( boreholes(no).X ));
	set(handles.edit_forage_Y,'String',num2str( boreholes(no).Y ));
	set(handles.edit_forage_Z,'String',num2str( boreholes(no).Z ));
	set(handles.edit_forage_Xmax,'String',num2str( boreholes(no).Xmax ));
	set(handles.edit_forage_Ymax,'String',num2str( boreholes(no).Ymax ));
	set(handles.edit_forage_Zmax,'String',num2str( boreholes(no).Zmax ));
	set(handles.edit_forage_surface,'String',num2str( boreholes(no).Z_surf ));
	set(handles.edit_forage_eau,'String',num2str( boreholes(no).Z_water ));
	set(handles.edit_forage_diam,'String',num2str( boreholes(no).diam ));
end
h = getappdata(handles.fig_bh_tomo_db,'h');
h.saved = false;
setappdata(handles.fig_bh_tomo_db,'h',h)


function edit_forage_X_Callback(hObject, eventdata, handles)
boreholes = getappdata(handles.fig_bh_tomo_db,'boreholes');
no = get(handles.listbox_forages,'Value');
boreholes(no).X = str2double( get(hObject,'String') );
boreholes(no).fdata = [boreholes(no).X boreholes(no).Y boreholes(no).Z; ...
	boreholes(no).Xmax boreholes(no).Ymax boreholes(no).Zmax];
setappdata(handles.fig_bh_tomo_db,'boreholes',boreholes);
h = getappdata(handles.fig_bh_tomo_db,'h');
h.saved = false;
setappdata(handles.fig_bh_tomo_db,'h',h)
update_CDM_coord(handles,true)


function edit_forage_X_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


function edit_forage_Y_Callback(hObject, eventdata, handles)
boreholes = getappdata(handles.fig_bh_tomo_db,'boreholes');
no = get(handles.listbox_forages,'Value');
boreholes(no).Y = str2double( get(hObject,'String') );
boreholes(no).fdata = [boreholes(no).X boreholes(no).Y boreholes(no).Z; ...
	boreholes(no).Xmax boreholes(no).Ymax boreholes(no).Zmax];
setappdata(handles.fig_bh_tomo_db,'boreholes',boreholes);
h = getappdata(handles.fig_bh_tomo_db,'h');
h.saved = false;
setappdata(handles.fig_bh_tomo_db,'h',h)
update_CDM_coord(handles,true)


function edit_forage_Y_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


function edit_forage_surface_Callback(hObject, eventdata, handles)
boreholes = getappdata(handles.fig_bh_tomo_db,'boreholes');
no = get(handles.listbox_forages,'Value');
boreholes(no).Z_surf = str2double( get(hObject,'String') );
setappdata(handles.fig_bh_tomo_db,'boreholes',boreholes);
h = getappdata(handles.fig_bh_tomo_db,'h');
h.saved = false;
setappdata(handles.fig_bh_tomo_db,'h',h)
update_CDM_coord(handles,true)


function edit_forage_surface_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


function edit_forage_Z_Callback(hObject, eventdata, handles)
boreholes = getappdata(handles.fig_bh_tomo_db,'boreholes');
no = get(handles.listbox_forages,'Value');
boreholes(no).Z = str2double( get(hObject,'String') );
boreholes(no).fdata = [boreholes(no).X boreholes(no).Y boreholes(no).Z; ...
	boreholes(no).Xmax boreholes(no).Ymax boreholes(no).Zmax];
setappdata(handles.fig_bh_tomo_db,'boreholes',boreholes);
h = getappdata(handles.fig_bh_tomo_db,'h');
h.saved = false;
setappdata(handles.fig_bh_tomo_db,'h',h)
update_CDM_coord(handles,true)


function edit_forage_Z_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


function edit_forage_eau_Callback(hObject, eventdata, handles)
boreholes = getappdata(handles.fig_bh_tomo_db,'boreholes');
no = get(handles.listbox_forages,'Value');
boreholes(no).Z_water = str2double( get(hObject,'String') );
setappdata(handles.fig_bh_tomo_db,'boreholes',boreholes);
h = getappdata(handles.fig_bh_tomo_db,'h');
h.saved = false;
setappdata(handles.fig_bh_tomo_db,'h',h)
update_CDM_coord(handles,true)


function edit_forage_eau_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


function edit_forage_diam_Callback(hObject, eventdata, handles)
boreholes = getappdata(handles.fig_bh_tomo_db,'boreholes');
no = get(handles.listbox_forages,'Value');
boreholes(no).diam = str2double( get(hObject,'String') );
setappdata(handles.fig_bh_tomo_db,'boreholes',boreholes);
h = getappdata(handles.fig_bh_tomo_db,'h');
h.saved = false;
setappdata(handles.fig_bh_tomo_db,'h',h)
update_CDM_coord(handles,true)


function edit_forage_diam_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


function pushbutton_ajouteCDM_Callback(hObject, eventdata, handles)
str = getappdata(handles.fig_bh_tomo_db,'str');
boreholes = getappdata(handles.fig_bh_tomo_db,'boreholes');

if isempty( boreholes )
	warndlg('Define boreholes first.')
	return
end

h = getappdata(handles.fig_bh_tomo_db,'h');
old_rep = pwd;
cd( h.data_rep );
[file, rep, filterindex] = uigetfile({
    '*.rad;*.RAD',str.s46;...
	'*.mat',str.s294;...
	'*.hd;*.HD',str.s265;...
	'*.mat',str.s76;...
	'*.dat',str.s47;...
	'*.su;*.SU',str.s268;...
	'*.sgy;*.segy;*.SGY;*.SEGY',str.s293},str.s50);
cd( old_rep );
if isequal(file,0) || isequal(rep,0)
%    disp(['Error: file: ',file,', rep: ',rep])
    return
end
h.data_rep = rep;
h.saved = false;
setappdata(handles.fig_bh_tomo_db,'h',h)

mogs = getappdata(handles.fig_bh_tomo_db,'mogs');
names_mog = getappdata(handles.fig_bh_tomo_db,'names_mog');
par_decime = getappdata(handles.fig_bh_tomo_db,'par_decime');
no = length(mogs)+1;

switch filterindex
	case 1
		name = file(1:end-4);
		data = lisRAMAC2([rep,name]);
	case 2
		name = file(1:end-4);
		load([rep,name]);
		if exist('data','var') ~= 1
			warndlg('The file does not contain a variable named ''data''.')
			return
		end
		if ~isfield(data, 'comment') %#ok<NODEF>
			data.comment = '';
        end
        if strcmp(data.tunits, 's') == 1
           data.timec = data.timec*1e9;
           data.timestp = data.timestp*1e9;
           data.tunits = 'ns';
        end
	case 3
		name = file(1:end-3);
		data = lisEKKO([rep,name]);
	case 4
		name = file(1:end-4);
		tmp = load([rep,name]);
		data = tmp.d;
		clear tmp
	case 5
		name = file(1:end-4);
		data = lisMSIS([rep,name]);
	case 6
		name = file(1:end-3);
		data = lisSU([rep,file]);
    case 7
        ind = strfind(file, '.');
        name = file(1:(ind(end)-1));
        data = lisSEGY([rep,file]);
end

if ~isstruct(data)
	loc_str = getappdata(handles.fig_bh_tomo_db,'str');
	warndlg([loc_str.s146,': ',rep,file])
	return
end
for n=1:length(names_mog)
	if strcmp(names_mog{n}, name) %#ok<ALIGN>
		loc_str = getappdata(handles.fig_bh_tomo_db,'str');
		ButtonName=questdlg(loc_str.s147);
		switch lower(ButtonName),
			case 'yes',
				name=[name,'_2']; %#ok<AGROW>
				break
			otherwise
				return
		end
    end
end
names_mog{no} = char( name );
mogs(no).name = char( name );
%2012-10-17  display date assoiated with selected mog or real date YH
if isfield(data,'date')
    mogs(no).date = data.date;
    set(handles.edit_date,'String',data.date);
else
    mogs(no).date = get(handles.edit_date,'String');  %date;
end


% added BG 2013-05-16
if ~isfield(data,'tdata')
    data.tdata = [];
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
mogs(no).data = data;    % donnees RAMAC
mogs(no).av = [];        % no du tir aerien pre acquisition
mogs(no).ap = [];        % no du tir aerien post acquisition
mogs(no).Tx = 1;         % no du Tx
mogs(no).Rx = 1;         % no du Rx
mogs(no).tt = -1*ones(1,data.ntrace);        % temps d'arrivee
mogs(no).et = -1*ones(1,data.ntrace);        % ecart-type du temps d'arrivee
mogs(no).tt_done = false(1,data.ntrace);     % temps d'arrivee determine (booleen)
if isempty(data.tdata)
	mogs(no).ttTx = [];
	mogs(no).ttTx_done = [];
else
	mogs(no).ttTx = zeros(1,data.ntrace);        % temps d'arrivee
	mogs(no).ttTx_done = false(1,data.ntrace);   % temps d'arrivee determine (booleen)
end
mogs(no).amp_tmin = -1*ones(1,data.ntrace);  % t min fenetre de determination Amplitude
mogs(no).amp_tmax = -1*ones(1,data.ntrace);  % t max fenetre de determination Amplitude
mogs(no).amp_done = false(1,data.ntrace);    % fenetre de determination Amplitude determinee (booleen)
mogs(no).App = zeros(1,data.ntrace);         % Amp pic a pic
mogs(no).fcentroid = zeros(1,data.ntrace);   % freq centroide
mogs(no).scentroid = zeros(1,data.ntrace);   % variance centroide
mogs(no).tauApp = -1*ones(1,data.ntrace);       % amplitudes corrigees - amplitude ratio
mogs(no).tauApp_et = -1*ones(1,data.ntrace);    % ecart-type des amplitudes corrigees
mogs(no).tauFce = -1*ones(1,data.ntrace);       % amplitudes corrigees - freq centroide
mogs(no).tauFce_et = -1*ones(1,data.ntrace);    % ecart-type des amplitudes corrigees
mogs(no).tauHyb = -1*ones(1,data.ntrace);       % amplitudes corrigees - meth. hybride
mogs(no).tauHyb_et = -1*ones(1,data.ntrace);    % ecart-type des amplitudes corrigees
mogs(no).tau_params = [];
mogs(no).fw = [];                            % donnees filtrees par transf. ondelettes
mogs(no).f_et = 1;
mogs(no).amp_name_Ldc = {};
mogs(no).type = 1;                           % X-hole (1) ou VRP (2)
mogs(no).Tx_z_orig = data.Tx_z;
mogs(no).Rx_z_orig = data.Rx_z;
mogs(no).fac_dt = 1;
mogs(no).user_fac_dt = 0;
mogs(no).in = true(1,data.ntrace);
if exist([rep,name,'_bh_ramac.mat'],'file')
	load([rep,name,'_bh_ramac.mat'])
	if exist('data_pick','var')==1
		mogs(no).tt = data_pick;
		mogs(no).et = data_pick_et;
	else
		ind = data_pick_min ~= -1;
		mogs(no).tt(ind) = 0.5*(data_pick_max(ind)+data_pick_min(ind));
		mogs(no).et(ind) = 0.5*(data_pick_max(ind)-data_pick_min(ind));
	end
	mogs(no).tt_done = data_done;
	if exist('data_amp_tmin','var')==1
		mogs(no).amp_tmin = data_amp_tmin;
		mogs(no).amp_tmax = data_amp_tmax;
		mogs(no).amp_done = logical(data_done_amp);
		if exist('data_fcentroid','var')==1
			mogs(no).fcentroid = data_fcentroid;
			mogs(no).scentroid = data_scentroid;
		end
		if exist('data_App','var')==1
			mogs(no).App = data_App;
		end
	end
end

debut=0;
if no>1
	debut = mogs(no-1).no_traces( length(mogs(no-1).no_traces) );
end
mogs(no).no_traces = debut + (1:data.ntrace);
mogs(no).sorted = false;

par_decime(no).sautTx = 0;
par_decime(no).sautRx = 0;
par_decime(no).arrondi = 0;
par_decime(no).use_SB = 0;
par_decime(no).seuil_SB= 0;
par_decime(no).zmin = min([mogs(no).data.Tx_z mogs(no).data.Rx_z]);
par_decime(no).zmax = max([mogs(no).data.Tx_z mogs(no).data.Rx_z]);

str = get(handles.listbox_CDM,'String');
str{no} = char( mogs(no).name );
set(handles.listbox_CDM,'String',str)
set(handles.listbox_CDM,'Value',no)
set(handles.text_av,'String','')
set(handles.text_ap,'String','')
set(handles.popupmenu_type_CDM,'Value',1)
set(handles.popupmenu_RxCDM,'Value',1)
set(handles.popupmenu_TxCDM,'Value',1)
set(handles.edit_Tx_offset,'String',num2str(mogs(no).data.TxOffset))
set(handles.edit_Rx_offset,'String',num2str(mogs(no).data.RxOffset))
set(handles.edit_fac_dt,'String',num2str(mogs(no).fac_dt))
set(handles.checkbox_fac_dt,'Value',mogs(no).user_fac_dt)
set(handles.edit_f_mul_et,'String',num2str(mogs(no).f_et))
set(handles.edit_fnom,'String',num2str(mogs(no).data.rnomfreq))
%set(handles.edit_date,'String',mogs(no).date)  %YH

setappdata(handles.fig_bh_tomo_db,'mogs',mogs);
setappdata(handles.fig_bh_tomo_db,'names_mog',names_mog);
setappdata(handles.fig_bh_tomo_db,'par_decime',par_decime);

update_CDM_coord(handles)
update_info(handles)


function pushbutton_supprimeCDM_Callback(hObject, eventdata, handles)
mogs = getappdata(handles.fig_bh_tomo_db,'mogs');
names_mog = getappdata(handles.fig_bh_tomo_db,'names_mog');
panels = getappdata(handles.fig_bh_tomo_db,'panels');

names_mog_panels = [];
for n=1:length(panels)
	for nn=1:length(panels(n).mogs)
		names_mog_panels(n).name{nn} = names_mog{ panels(n).mogs(nn) }; %#ok<AGROW>
	end
end

no = get(handles.listbox_CDM,'Value');
ind=1:length(mogs);
ind = ind~=no;

mogs = mogs(ind);
str2 = get(handles.listbox_CDM,'String');
n1=1;
names_mog = cell(0);
for n=1:length(ind)
	if ind(n)==1
		names_mog{n1} = str2{n};
		n1 = n1+1;
	end
end

for n=1:length(panels)
	new_mogs = [];
	for nn=1:length(panels(n).mogs)
		for nnn=1:length(names_mog)
			if strcmp(names_mog_panels(n).name{nn}, names_mog{nnn})
				new_mogs = [new_mogs nnn]; %#ok<AGROW>
				break;
			end
		end
	end
	panels(n).mogs = new_mogs;
end

setappdata(handles.fig_bh_tomo_db,'mogs',mogs)
setappdata(handles.fig_bh_tomo_db,'names_mog',names_mog)
setappdata(handles.fig_bh_tomo_db,'panels',panels)
no = no-1;
if no==0
	no=1;
end
if isempty(mogs)
	set(handles.listbox_CDM,'String','')
	set(handles.text_av,'String','');
	set(handles.text_ap,'String','');
else
	if length(mogs)==1
		set(handles.listbox_CDM,'Value',no,'String',names_mog{1})
	else
		set(handles.listbox_CDM,'Value',no,'String',names_mog)
	end
	air = getappdata(handles.fig_bh_tomo_db,'air');
	if ~isempty(mogs(no).av), set(handles.text_av,'String',air( mogs(no).av ).name), end
	if ~isempty(mogs(no).ap), set(handles.text_ap,'String',air( mogs(no).ap ).name), end
	set(handles.edit_Tx_offset,'String',num2str(mogs(no).data.TxOffset))
	set(handles.edit_Rx_offset,'String',num2str(mogs(no).data.RxOffset))
	set(handles.edit_fac_dt,'String',num2str(mogs(no).fac_dt))
	set(handles.checkbox_fac_dt,'Value',mogs(no).user_fac_dt)
	set(handles.edit_f_mul_et,'String',num2str(mogs(no).f_et))
	set(handles.edit_fnom,'String',num2str(mogs(no).data.rnomfreq))
%    set(handles.edit_date,'String',mogs(no).date)  %YH
end
h = getappdata(handles.fig_bh_tomo_db,'h');
h.saved = false;
setappdata(handles.fig_bh_tomo_db,'h',h)


function listbox_CDM_Callback(hObject, eventdata, handles)
%no = get(hObject,'Value');
no = get(handles.listbox_CDM,'Value');  %YH
mogs = getappdata(handles.fig_bh_tomo_db,'mogs');
air = getappdata(handles.fig_bh_tomo_db,'air');
if isempty(mogs(no).av)
	set(handles.text_av,'String','')
else
	set(handles.text_av,'String',air( mogs(no).av ).name)
end
if isempty(mogs(no).ap)
	set(handles.text_ap,'String','')
else
	set(handles.text_ap,'String',air( mogs(no).ap ).name)
end
if isempty(mogs(no).Tx)
	set(handles.popupmenu_TxCDM,'Value',1)
else
	set(handles.popupmenu_TxCDM,'Value',mogs(no).Tx)
end
if isempty(mogs(no).Rx)
	set(handles.popupmenu_RxCDM,'Value',1)
else
	set(handles.popupmenu_RxCDM,'Value',mogs(no).Rx)
end
set(handles.edit_f_mul_et,'String',num2str(mogs(no).f_et))
set(handles.popupmenu_type_CDM,'Value', mogs(no).type)
set(handles.edit_fac_dt,'String',num2str(mogs(no).fac_dt))
set(handles.checkbox_fac_dt,'Value',mogs(no).user_fac_dt)
set(handles.edit_Tx_offset,'String',num2str(mogs(no).data.TxOffset))
set(handles.edit_Rx_offset,'String',num2str(mogs(no).data.RxOffset))
set(handles.edit_fnom,'String',num2str(mogs(no).data.rnomfreq))
%%2012-10-17  YH
set(handles.edit_date,'String',mogs(no).date)
%%%%%%%%%%%%%%%%%%%%%%%%%%%

function listbox_CDM_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


function popupmenu_TxCDM_Callback(hObject, eventdata, handles)
mogs = getappdata(handles.fig_bh_tomo_db,'mogs');
Tx = get(hObject,'Value');
no = get(handles.listbox_CDM,'Value');
mogs(no).Tx = Tx;
setappdata(handles.fig_bh_tomo_db,'mogs',mogs);
update_CDM_coord(handles)


function popupmenu_TxCDM_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


function popupmenu_RxCDM_Callback(hObject, eventdata, handles)
mogs = getappdata(handles.fig_bh_tomo_db,'mogs');
Rx = get(hObject,'Value');
no = get(handles.listbox_CDM,'Value');
mogs(no).Rx = Rx;
setappdata(handles.fig_bh_tomo_db,'mogs',mogs);
update_CDM_coord(handles)


function popupmenu_RxCDM_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


function pushbutton_av_Callback(hObject, eventdata, handles)
str = getappdata(handles.fig_bh_tomo_db,'str');

h = getappdata(handles.fig_bh_tomo_db,'h');
old_rep = pwd;
cd( h.data_rep );
[file, rep] = uigetfile({'*.rad;*.RAD',str.s46;'*.dat',str.s47},str.s48);
cd( old_rep );
if isequal(file,0) || isequal(rep,0)
    return
end
h.data_rep = rep;
h.saved = false;
setappdata(handles.fig_bh_tomo_db,'h',h)

mogs = getappdata(handles.fig_bh_tomo_db,'mogs');
air = getappdata(handles.fig_bh_tomo_db,'air');
no = get(handles.listbox_CDM,'Value');
name = file(1:end-4);
found = false;
for n=1:length(air)
	if strcmp( char(air(n).name), char( name ) )
		mogs(no).av = n;
		found = true;
		break
	end
end
if ~found
	n=length(air)+1;
	data = lisRAMAC2([rep,name]);
	if ~isstruct(data)
		loc_str = getappdata(handles.fig_bh_tomo_db,'str');
		warndlg([loc_str.s146,': ',rep,file])
		return
	end
	
	loc_str = getappdata(handles.fig_bh_tomo_db,'str');
	l = inputdlg(loc_str.s148);
	try
		d = eval(l{1});
	catch 
		try
			d = eval(['[',l{1},']']);
			if length(d) ~= data.ntrace
				ME = MException('MATLAB:INVFMT', 'Number of positions inconsistent with number of traces');
				throw(ME);
			end
		catch Error
			rethrow(Error);
		end
	end
	
	air(n).data = data;
	air(n).name = char( name );
	air(n).tt = -1*ones(1,data.ntrace);        % temps d'arrivee
	air(n).et = -1*ones(1,data.ntrace);        % ecart-type du temps d'arrivee
	air(n).tt_done = false(1,data.ntrace);     % temps d'arrivee determine (booleen)
%	air(n).d_TxRx = getDistLog([rep,name]);
    air(n).d_TxRx = d;
	air(n).fac_dt = 1;
	air(n).in = true(1,data.ntrace);
	if length(d)==1
		air(n).method = 'fixed_antenna';           % 'fixed_antenna' ou 'walkaway'
	else
		air(n).method = 'walkaway';
	end
%	if isempty( air(n).d_TxRx ) %#ok<ALIGN>
%		loc_str = getappdata(handles.fig_bh_tomo_db,'str');
%		l = inputdlg(loc_str.s148);
%		if isempty(l), return, end
%       if isempty(l{1}), return, end
%       air(n).d_TxRx = str2double( l{1} );
%    end

	mogs(no).av = n;
	setappdata(handles.fig_bh_tomo_db,'air',air);
end
set(handles.text_av,'String',name)
setappdata(handles.fig_bh_tomo_db,'mogs',mogs);


function pushbutton_ap_Callback(hObject, eventdata, handles)
str = getappdata(handles.fig_bh_tomo_db,'str');
h = getappdata(handles.fig_bh_tomo_db,'h');
old_rep = pwd;
cd( h.data_rep );
[file, rep] = uigetfile({'*.rad;*.RAD',str.s46;'*.dat',str.s47},str.s51);
cd( old_rep );
if isequal(file,0) || isequal(rep,0)
    return
end
h.data_rep = rep;
h.saved = false;
setappdata(handles.fig_bh_tomo_db,'h',h)

mogs = getappdata(handles.fig_bh_tomo_db,'mogs');
air = getappdata(handles.fig_bh_tomo_db,'air');
no = get(handles.listbox_CDM,'Value');
name = file(1:end-4);
found = false;
for n=1:length(air)
	if strcmp( char(air(n).name), char( name ) )
		mogs(no).ap = n;
		found = true;
		break
	end
end
if ~found
	n=length(air)+1;
	data = lisRAMAC2([rep,name]);
	if ~isstruct(data) %#ok<ALIGN>
		loc_str = getappdata(handles.fig_bh_tomo_db,'str');
		warndlg([loc_str.s146,': ',rep,file])
		return
    end

    loc_str = getappdata(handles.fig_bh_tomo_db,'str');
	l = inputdlg(loc_str.s148);
	try
		d = eval(l{1});
	catch 
		try
			d = eval(['[',l{1},']']);
			if length(d) ~= data.ntrace
				ME = MException('MATLAB:INVFMT', 'Number of positions inconsistent with number of traces');
				throw(ME);
			end
		catch Error
			rethrow(Error);
		end
	end
	
	air(n).data = data;
	air(n).name = char( name );
	air(n).tt = -1*ones(1,data.ntrace);        % temps d'arrivee
	air(n).et = -1*ones(1,data.ntrace);        % ecart-type du temps d'arrivee
	air(n).tt_done = false(1,data.ntrace);     % temps d'arrivee determine (booleen)
%	air(n).d_TxRx = getDistLog([rep,name]);
    air(n).d_TxRx = d;
	air(n).fac_dt = 1;
	air(n).in = true(1,data.ntrace);
	if length(d)==1
		air(n).method = 'fixed_antenna';           % 'fixed_antenna' ou 'walkaway'
	else
		air(n).method = 'walkaway';
	end
%	if isempty( air(n).d_TxRx )
%		loc_str = getappdata(handles.fig_bh_tomo_db,'str');
%		l = inputdlg(loc_str.s148);
%		if isempty(l), return, end
%        if isempty(l{1}), return, end
%        air(n).d_TxRx = str2double( l{1} );
%	end
	mogs(no).ap = n;
	setappdata(handles.fig_bh_tomo_db,'air',air);
end
set(handles.text_ap,'String',name)
setappdata(handles.fig_bh_tomo_db,'mogs',mogs);


function FileMenu_Callback(hObject, eventdata, handles)


function OpenMenuItem_Callback(hObject, eventdata, handles)
loc_str = getappdata(handles.fig_bh_tomo_db,'str');
[file, rep] = uigetfile('*.mat',loc_str.s27);
if isequal(file,0)
	return
end
h = getappdata(handles.fig_bh_tomo_db,'h');
h.db_file = [rep,file];

set(handles.fig_bh_tomo_db, 'Name', ['bh_tomo_db - ', h.db_file])
load(h.db_file)
if exist('boreholes','var')~=1
	return
end
h.saved = true;
setappdata(handles.fig_bh_tomo_db,'h',h);
setappdata(handles.fig_bh_tomo_db,'mogs',mogs);
setappdata(handles.fig_bh_tomo_db,'names_mog',names_mog);
setappdata(handles.fig_bh_tomo_db,'air',air);
setappdata(handles.fig_bh_tomo_db,'boreholes',boreholes);
setappdata(handles.fig_bh_tomo_db,'panels',panels);
setappdata(handles.fig_bh_tomo_db,'par_decime',par_decime);
setappdata(handles.fig_bh_tomo_db,'model3d',model3d);  %YH
update_tout(hObject, eventdata,handles)  %YH

handles.output = h.db_file;
guidata(hObject, handles);

function SaveMenuItem_Callback(hObject, eventdata, handles)
h = getappdata(handles.fig_bh_tomo_db,'h');
if isempty(h.db_file)
	[file, rep] = uiputfile('*.mat','Enregistrer DB');
	if isequal(file,0)
		return
	end
	h.db_file = [rep,file];
	setappdata(handles.fig_bh_tomo_db,'h',h)
	update_info(handles)
end
mogs = getappdata(handles.fig_bh_tomo_db,'mogs'); %#ok<NASGU>
names_mog = getappdata(handles.fig_bh_tomo_db,'names_mog'); %#ok<NASGU>
air = getappdata(handles.fig_bh_tomo_db,'air'); %#ok<NASGU>
boreholes = getappdata(handles.fig_bh_tomo_db,'boreholes'); %#ok<NASGU>
panels = getappdata(handles.fig_bh_tomo_db,'panels'); %#ok<NASGU>
par_decime = getappdata(handles.fig_bh_tomo_db,'par_decime'); %#ok<NASGU>
auto_pick = getappdata(handles.fig_bh_tomo_db,'auto_pick'); %#ok<NASGU>
model3d = getappdata(handles.fig_bh_tomo_db,'model3d');
save(h.db_file,'names_mog','mogs','air','boreholes','panels','par_decime','auto_pick','model3d'); %YH
h.saved = true;
setappdata(handles.fig_bh_tomo_db,'h',h)


function pushbutton_ajoute_Panneau_Callback(hObject, eventdata, handles)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  YH
panels = getappdata(handles.fig_bh_tomo_db, 'panels');
if ~isempty(panels)
    n = length(panels)+1;
else
    n = 1;
end
loc_str = getappdata(handles.fig_bh_tomo_db,'str');
prompt = loc_str.s149;
nameDefault = ['panel',num2str(n)]; 
names = get(handles.listbox_model3d,'String');
for i =1:length(names)
    if strcmp(names(i),nameDefault)
        nameDefault = ['panel',num2str(n+1)];
        break;
    end
end
title = 'Add panel';
nblines = 1;
answer = inputdlg(prompt,title,nblines,{nameDefault});
name = [];
if ~isempty(answer)
    name=answer{1};
end

% loc_str = getappdata(handles.fig_bh_tomo_db,'str');
% name = inputdlg(loc_str.s149);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if isempty(name)
	return
end
panels = getappdata(handles.fig_bh_tomo_db,'panels');
no = length(panels)+1;
panels(no).name = char( name );
panels(no).mogs = [];
panels(no).boreholes = [];
panels(no).type = 'normal';

str = get(handles.listbox_Panneaux,'String');
if iscell( str )
	no2 = 1+length( str );
	str{no2} = char( panels(no).name );
else
	str = { str; panels(no).name };
	no2 = 2;
end
set(handles.listbox_Panneaux,'String',str)
set(handles.listbox_Panneaux,'Value',no2)
set(handles.listbox_CDM_Panneaux,'String','')
setappdata(handles.fig_bh_tomo_db,'panels',panels);
update_info(handles)
h = getappdata(handles.fig_bh_tomo_db,'h');
h.saved = false;
setappdata(handles.fig_bh_tomo_db,'h',h)


function listbox_Panneaux_Callback(hObject, eventdata, handles)
panels = getappdata(handles.fig_bh_tomo_db,'panels');
names_mog = getappdata(handles.fig_bh_tomo_db,'names_mog');

%no = get_no_panneau(handles, hObject);     %YH
no = get(handles.listbox_Panneaux,'Value');  %YH
str = cell(1,length(panels(no).mogs));
for n=1:length(panels(no).mogs)
	str{n} = char( names_mog{ panels(no).mogs(n) } );
end
set(handles.listbox_CDM_Panneaux,'Value',1)
set(handles.listbox_CDM_Panneaux,'String',str)


function listbox_Panneaux_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


function pushbutton_supprime_Panneau_Callback(hObject, eventdata, handles)
panels = getappdata(handles.fig_bh_tomo_db,'panels');
no = get_no_panneau(handles, handles.listbox_Panneaux);
ind=1:length(panels);
ind = ind~=no;
panels = panels(ind);
liste = cell(0);
nn=1;
for n=1:length(panels)
	if strcmp( panels(n).type, 'normal' )
		liste{nn} = char( panels(n).name );
		nn=nn+1;
	end
end
set(handles.listbox_Panneaux,'String',liste,'Value',length(liste))
setappdata(handles.fig_bh_tomo_db,'panels',panels)
no=no-1;
if isempty(panels)
    set(handles.listbox_CDM_Panneaux,'String','')
    return
elseif no==0
    no=1;
end
names_mog = getappdata(handles.fig_bh_tomo_db,'names_mog');
str = cell(1,length(panels(no).mogs));
for n=1:length(panels(no).mogs)
	str{n} = char( names_mog{ panels(no).mogs(n) } );
end
%set(handles.listbox_CDM_Panneaux,'Value',no)  %%YH
set(handles.listbox_CDM_Panneaux,'Value',1,'String',str) %%YH
h = getappdata(handles.fig_bh_tomo_db,'h');
h.saved = false;
setappdata(handles.fig_bh_tomo_db,'h',h)


function pushbutton_ajoute_CDM_Panneaux_Callback(hObject, eventdata, handles)
names_mog = getappdata(handles.fig_bh_tomo_db,'names_mog');
mogs = getappdata(handles.fig_bh_tomo_db,'mogs');
[no_mog,h] = ajouteMOG(names_mog);
delete(h)
if no_mog==0
    return
end
panels = getappdata(handles.fig_bh_tomo_db,'panels');
no = get_no_panneau(handles, handles.listbox_Panneaux);
nn = 1+length(panels(no).mogs);
if nn==1
    panels(no).boreholes(nn) = mogs(no_mog).Tx;
    panels(no).boreholes(nn+1) = mogs(no_mog).Rx;
else
    ii=find(panels(no).boreholes==mogs(no_mog).Tx,1);
    if isempty(ii) %#ok<ALIGN>
		str = getappdata(handles.fig_bh_tomo_db,'str');
		uiwait(errordlg(str.s181))
		return
	end
    ii=find(panels(no).boreholes==mogs(no_mog).Rx,1);
    if isempty(ii) %#ok<ALIGN>
		str = getappdata(handles.fig_bh_tomo_db,'str');
		uiwait(errordlg(str.s181))
		return
	end
end
panels(no).mogs(nn) = no_mog;
str = get(handles.listbox_CDM_Panneaux,'String');
str{nn} = char( names_mog(no_mog) );
set(handles.listbox_CDM_Panneaux,'String',str)
set(handles.listbox_CDM_Panneaux,'Value',nn)
setappdata(handles.fig_bh_tomo_db,'panels',panels);
h = getappdata(handles.fig_bh_tomo_db,'h');
h.saved = false;
setappdata(handles.fig_bh_tomo_db,'h',h)

%
% update grid if it exists
%
if isfield(panels(no), 'grid')
    if ~isempty(panels(no).grid)
        pushbutton_edit_grille_Callback(hObject, eventdata, handles)
    end
end

function listbox_CDM_Panneaux_Callback(hObject, eventdata, handles)


function listbox_CDM_Panneaux_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


function pushbutton_supprime_CDM_Panneaux_Callback(hObject, eventdata, handles)
panels = getappdata(handles.fig_bh_tomo_db,'panels');
no = get_no_panneau(handles, handles.listbox_Panneaux);
no_CDM = get(handles.listbox_CDM_Panneaux,'Value');
ind=1:length(panels(no).mogs);
ind = ind~=no_CDM;
panels(no).mogs = panels(no).mogs(ind);
setappdata(handles.fig_bh_tomo_db,'panels',panels);
str = get(handles.listbox_CDM_Panneaux,'String');
nn=1;
str2=cell(0);
for n=1:length(ind)
	if ind(n)
		str2{nn} = str{n};
		nn=nn+1;
	end
end
no_CDM = no_CDM-1;
if no_CDM==0
	no_CDM=1;
end
set(handles.listbox_CDM_Panneaux,'Value',no_CDM,'String',str2);
h = getappdata(handles.fig_bh_tomo_db,'h');
h.saved = false;
setappdata(handles.fig_bh_tomo_db,'h',h)

function update_tout(hObject, eventdata,handles) 
names_mog = getappdata(handles.fig_bh_tomo_db,'names_mog');
boreholes = getappdata(handles.fig_bh_tomo_db,'boreholes');
panels = getappdata(handles.fig_bh_tomo_db,'panels');
mogs = getappdata(handles.fig_bh_tomo_db,'mogs');
model3d = getappdata(handles.fig_bh_tomo_db,'model3d');  %YH
str = cell(0);
for n=1:length(boreholes)
	str{n} = char( boreholes(n).name );
    listbox_forages_Callback(hObject, eventdata, handles);  %YH
end
if ~isempty(str)
	set(handles.listbox_forages,'String',str)
	set(handles.popupmenu_TxCDM,'String',str)
	set(handles.popupmenu_RxCDM,'String',str)
end
set(handles.listbox_CDM,'String',names_mog)
%%%%%%%%%%%%%%%%%%%%%  YH
%display date 
if ~isempty(mogs)
    set(handles.edit_date,'String',mogs(1).date);
    listbox_CDM_Callback(hObject, eventdata, handles);   
else
    set(handles.edit_date,'String',date);
end
%display model3d
str = cell(0);
nn=1;
for n=1:length(model3d)
	if strcmp( model3d(n).type, 'normal' )
		str{nn} = char( model3d(n).name );
		nn=nn+1;
	end
end
set(handles.listbox_model3d,'String',str)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

str = cell(0);
nn=1;
for n=1:length(panels)
	if strcmp( panels(n).type, 'normal' )
		str{nn} = char( panels(n).name );
		nn=nn+1;
	end
end
set(handles.listbox_Panneaux,'String',str)
str = cell(0);
nn=1;
for n=1:length(panels)
	if strcmp( panels(n).type, 'super' )
		str{nn} = char( panels(n).name );
		nn=nn+1;
	end
end
set(handles.listbox_SPanneaux,'String',str)
update_info(handles)
%update listbox_CDM_panels and Slistbox_CDM_panels 
if ~isempty(panels)
    listbox_Panneaux_Callback(hObject, eventdata, handles);
    listbox_SPanneaux_Callback(hObject, eventdata, handles);
else
    set(handles.listbox_CDM_Panneaux,'String','')
    set(handles.listbox_CDM_SPanneaux,'String','')
end
if ~isempty(model3d)
    listbox_model3d_Callback(hObject, eventdata, handles);
else
    set(handles.listbox_CDM_model3d,'String','')
end

function QuitMenuItem_Callback(hObject, eventdata, handles)
h = getappdata(handles.fig_bh_tomo_db,'h');
if h.saved == false
	str = getappdata(handles.fig_bh_tomo_db,'str');
	ButtonName=questdlg(str.s236);
	switch ButtonName,
		case 'Yes',
			SaveMenuItem_Callback(hObject, eventdata, handles)
		case 'No',
		case 'Cancel',
			return
	end % switch
end
%guidata(hObject, handles);
%uiresume;
delete(handles.fig_bh_tomo_db)

function update_info(handles)
h = getappdata(handles.fig_bh_tomo_db,'h');
mogs = getappdata(handles.fig_bh_tomo_db,'mogs');
boreholes = getappdata(handles.fig_bh_tomo_db,'boreholes');
panels = getappdata(handles.fig_bh_tomo_db,'panels');
model3d = getappdata(handles.fig_bh_tomo_db,'model3d');   %YH
loc_str = getappdata(handles.fig_bh_tomo_db,'str');
str{1} = '';
if isempty(h.db_file)
	str{2} = 'DB: ';
else
	tmp = findstr('/',h.db_file);
	if ~isempty(tmp)
		bn = h.db_file((tmp(length(tmp))+1):length(h.db_file));
	else
		tmp = findstr('\',h.db_file); %windoze
		if ~isempty(tmp)
			bn = h.db_file((tmp(length(tmp))+1):length(h.db_file));
		else
			bn=h.db_file;
		end
	end
	
	str{2} = ['DB: ',bn];
end
str{3} = '';
str{4} = [num2str(length(boreholes)),' ',loc_str.s132];
nm = length(mogs);
str{5} = [num2str(nm),' ',loc_str.s136];
if nm > 0
	str{6} = [num2str(mogs(nm).no_traces(length( mogs(nm).no_traces))),' Traces'];
else
	str{6} = '0 Traces';
end
str{7} = [num2str(length(panels)),' ',loc_str.s135];
%%%%%%%YH
if length(model3d)==1
 str{8} = [num2str(length(model3d)),' ',loc_str.s296];  
else
    str{8} = [num2str(length(model3d)),' ',loc_str.s296 's']; 
end
%%%%%%%%%%%%%%%%%%
set(handles.text_info,'String',str)



function SaveAsMenuItem_Callback(hObject, eventdata, handles)
h = getappdata(handles.fig_bh_tomo_db,'h');
[file, rep] = uiputfile('*.mat','Enregistrer DB');
if isequal(file,0)
	return
end
h.db_file = [rep,file];
h.saved = true;
setappdata(handles.fig_bh_tomo_db,'h',h)
SaveMenuItem_Callback(hObject, eventdata, handles)
update_info(handles)


function pushbutton_creer_grille_Callback(hObject, eventdata, handles)
no = get_no_panneau(handles, handles.listbox_Panneaux);
% if check_forages_verticaux(handles, no )
%     str = getappdata(handles.fig_bh_tomo_db,'str');
%     uiwait(errordlg(str.s189))
%     return
% end
tmp = prepare_grille_data(handles, no);
[hh,grid] = bh_tomo_grille(tmp);
if ishandle(hh), delete(hh), end
panels = getappdata(handles.fig_bh_tomo_db,'panels');
no = get_no_panneau(handles, handles.listbox_Panneaux);
if ~isempty(grid.grx)
    panels(no).grid = grid;
    setappdata(handles.fig_bh_tomo_db,'panels', panels)
end
h = getappdata(handles.fig_bh_tomo_db,'h');
h.saved = false;
setappdata(handles.fig_bh_tomo_db,'h',h)


function pushbutton_edit_grille_Callback(hObject, eventdata, handles)
panels = getappdata(handles.fig_bh_tomo_db,'panels');
no = get_no_panneau(handles, handles.listbox_Panneaux);
tmp = prepare_grille_data(handles, no);
%if isfield(panels(no), 'grid')
if isfield(panels(no), 'grid') && ~isempty(panels(no).grid)  %2012-10-17  YH
    [hh,grid] = bh_tomo_grille(tmp,panels(no).grid);
else
    [hh,grid] = bh_tomo_grille(tmp);
end
if ishandle(hh), delete(hh), end
if ~isempty(grid.grx)
    panels(no).grid = grid;
    setappdata(handles.fig_bh_tomo_db,'panels', panels)
end
h = getappdata(handles.fig_bh_tomo_db,'h');
h.saved = false;
setappdata(handles.fig_bh_tomo_db,'h',h)


function edit_forage_Xmax_Callback(hObject, eventdata, handles)
boreholes = getappdata(handles.fig_bh_tomo_db,'boreholes');
no = get(handles.listbox_forages,'Value');
boreholes(no).Xmax = str2double( get(hObject,'String') );
boreholes(no).fdata = [boreholes(no).X boreholes(no).Y boreholes(no).Z; ...
	boreholes(no).Xmax boreholes(no).Ymax boreholes(no).Zmax];
setappdata(handles.fig_bh_tomo_db,'boreholes',boreholes);
h = getappdata(handles.fig_bh_tomo_db,'h');
h.saved = false;
setappdata(handles.fig_bh_tomo_db,'h',h)
update_CDM_coord(handles,true)


function edit_forage_Xmax_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


function edit_forage_Ymax_Callback(hObject, eventdata, handles)
boreholes = getappdata(handles.fig_bh_tomo_db,'boreholes');
no = get(handles.listbox_forages,'Value');
boreholes(no).Ymax = str2double( get(hObject,'String') );
boreholes(no).fdata = [boreholes(no).X boreholes(no).Y boreholes(no).Z; ...
	boreholes(no).Xmax boreholes(no).Ymax boreholes(no).Zmax];
setappdata(handles.fig_bh_tomo_db,'boreholes',boreholes);
h = getappdata(handles.fig_bh_tomo_db,'h');
h.saved = false;
setappdata(handles.fig_bh_tomo_db,'h',h)
update_CDM_coord(handles,true)


function edit_forage_Ymax_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


function edit_forage_Zmax_Callback(hObject, eventdata, handles)
boreholes = getappdata(handles.fig_bh_tomo_db,'boreholes');
no = get(handles.listbox_forages,'Value');
boreholes(no).Zmax = str2double( get(hObject,'String') );
boreholes(no).fdata = [boreholes(no).X boreholes(no).Y boreholes(no).Z; ...
	boreholes(no).Xmax boreholes(no).Ymax boreholes(no).Zmax];
setappdata(handles.fig_bh_tomo_db,'boreholes',boreholes);
h = getappdata(handles.fig_bh_tomo_db,'h');
h.saved = false;
setappdata(handles.fig_bh_tomo_db,'h',h)
update_CDM_coord(handles,true)


function edit_forage_Zmax_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


%
% -------------------------------------------------------------------------
%
function set_String_locale(handles, str)

set(handles.pushbutton_ajouteForage,          'String', str.s128)
set(handles.pushbutton_ajouteCDM,             'String', str.s128)
set(handles.pushbutton_ajoute_Panneau,        'String', str.s128)
set(handles.pushbutton_ajoute_CDM_Panneaux,   'String', str.s128)
set(handles.pushbutton_ajoute_SPanneau,       'String', str.s128)
set(handles.pushbutton_ajoute_CDM_SPanneaux,  'String', str.s128)
set(handles.pushbutton_supprimeForage,        'String', str.s129)
set(handles.pushbutton_supprimeCDM,           'String', str.s129)
set(handles.pushbutton_supprime_Panneau,      'String', str.s129)
set(handles.pushbutton_supprime_CDM_Panneaux, 'String', str.s129)
set(handles.pushbutton_supprime_SPanneau,     'String', str.s129)
set(handles.pushbutton_supprime_CDM_SPanneaux,'String', str.s129)
set(handles.uipanel_forages,                  'Title', str.s132)
set(handles.uipanel_info,                     'Title', str.s133)
set(handles.uipanel_CDM,                      'Title', str.s134)
set(handles.uipanel_panneaux,                 'Title', str.s135)
set(handles.uipanel_CDMs_Panneaux,            'Title', str.s136)
set(handles.uipanel_grille_Panneaux,          'Title', str.s92)
set(handles.uipanel_S_panneaux,               'Title', str.s137)
set(handles.uipanel_CDMs_SPanneaux,           'Title', str.s136)
set(handles.uipanel_grille_SPanneaux,         'Title', str.s92)
set(handles.text_coord,                       'String', str.s138)
set(handles.text_forage_collet,               'String', str.s139)
set(handles.text_forage_fond,                 'String', str.s140)
set(handles.pushbutton_av,                    'String', str.s141)
set(handles.pushbutton_ap,                    'String', str.s142)
set(handles.pushbutton_creer_grille,          'String', str.s143)
set(handles.pushbutton_edit_grille,           'String', str.s122)
set(handles.pushbutton_creer_Sgrille,         'String', str.s143)
set(handles.pushbutton_edit_Sgrille,          'String', str.s122)
set(handles.FileMenu,                         'Label', str.s25)
set(handles.OpenMenuItem,                     'Label', str.s118)
set(handles.SaveMenuItem,                     'Label', str.s29)
set(handles.SaveAsMenuItem,                   'Label', str.s144)
set(handles.QuitMenuItem,                     'Label', str.s31)
set(handles.text_f_mul_et,                    'String', str.s152)
set(handles.popupmenu_type_CDM,               'String', {str.s202; str.s203})
set(handles.pushbutton_export_tt,             'String', str.s66)
set(handles.pushbutton_export_tau,            'String', str.s263)
set(handles.pushbutton_vel_con,               'String', str.s250)
set(handles.pushbutton_att_con,               'String', str.s251)
set(handles.pushbutton_couverture,            'String', str.s179)
set(handles.pushbutton_renameCDM,             'String', str.s260)
set(handles.pushbutton_importerCDM,           'String', str.s215)
set(handles.pushbutton_importerBH,            'String', str.s215)
set(handles.pushbutton_plotBH,                'String', str.s262)
set(handles.text_Tx_offset,                   'String', str.s269)
set(handles.text_Rx_offset,                   'String', str.s274)
set(handles.checkbox_fac_dt,                  'String', str.s270)
set(handles.pushbutton_elaguerMOG,            'String', str.s275)
set(handles.edit_date,'String',date)
set(handles.uipanel_CDMs_model3d,             'Title', str.s136)  %YH

function pushbutton_ajoute_CDM_SPanneaux_Callback(hObject, eventdata, handles)
names_mog = getappdata(handles.fig_bh_tomo_db,'names_mog');
mogs = getappdata(handles.fig_bh_tomo_db,'mogs');
[no_mog,h] = ajouteMOG(names_mog);
delete(h)
if no_mog==0
    return
end
panels = getappdata(handles.fig_bh_tomo_db,'panels');
no = get_no_panneau(handles, handles.listbox_SPanneaux);
nn = 1+length(panels(no).mogs);
panels(no).mogs(nn) = no_mog;
if nn==1
    panels(no).boreholes(nn) = mogs(no_mog).Tx;
    panels(no).boreholes(nn+1) = mogs(no_mog).Rx;
else
    ii=find(panels(no).boreholes==mogs(no_mog).Tx,1);
    if isempty(ii), panels(no).boreholes = [panels(no).boreholes mogs(no_mog).Tx]; end
    ii=find(panels(no).boreholes==mogs(no_mog).Rx,1);
    if isempty(ii), panels(no).boreholes = [panels(no).boreholes mogs(no_mog).Rx]; end
end
        
str = get(handles.listbox_CDM_SPanneaux,'String');
str{nn} = char( names_mog(no_mog) );
set(handles.listbox_CDM_SPanneaux,'String',str)
set(handles.listbox_CDM_SPanneaux,'Value',nn)
setappdata(handles.fig_bh_tomo_db,'panels',panels);
h = getappdata(handles.fig_bh_tomo_db,'h');
h.saved = false;
setappdata(handles.fig_bh_tomo_db,'h',h)


function pushbutton_supprime_CDM_SPanneaux_Callback(hObject, eventdata, handles)
panels = getappdata(handles.fig_bh_tomo_db,'panels');
no = get_no_panneau(handles, handles.listbox_SPanneaux);
no_CDM = get(handles.listbox_CDM_SPanneaux,'Value');
ind=1:length(panels(no).mogs);
ind = ind~=no_CDM;
panels(no).mogs = panels(no).mogs(ind);
setappdata(handles.fig_bh_tomo_db,'panels',panels);
str = get(handles.listbox_CDM_SPanneaux,'String');
nn=1;
str2=cell(0);
for n=1:length(ind)
	if ind(n)
		str2{nn} = str{n};
		nn=nn+1;
	end
end
no_CDM = no_CDM-1;
if no_CDM==0
	no_CDM=1;
end
set(handles.listbox_CDM_SPanneaux,'Value',no_CDM,'String',str2);
h = getappdata(handles.fig_bh_tomo_db,'h');
h.saved = false;
setappdata(handles.fig_bh_tomo_db,'h',h)



function listbox_CDM_SPanneaux_Callback(hObject, eventdata, handles)


function listbox_CDM_SPanneaux_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


function pushbutton_creer_Sgrille_Callback(hObject, eventdata, handles)
no = get_no_panneau(handles, handles.listbox_SPanneaux);
if check_boreholes_verticaux(handles, no )
    str = getappdata(handles.fig_bh_tomo_db,'str');
    uiwait(errordlg(str.s189))
    return
end
tmp = prepare_grille_data(handles, no);
[hh,grid] = bh_tomo_Sgrille(tmp);
if ishandle(hh), delete(hh), end
panels = getappdata(handles.fig_bh_tomo_db,'panels');
no = get_no_panneau(handles, handles.listbox_SPanneaux);
if ~isempty(grid.grx)
    panels(no).grid = grid;
    setappdata(handles.fig_bh_tomo_db,'panels', panels)
end
h = getappdata(handles.fig_bh_tomo_db,'h');
h.saved = false;
setappdata(handles.fig_bh_tomo_db,'h',h)


function pushbutton_edit_Sgrille_Callback(hObject, eventdata, handles)
panels = getappdata(handles.fig_bh_tomo_db,'panels');
no = get_no_panneau(handles, handles.listbox_SPanneaux);
tmp = prepare_grille_data(handles, no);
if isfield(panels(no), 'grid')
    [hh,grid] = bh_tomo_Sgrille(tmp,panels(no).grid);
else
    [hh,grid] = bh_tomo_Sgrille(tmp);
end
if ishandle(hh), delete(hh), end
if ~isempty(grid.grx)
    panels(no).grid = grid;
    setappdata(handles.fig_bh_tomo_db,'panels', panels)
end
h = getappdata(handles.fig_bh_tomo_db,'h');
h.saved = false;
setappdata(handles.fig_bh_tomo_db,'h',h)


function pushbutton_ajoute_SPanneau_Callback(hObject, eventdata, handles)
%%%%%%%%%%%%%%%%%%%%%%%% YH
panels = get(handles.listbox_SPanneaux,'String');%getappdata(handles.fig_bh_tomo_db, 'panels');
if ~isempty(panels)
    n = length(panels)+1;
else
    n = 1;
end
loc_str = getappdata(handles.fig_bh_tomo_db,'str');
prompt = loc_str.s149;
nameDefault = ['panel',num2str(n)]; 
names = panels;%get(handles.listbox_model3d,'String');
for i =1:length(names)
    if strcmp(names(i),nameDefault)
        nameDefault = ['panel',num2str(n+1)];
        break;
    end
end
title = 'Add panel';
nblines = 1;
answer = inputdlg(prompt,title,nblines,{nameDefault});
name = [];
if ~isempty(answer)
    name=answer{1};
end

% loc_str = getappdata(handles.fig_bh_tomo_db,'str');
% name = inputdlg(loc_str.s149);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if isempty(name)
	return
end
panels = getappdata(handles.fig_bh_tomo_db,'panels');
no = length(panels)+1;
panels(no).name = char( name );
panels(no).mogs = [];
panels(no).boreholes = [];
panels(no).type = 'super';

str = get(handles.listbox_SPanneaux,'String');
if iscell( str )
	no2 = 1+length( str );
	str{no2} = char( panels(no).name );
else
	str = { str; panels(no).name };
	no2 = 2;
end
set(handles.listbox_SPanneaux,'String',str)
set(handles.listbox_SPanneaux,'Value',no2)
set(handles.listbox_CDM_SPanneaux,'String','')
setappdata(handles.fig_bh_tomo_db,'panels',panels);
update_info(handles)
h = getappdata(handles.fig_bh_tomo_db,'h');
h.saved = false;
setappdata(handles.fig_bh_tomo_db,'h',h)


function pushbutton_supprime_SPanneau_Callback(hObject, eventdata, handles)
panels = getappdata(handles.fig_bh_tomo_db,'panels');
no = get_no_panneau(handles, handles.listbox_SPanneaux);
ind=1:length(panels);
ind = ind~=no;
panels = panels(ind);
liste = cell(0);
nn=1;
for n=1:length(panels)
	if strcmp( panels(n).type, 'super' )
		liste{nn} = char( panels(n).name );
		nn=nn+1;
	end
end
set(handles.listbox_SPanneaux,'String',liste,'Value',length(liste))
setappdata(handles.fig_bh_tomo_db,'panels',panels)
no=no-1;
if isempty(panels)
    set(handles.listbox_CDM_SPanneaux,'String','')
    return
elseif no==0
    no=1;
end
names_mog = getappdata(handles.fig_bh_tomo_db,'names_mog');
str = cell(1,length(panels(no).mogs));
for n=1:length(panels(no).mogs)
	str{n} = char( names_mog{ panels(no).mogs(n) } );
end
set(handles.listbox_CDM_SPanneaux,'Value',no)
set(handles.listbox_CDM_SPanneaux,'String',str)
h = getappdata(handles.fig_bh_tomo_db,'h');
h.saved = false;
setappdata(handles.fig_bh_tomo_db,'h',h)


function listbox_SPanneaux_Callback(hObject, eventdata, handles)
panels = getappdata(handles.fig_bh_tomo_db,'panels');
names_mog = getappdata(handles.fig_bh_tomo_db,'names_mog');
%no = get_no_panneau(handles, hObject);
no = get(handles.listbox_SPanneaux,'Value');
str = cell(1,length(panels(no).mogs));
if strcmp(panels(no).type,'super')   %YH
for n=1:length(panels(no).mogs)
	str{n} = char( names_mog{ panels(no).mogs(n) } );
end
set(handles.listbox_CDM_SPanneaux,'Value',1)
set(handles.listbox_CDM_SPanneaux,'String',str)
end

function listbox_SPanneaux_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function tmp = prepare_grille_data(handles, no)
boreholes = getappdata(handles.fig_bh_tomo_db,'boreholes');
mogs = getappdata(handles.fig_bh_tomo_db,'mogs');
panels = getappdata(handles.fig_bh_tomo_db,'panels');

%no = get_no_panneau(handles, handles.listbox_panels);
tmp.boreholes = [];
for n=1:length(panels(no).boreholes)
    tmp.boreholes = [tmp.boreholes boreholes(panels(no).boreholes(n))];
end
tmp.Tx = [];
tmp.Rx = [];
tmp.TxCosDir = [];
tmp.RxCosDir = [];
tmp.Tx_Z_water = [];
tmp.Rx_Z_water = [];
tmp.in = [];
for n=1:length(panels(no).mogs)
    tmp.in = [tmp.in; mogs(panels(no).mogs(n)).in'];
    tmp.Tx = [tmp.Tx; [mogs(panels(no).mogs(n)).data.Tx_x' ...
        mogs(panels(no).mogs(n)).data.Tx_y' ...
        mogs(panels(no).mogs(n)).data.Tx_z'] ];
    tmp.Rx = [tmp.Rx; [mogs(panels(no).mogs(n)).data.Rx_x' ...
        mogs(panels(no).mogs(n)).data.Rx_y' ...
        mogs(panels(no).mogs(n)).data.Rx_z'] ];
    
    tmp.TxCosDir = [tmp.TxCosDir; mogs(panels(no).mogs(n)).TxCosDir];
    tmp.RxCosDir = [tmp.RxCosDir; mogs(panels(no).mogs(n)).RxCosDir];
    
    trou = boreholes(mogs(panels(no).mogs(n)).Tx);
    
    if ~isnan(trou.Z_water) 
        x = interp1(trou.fdata(:,3), trou.fdata(:,1), trou.Z_water);
        y = interp1(trou.fdata(:,3), trou.fdata(:,2), trou.Z_water);
    else
        x = NaN;
        y = NaN;
    end
    tmp.Tx_Z_water = [tmp.Tx_Z_water; [x y trou.Z_water] ];
        
    trou = boreholes(mogs(panels(no).mogs(n)).Rx);
    if ~isnan(trou.Z_water)
        x = interp1(trou.fdata(:,3), trou.fdata(:,1), trou.Z_water);
        y = interp1(trou.fdata(:,3), trou.fdata(:,2), trou.Z_water);
    else
        x = NaN;
        y = NaN;
    end
    tmp.Rx_Z_water = [tmp.Rx_Z_water; [x y trou.Z_water] ];
    
end
tmp.in = logical(tmp.in);

% function tmp = prepare_Sgrille_data(handles)
% boreholes = getappdata(handles.fig_bh_tomo_db,'boreholes');
% mogs = getappdata(handles.fig_bh_tomo_db,'mogs');
% panels = getappdata(handles.fig_bh_tomo_db,'panels');
% 
% no = get_no_panneau(handles, handles.listbox_Spanels);
% tmp.boreholes = [];
% for n=1:length(panels(no).boreholes)
%     tmp.boreholes = [tmp.boreholes boreholes(panels(no).boreholes(n))];
% end
% tmp.Tx = [];
% tmp.Rx = [];
% for n=1:length(panels(no).mogs)
%     tmp.Tx = [tmp.Tx; [mogs(panels(no).mogs(n)).data.Tx_x' ...
%     mogs(panels(no).mogs(n)).data.Tx_y' ...
%     mogs(panels(no).mogs(n)).data.Tx_z'] ];
%     tmp.Rx = [tmp.Rx; [mogs(panels(no).mogs(n)).data.Rx_x' ...
%     mogs(panels(no).mogs(n)).data.Rx_y' ...
%     mogs(panels(no).mogs(n)).data.Rx_z'] ];
% end


function pushbutton_tr_ZOP_Callback(hObject, eventdata, handles)
no = get(handles.listbox_CDM,'Value');
mogs = getappdata(handles.fig_bh_tomo_db,'mogs');
if isempty( strfind( mogs(no).data.csurvmod, 'VRP' ) )
	air = getappdata(handles.fig_bh_tomo_db,'air');
    boreholes = getappdata(handles.fig_bh_tomo_db,'boreholes');
    bh_tomo_ZOP(mogs(no), air, boreholes);
else
	str = getappdata(handles.fig_bh_tomo_db,'str');
	uiwait(warndlg(str.s180))
end

function edit_f_mul_et_Callback(hObject, eventdata, handles)
val = str2double(get(hObject,'String'));
if isnan(val)
    str = getappdata(handles.fig_bh_tomo_db, 'str');
    errordlg(str.s54, str.s45,'modal')
    return
end
mogs = getappdata(handles.fig_bh_tomo_db,'mogs');
no = get(handles.listbox_CDM,'Value');
mogs(no).f_et = val;
setappdata(handles.fig_bh_tomo_db,'mogs', mogs)
h = getappdata(handles.fig_bh_tomo_db,'h');
h.saved = false;
setappdata(handles.fig_bh_tomo_db,'h',h)


function edit_f_mul_et_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


function no = get_no_panneau(handles,hObject)
panels = getappdata(handles.fig_bh_tomo_db,'panels');
str = get(hObject,'String');
no=0;
if iscell( str )
	name = str{ get(hObject,'Value') };
else
	name = str;
end
for n=1:length(panels)
	if strcmp(name, panels(n).name)
		no=n;
		return
	end
end

function flag = check_boreholes_verticaux(handles, no_panneau)
flag = false;
panels = getappdata(handles.fig_bh_tomo_db,'panels');
boreholes = getappdata(handles.fig_bh_tomo_db,'boreholes');
for n=1:length(panels(no_panneau).boreholes )
    nf = panels(no_panneau).boreholes(n);
    if abs(boreholes(nf).X-boreholes(nf).Xmax)>1.0e-5 || ...
            abs(boreholes(nf).Y-boreholes(nf).Ymax)>1.0e-5
        flag = true;
        return
    end
end


function pushbutton_stats_tt_Callback(hObject, eventdata, handles)
no = get(handles.listbox_CDM,'Value');
mogs = getappdata(handles.fig_bh_tomo_db,'mogs');
air = getappdata(handles.fig_bh_tomo_db,'air');
str = getappdata(handles.fig_bh_tomo_db,'str');

tt = mogs(no).tt;
et = mogs(no).et;
done = mogs(no).tt_done & mogs(no).in;

ind = tt~=-1 & mogs(no).in;

if sum(done)==0
	uiwait(errordlg('Data not processed!!!'))
	return
end

if strcmp( mogs(no).data.tunits,'ns' )==1
    av = air( mogs(no).av );
    ap = air( mogs(no).ap );
    
    if mogs(no).data.synthetique==1
        t0 = zeros(1,length(mogs(no).tt));
        fac_dt_av = 1;
        fac_dt_ap = 1;
    else
        [t0, fac_dt_av, fac_dt_ap] = corr_t0(length(mogs(no).tt), av, ap, true);
    end
    
    if ~isempty(av), air( mogs(no).av ).fac_dt = fac_dt_av; end
    if ~isempty(ap), air( mogs(no).ap ).fac_dt = fac_dt_ap; end
    if mogs(no).user_fac_dt == 0
        if fac_dt_av~=1 && fac_dt_ap ~= 1
            mogs(no).fac_dt = 0.5*(fac_dt_av+fac_dt_ap);
        elseif fac_dt_av~=1
            mogs(no).fac_dt = fac_dt_av;
        elseif fac_dt_ap~=1
            mogs(no).fac_dt = fac_dt_ap;
        else
            mogs(no).fac_dt = 1;
        end
        set(handles.edit_fac_dt,'String',num2str(mogs(no).fac_dt));
    end
end
setappdata(handles.fig_bh_tomo_db,'air', air);
setappdata(handles.fig_bh_tomo_db,'mogs', mogs);

hyp = sqrt( (mogs(no).data.Tx_x(ind)-mogs(no).data.Rx_x(ind)).^2 + ...
	(mogs(no).data.Tx_y(ind)-mogs(no).data.Rx_y(ind)).^2 + ...
	(mogs(no).data.Tx_z(ind)-mogs(no).data.Rx_z(ind)).^2 );
dz = mogs(no).data.Rx_z(ind)-mogs(no).data.Tx_z(ind);
theta = 180/pi*asin(dz./hyp);

tt = tt(ind);
et = et(ind);

tt = tt*mogs(no).fac_dt;
et = et*mogs(no).fac_dt;
t0 = t0*mogs(no).fac_dt;

h_stat = figure;
set(h_stat, 'Position',[256 71 512 620]);

subplot(321)
plot(hyp, tt, 'o')
xlabel(str.s35)
ylabel([str.s22,' [',mogs(no).data.tunits,']'])

subplot(322)
plot(theta, hyp./tt, 'o')
xlabel(str.s36)
ylabel([str.s37, '[',mogs(no).data.cunits,'/',mogs(no).data.tunits,']'])
title(str.s38)

vapp = hyp./(tt-t0(ind));
n = 1:length(ind);
n=n(ind);
ind2 = vapp<0;
if ~isempty(n(ind2))
  disp([str.s40, num2str(n(ind2))])
end

subplot(323)
plot(theta, hyp./(tt-t0(ind)), 'o')
xlabel(str.s36)
ylabel([str.s37, '[',mogs(no).data.cunits,'/',mogs(no).data.tunits,']'])
title(str.s39)

subplot(324)
plot(t0)
xlabel(str.s41)
ylabel([str.s22,' [',mogs(no).data.tunits,']'])
title(str.s42)

subplot(325)
plot(hyp,et,'o')
xlabel(str.s35)
ylabel(str.s43)

subplot(326)
plot(theta,et,'o')
xlabel(str.s36)
ylabel(str.s43)

suptitle(mogs(no).name,'Interpreter','none')

h_stat2 = figure;
set(h_stat2, 'Position',[356 71 512 620]);
drawnow

vapp = hyp./tt;
%lapp = tt./hyp;
Tx = [mogs(no).data.Tx_x(ind)' mogs(no).data.Tx_y(ind)' mogs(no).data.Tx_z(ind)'];
Rx = [mogs(no).data.Rx_x(ind)' mogs(no).data.Rx_y(ind)' mogs(no).data.Rx_z(ind)'];

vmin = min(vapp);
vmax = max(vapp);
c=cmr;

[x0, a]=lsplane([Tx; Rx]);
el = (pi-a(3))*180/pi;
az = atan( cos(a(2))/cos(a(1)) )*180/pi;


m = (size(c,1)-1)/(vmax-vmin);
b = 1-vmin*m;
p = m*vapp(1)+b;
couleur = interp1(c,p);

plot3([Tx(1,1) Rx(1,1)], [Tx(1,2) Rx(1,2)], [Tx(1,3) Rx(1,3)], 'Color',couleur)
hold on
for n=2:length(vapp)
    p = m*vapp(n)+b;
    couleur = interp1(c,p);
    plot3([Tx(n,1) Rx(n,1)], [Tx(n,2) Rx(n,2)], [Tx(n,3) Rx(n,3)], 'Color',couleur)
end
hold off
set(gca, 'DataAspectRatio',[1 1 1],'Units','normalized')
axis tight
view(az,el)
caxis([vmin vmax])
colormap(cmr)
colorbar()
title([mogs(no).name,' - apparent velocity'],'Interpreter','none')
%%%%%%%YH
tlabh = get(gca,'Title');
set(tlabh,'FontSize',12,'Position',get(tlabh,'Position') + [0 0 .4])
%%%%%%%%%%%%%%%%%%%%

function pushbutton_stats_amp_Callback(hObject, eventdata, handles)
no = get(handles.listbox_CDM,'Value');
mogs = getappdata(handles.fig_bh_tomo_db,'mogs');
str = getappdata(handles.fig_bh_tomo_db,'str');

hyp = sqrt( (mogs(no).data.Tx_x-mogs(no).data.Rx_x).^2 + ...
						(mogs(no).data.Tx_y-mogs(no).data.Rx_y).^2 + ...
						(mogs(no).data.Tx_z-mogs(no).data.Rx_z).^2 );
dz = mogs(no).data.Rx_z-mogs(no).data.Tx_z;
theta = 180/pi*asin(dz./hyp);

h_stat = figure;
set(h_stat, 'Position',[256 71 512 620]);

ind = mogs(no).amp_tmax ~= -1 & mogs(no).tauApp ~= -1 & mogs(no).amp_done;
ind = ind & mogs(no).in;

subplot(321)
plot(hyp(ind), mogs(no).tauApp(ind),'o')
xlabel(str.s35)
ylabel('\tau_a')
title(str.s156{2})

subplot(322)
plot(theta(ind), mogs(no).tauApp(ind)./hyp(ind),'o')
xlabel(str.s36)
ylabel('\alpha_a')
title(str.s156{2})

ind = mogs(no).amp_tmax ~= -1 & mogs(no).tauFce ~= -1 & mogs(no).amp_done;
ind = ind & mogs(no).in;

subplot(323)
plot(hyp(ind), mogs(no).tauFce(ind),'o')
xlabel(str.s35)
ylabel('\tau_a')
title(str.s156{3})

subplot(324)
plot(theta(ind), mogs(no).tauFce(ind)./hyp(ind),'o')
xlabel(str.s36)
ylabel('\alpha_a')
title(str.s156{3})

ind = mogs(no).amp_tmax ~= -1 & mogs(no).tauHyb ~= -1 & mogs(no).amp_done;
ind = ind & mogs(no).in;

subplot(325)
plot(hyp(ind), mogs(no).tauHyb(ind),'o')
xlabel(str.s35)
ylabel('\tau_a')
title(str.s156{4})

subplot(326)
plot(theta(ind), mogs(no).tauHyb(ind)./hyp(ind),'o')
xlabel(str.s36)
ylabel('\alpha_a')
title(str.s156{4})



function popupmenu_type_CDM_Callback(hObject, eventdata, handles)
mogs = getappdata(handles.fig_bh_tomo_db,'mogs');
type = get(hObject,'Value');
no = get(handles.listbox_CDM,'Value');
mogs(no).type = type;
setappdata(handles.fig_bh_tomo_db,'mogs',mogs);
update_CDM_coord(handles)
h = getappdata(handles.fig_bh_tomo_db,'h');
h.saved = false;
setappdata(handles.fig_bh_tomo_db,'h',h)


function popupmenu_type_CDM_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


function update_CDM_coord(handles,varargin)
mogs = getappdata(handles.fig_bh_tomo_db,'mogs');
boreholes = getappdata(handles.fig_bh_tomo_db,'boreholes');

updateAll = false;
if nargin >=2
	if islogical( varargin{1} )
		updateAll = varargin{1};
	else
		nos = varargin{1};
	end
end
if exist('nos','var')==0
	if updateAll
		nos = 1:length(mogs);
	else
		nos = get(handles.listbox_CDM,'Value');
	end
end

for no=nos
	if strcmp(mogs(no).data.comment, 'true positions')
        Tx = [mogs(no).data.Tx_x(:) mogs(no).data.Tx_y(:) mogs(no).data.Tx_z(:)];
        mogs(no).TxCosDir = zeros(size(Tx));
        tmp = unique(Tx,'rows');
        tmp = sort(tmp,1,'descend');
        v = -diff(tmp);
        d = sqrt(sum(v.^2,2));
        l = v./kron(d,[1 1 1]);
        l = [l; l(end,:)];
        
        for n=1:size(tmp,1)
            ind = Tx(:,1)==tmp(n,1) & Tx(:,2)==tmp(n,2) & Tx(:,3)==tmp(n,3);
            mogs(no).TxCosDir(ind,1) = l(n,1);
            mogs(no).TxCosDir(ind,2) = l(n,2);
            mogs(no).TxCosDir(ind,3) = l(n,3);
        end
        
        Rx = [mogs(no).data.Rx_x(:) mogs(no).data.Rx_y(:) mogs(no).data.Rx_z(:)];
        mogs(no).RxCosDir = zeros(size(Rx));
        tmp = unique(Rx,'rows');
        tmp = sort(tmp,1,'descend');
        v = -diff(tmp);
        d = sqrt(sum(v.^2,2));
        l = v./kron(d,[1 1 1]);
        l = [l; l(end,:)];
        
        for n=1:size(tmp,1)
            ind = Tx(:,1)==tmp(n,1) & Tx(:,2)==tmp(n,2) & Tx(:,3)==tmp(n,3);
            mogs(no).RxCosDir(ind,1) = l(n,1);
            mogs(no).RxCosDir(ind,2) = l(n,2);
            mogs(no).RxCosDir(ind,3) = l(n,3);
        end
        
		continue
	end
	if isempty( mogs(no).Tx ) || isempty( mogs(no).Rx )
		if length(nos)==1
			return
		else
			continue
		end
	end

	if mogs(no).Tx == mogs(no).Rx
		if length(nos)==1
			uiwait(warndlg(['Tx et Rx are in the same well: ', boreholes(mogs(no).Rx).name]))
			return
		else
			continue
		end
	end

	if mogs(no).type == 1 % cross hole
		mogs(no).data.csurvmod = 'SURVEY MODE        = Trans. - MOG';
		% Tx
		if abs(boreholes(mogs(no).Tx).X-boreholes(mogs(no).Tx).Xmax)<1.0e-5 && ...
				abs(boreholes(mogs(no).Tx).Y-boreholes(mogs(no).Tx).Ymax)<1.0e-5
			% forage vertical
			mogs(no).data.Tx_x(:) = boreholes( mogs(no).Tx ).X;
			mogs(no).data.Tx_y(:) = boreholes( mogs(no).Tx ).Y;
			mogs(no).data.Tx_z = boreholes( mogs(no).Tx ).Z - mogs(no).data.TxOffset - ...
				mogs(no).Tx_z_orig;
            mogs(no).TxCosDir = repmat([0 0 1],mogs(no).data.ntrace,1);
			boreholes( mogs(no).Tx ).fdata = [boreholes( mogs(no).Tx ).X boreholes( mogs(no).Tx ).Y boreholes( mogs(no).Tx ).Z;
				boreholes( mogs(no).Tx ).Xmax boreholes( mogs(no).Tx ).Ymax boreholes( mogs(no).Tx ).Zmax];
		else
			[mogs(no).data.Tx_x, mogs(no).data.Tx_y, mogs(no).data.Tx_z mogs(no).TxCosDir] = ...
				projeteForage(boreholes(mogs(no).Tx).fdata, ...
				mogs(no).Tx_z_orig+mogs(no).data.TxOffset, ['Tx - ',boreholes(mogs(no).Tx).name]);
		end
		% Rx
		if abs(boreholes(mogs(no).Rx).X-boreholes(mogs(no).Rx).Xmax)<1.0e-5 && ...
				abs(boreholes(mogs(no).Rx).Y-boreholes(mogs(no).Rx).Ymax)<1.0e-5
			mogs(no).data.Rx_x(:) = boreholes( mogs(no).Rx ).X;
			mogs(no).data.Rx_y(:) = boreholes( mogs(no).Rx ).Y;
			mogs(no).data.Rx_z = boreholes( mogs(no).Rx ).Z - mogs(no).data.RxOffset - ...
				mogs(no).Rx_z_orig;
            mogs(no).RxCosDir = repmat([0 0 1],mogs(no).data.ntrace,1);
			boreholes( mogs(no).Rx ).fdata = [boreholes( mogs(no).Rx ).X boreholes( mogs(no).Rx ).Y boreholes( mogs(no).Rx ).Z;
				boreholes( mogs(no).Rx ).Xmax boreholes( mogs(no).Rx ).Ymax boreholes( mogs(no).Rx ).Zmax];
		else
			[mogs(no).data.Rx_x, mogs(no).data.Rx_y, mogs(no).data.Rx_z mogs(no).RxCosDir] = ...
				projeteForage(boreholes(mogs(no).Rx).fdata, ...
				mogs(no).Rx_z_orig+mogs(no).data.RxOffset, ['Rx - ',boreholes(mogs(no).Rx).name]);
		end
	else % VRP
		mogs(no).data.csurvmod = 'SURVEY MODE        = Trans. - VRP';
        % Rx
		if abs(boreholes(mogs(no).Rx).X-boreholes(mogs(no).Rx).Xmax)<1.0e-5 || ...
                abs(boreholes(mogs(no).Rx).Y-boreholes(mogs(no).Rx).Ymax)<1.0e-5 %#ok<ALIGN>
    		mogs(no).data.Rx_x(:) = boreholes( mogs(no).Rx ).X;
        	mogs(no).data.Rx_y(:) = boreholes( mogs(no).Rx ).Y;
            mogs(no).data.Rx_z = boreholes( mogs(no).Rx ).Z - mogs(no).data.RxOffset - ...
                mogs(no).Rx_z_orig;
            mogs(no).RxCosDir = repmat([0 0 1],mogs(no).data.ntrace,1);
			boreholes( mogs(no).Rx ).fdata = [boreholes( mogs(no).Rx ).X boreholes( mogs(no).Rx ).Y boreholes( mogs(no).Rx ).Z;
                boreholes( mogs(no).Rx ).Xmax boreholes( mogs(no).Rx ).Ymax boreholes( mogs(no).Rx ).Zmax];
		else
			[mogs(no).data.Rx_x, mogs(no).data.Rx_y, mogs(no).data.Rx_z mogs(no).RxCosDir] = ...
				projeteForage(boreholes(mogs(no).Rx).fdata, ...
				mogs(no).Rx_z_orig+mogs(no).data.RxOffset, ['Rx - ',boreholes(mogs(no).Rx).name]);
        end
        % Tx en surface
		theta = atan2( boreholes( mogs(no).Tx ).Y-boreholes( mogs(no).Rx ).Y,...
			boreholes( mogs(no).Tx ).X-boreholes( mogs(no).Rx ).X );
		mogs(no).data.Tx_x = boreholes( mogs(no).Rx ).X + mogs(no).Tx_z_orig*cos(theta);
		mogs(no).data.Tx_y = boreholes( mogs(no).Rx ).Y + mogs(no).Tx_z_orig*sin(theta);
		% z -> on assume que z varie lineairement entre les deux trous
		l = sqrt( (boreholes( mogs(no).Tx ).Y-boreholes( mogs(no).Rx ).Y)^2 + ...
			(boreholes( mogs(no).Tx ).X-boreholes( mogs(no).Rx ).X)^2 );
		dz = boreholes( mogs(no).Tx ).Z_surf-boreholes( mogs(no).Rx ).Z_surf;
		mogs(no).data.Tx_z = boreholes( mogs(no).Rx ).Z_surf + dz*mogs(no).Tx_z_orig/l;
        
        boreholes( mogs(no).Tx ).fdata = [mogs(no).data.Tx_x(1) mogs(no).data.Tx_y(1) mogs(no).data.Tx_z(1); ...
            mogs(no).data.Tx_x(end) mogs(no).data.Tx_y(end) mogs(no).data.Tx_z(end)];
        
        d = sqrt(sum((boreholes( mogs(no).Tx ).fdata(2,1:3)-boreholes( mogs(no).Tx ).fdata(1,1:3)).^2));
        % cosinus directeurs
        l = (boreholes( mogs(no).Tx ).fdata(2,1:3)-boreholes( mogs(no).Tx ).fdata(1,1:3))./d;
        mogs(no).TxCosDir = repmat(l,mogs(no).data.ntrace,1);
	end
end

setappdata(handles.fig_bh_tomo_db,'mogs',mogs);



function pushbutton_export_tt_Callback(hObject, eventdata, handles)
str=getappdata(handles.fig_bh_tomo_db,'str');
[filename, pathname] = uiputfile({'*.dat';'*.*'}, str.s74);
if isequal(filename,0) || isequal(pathname,0)
	return;
end
no = get(handles.listbox_CDM,'Value');
mogs = getappdata(handles.fig_bh_tomo_db,'mogs');
air = getappdata(handles.fig_bh_tomo_db,'air');
%str = getappdata(handles.fig_bh_tomo_db,'str');

tt = mogs(no).tt;
et = mogs(no).et;
done = mogs(no).tt_done & mogs(no).in;

ind = tt~=-1 & mogs(no).in;

if sum(done)==0
	uiwait(errordlg('Data not processed!!!'))
	return
end

if strcmp( mogs(no).data.tunits,'ns' )==1
    av = air( mogs(no).av );
    ap = air( mogs(no).ap );
    if mogs(no).data.synthetique==1
        t0 = zeros(1,length(mogs(no).tt));
        fac_dt_av = 1;
        fac_dt_ap = 1;
    else
        [t0, fac_dt_av, fac_dt_ap] = corr_t0(length(mogs(no).tt), av, ap, true);
    end
    if ~isempty(av), air( mogs(no).av ).fac_dt = fac_dt_av; end
    if ~isempty(ap), air( mogs(no).ap ).fac_dt = fac_dt_ap; end
    if mogs(no).user_fac_dt==0
        if fac_dt_av~=1 && fac_dt_ap ~= 1
            mogs(no).fac_dt = 0.5*(fac_dt_av+fac_dt_ap);
        elseif fac_dt_av~=1
            mogs(no).fac_dt = fac_dt_av;
        elseif fac_dt_ap~=1
            mogs(no).fac_dt = fac_dt_ap;
        else
            mogs(no).fac_dt = 1;
        end
        set(handles.edit_fac_dt,'String',num2str(mogs(no).fac_dt));
    end
end
setappdata(handles.fig_bh_tomo_db,'air', air);
setappdata(handles.fig_bh_tomo_db,'mogs', mogs);


tt = mogs(no).fac_dt*(tt(ind)-t0(ind));
et = et(ind);

data = [mogs(no).data.Tx_x(ind)' mogs(no).data.Tx_y(ind)' mogs(no).data.Tx_z(ind)' ...
	mogs(no).data.Rx_x(ind)' mogs(no).data.Rx_y(ind)' mogs(no).data.Rx_z(ind)' ...
	tt' et' mogs(no).no_traces(ind)']; %#ok<NASGU>

eval(['save ',pathname,filename,' data -ascii'])


function pushbutton_vel_con_Callback(hObject, eventdata, handles)
str=getappdata(handles.fig_bh_tomo_db,'str');
[file, rep] = uigetfile('*.con',[str.s252,' - ',str.s253]);
if file==0
    return
end
cont = load([rep,file]);
boreholes = getappdata(handles.fig_bh_tomo_db,'boreholes');
no = get(handles.listbox_forages,'Value');
scont.x = boreholes(no).X;    % boreholes droits
scont.y = boreholes(no).Y;
scont.z = boreholes(no).Z - cont(:,1);
scont.valeur = cont(:,2);
if size(cont,2)==3
	scont.variance = cont(:,3);
else
	scont.variance = zeros(size(cont(:,2)));
end
boreholes(no).scont = scont;
setappdata(handles.fig_bh_tomo_db,'boreholes',boreholes);
h = getappdata(handles.fig_bh_tomo_db,'h');
h.saved = false;
setappdata(handles.fig_bh_tomo_db,'h',h)


function pushbutton_att_con_Callback(hObject, eventdata, handles)
str=getappdata(handles.fig_bh_tomo_db,'str');
[file, rep] = uigetfile('*.con',[str.s252,' - ',str.s177]);
if file==0
    return
end
cont = load([rep,file]);
boreholes = getappdata(handles.fig_bh_tomo_db,'boreholes');
no = get(handles.listbox_forages,'Value');
acont.x = boreholes(no).X;    % boreholes droits
acont.y = boreholes(no).Y;
acont.z = boreholes(no).Z - cont(:,1);
acont.valeur = cont(:,2);
if size(cont,2)==3
	acont.variance = cont(:,3);
else
	acont.variance = zeros(size(cont(:,2)));
end
boreholes(no).acont = acont;
setappdata(handles.fig_bh_tomo_db,'boreholes',boreholes);
h = getappdata(handles.fig_bh_tomo_db,'h');
h.saved = false;
setappdata(handles.fig_bh_tomo_db,'h',h)


function pushbutton_export_tau_Callback(hObject, eventdata, handles)
str=getappdata(handles.fig_bh_tomo_db,'str');
[filename, pathname] = uiputfile({'*.dat';'*.*'}, str.s74);
if isequal(filename,0) || isequal(pathname,0)
	return;
end
no = get(handles.listbox_CDM,'Value');
mogs = getappdata(handles.fig_bh_tomo_db,'mogs');
%air = getappdata(handles.fig_bh_tomo_db,'air');
%str = getappdata(handles.fig_bh_tomo_db,'str');

tt = mogs(no).tauApp;
et = mogs(no).tauApp_et;
%done = mogs(no).amp_done;

ind = tt~=-1 & mogs(no).in;
tt = tt(ind);
et = et(ind);

if isempty(tt)
	uiwait(errordlg('Data not processed!!!'))
	return
end

data = [mogs(no).data.Tx_x(ind)' mogs(no).data.Tx_y(ind)' mogs(no).data.Tx_z(ind)' ...
	mogs(no).data.Rx_x(ind)' mogs(no).data.Rx_y(ind)' mogs(no).data.Rx_z(ind)' ...
	tt' et' mogs(no).no_traces(ind)']; %#ok<NASGU>

eval(['save ',pathname,filename,' data -ascii'])


function pushbutton_renameCDM_Callback(hObject, eventdata, handles)
no = get(handles.listbox_CDM,'Value');
mogs = getappdata(handles.fig_bh_tomo_db,'mogs');
names_mog = getappdata(handles.fig_bh_tomo_db,'names_mog');
old{1} = names_mog{no};
answer = inputdlg('New name','',1,old);
if isempty(answer)
    return
end
if isempty(answer{1})
    return
end

mogs(no).name = answer{1};
names_mog{no} = answer{1};
if length(mogs)==1
   	set(handles.listbox_CDM,'Value',no,'String',names_mog{1})
else
	set(handles.listbox_CDM,'Value',no,'String',names_mog)
end
    
setappdata(handles.fig_bh_tomo_db,'mogs',mogs);
setappdata(handles.fig_bh_tomo_db,'names_mog',names_mog);
h = getappdata(handles.fig_bh_tomo_db,'h');
h.saved = false;
setappdata(handles.fig_bh_tomo_db,'h',h)


function pushbutton_couverture_Callback(hObject, eventdata, handles)
no = get(handles.listbox_CDM,'Value');
mogs = getappdata(handles.fig_bh_tomo_db,'mogs');
boreholes = getappdata(handles.fig_bh_tomo_db,'boreholes');

ind1 = mogs(no).tt==-1 & mogs(no).in;
ind2 = mogs(no).tt~=-1 & mogs(no).in;

[~, a]=lsplane([boreholes(mogs(no).Tx).fdata; boreholes(mogs(no).Rx).fdata]);
el = (pi-a(3))*180/pi - 15;
az = atan( cos(a(2))/cos(a(1)) )*180/pi - 10;

figure
plot3([mogs(no).data.Tx_x(ind1); mogs(no).data.Rx_x(ind1)], ...
	[mogs(no).data.Tx_y(ind1); mogs(no).data.Rx_y(ind1)], ...
	[mogs(no).data.Tx_z(ind1); mogs(no).data.Rx_z(ind1)],'r')
hold on
plot3([mogs(no).data.Tx_x(ind2); mogs(no).data.Rx_x(ind2)], ...
	[mogs(no).data.Tx_y(ind2); mogs(no).data.Rx_y(ind2)], ...
	[mogs(no).data.Tx_z(ind2); mogs(no).data.Rx_z(ind2)],'g')
hold off
set(gca, 'DataAspectRatio',[1 1 1],'Units','normalized')
title([num2str(100*sum(ind2)/length(ind2)), '%'])
axis tight
xlabel('Tx-Rx Distance [m]','FontSize',8)
ylabel('Elevation [m]','FontSize',8)
view(az,el)
%%%%%%%YH
set(gca,'FontSize',8);
%colorbar('peer',gca);
% xlabh = get(gca,'XLabel');
% set(xlabh,'Position',get(xlabh,'Position') + [-2 0 -.7])
% ylabh = get(gca,'YLabel');
% set(ylabh,'Position',get(ylabh,'Position') + [2 0 -.7])
tlabh = get(gca,'Title');
set(tlabh,'FontSize',8,'Position',get(tlabh,'Position') + [0 0 .3])
%%%%%%%%%%%%%%%%%%%%

function pushbutton_show_raw_Callback(hObject, eventdata, handles)
no = get(handles.listbox_CDM,'Value');
mogs = getappdata(handles.fig_bh_tomo_db,'mogs');
in = mogs(no).in;

figure
imagesc(1:mogs(no).no_traces(in),mogs(no).data.timestp,mogs(no).data.rdata(:,in))
clim = caxis;
cmax = max(abs(clim));
caxis([-cmax cmax])
colormap ramac_cmap(16)
colorbar('peer',gca);
xlabel('Trace no')
ylabel('Time [ns]')

function pushbutton_importerCDM_Callback(hObject, eventdata, handles)
% h = getappdata(handles.fig_bh_tomo_db,'h');  %%YH
% [no_mog,db_file,hh] = choisirMOG('UserData', h.db_file);  %%YH
[no_mog,db_file,hh] = choisirMOG();
delete(hh)
if no_mog==0
	return
end
tmp = load(db_file,'mogs','boreholes','air');
mog = tmp.mogs(no_mog);

mogs = getappdata(handles.fig_bh_tomo_db,'mogs');
names_mog = getappdata(handles.fig_bh_tomo_db,'names_mog');
air = getappdata(handles.fig_bh_tomo_db,'air');
boreholes = getappdata(handles.fig_bh_tomo_db,'boreholes');

for n=1:length(names_mog)
	if strcmp(names_mog{n}, mog.name) %#ok<ALIGN>
		loc_str = getappdata(handles.fig_bh_tomo_db,'str');
		ButtonName=questdlg(loc_str.s147);
		switch lower(ButtonName),
			case 'yes',
				mog.name=[mog.name,'_2'];
				break
			otherwise
				return
		end
    end
end

for n=1:length(air)
	if strcmp(air(n).name, tmp.air(mog.av).name)
		tmp.air(mog.av).name = [tmp.air(mog.av).name,'_2'];
	end
	if strcmp(air(n).name, tmp.air(mog.ap).name)
		tmp.air(mog.ap).name = [tmp.air(mog.ap).name,'_2'];
	end
end

air = [air tmp.air(mog.av)];
mog.av = length(air);
air = [air tmp.air(mog.ap)];
mog.ap = length(air);

addTx = true;
addRx = true;
mogTx = mog.Tx;
mogRx = mog.Rx;
for n=1:length(boreholes)
	if strcmp(boreholes(n).name, tmp.boreholes(mogTx).name)
		addTx = false;
		mog.Tx = n;
	end
	if strcmp(boreholes(n).name, tmp.boreholes(mogRx).name)
		addRx = false;
		mog.Rx = n;
	end
end

if addTx
	if ~isfield(tmp.boreholes(mog.Tx),'fdata')
		tmp.boreholes(mog.Tx).fdata = [tmp.boreholes(mog.Tx).X tmp.boreholes(mog.Tx).Y tmp.boreholes(mog.Tx).Z; ...
				tmp.boreholes(mog.Tx).Xmax tmp.boreholes(mog.Tx).Ymax tmp.boreholes(mog.Tx).Zmax];
	end
	boreholes = [boreholes tmp.boreholes(mog.Tx)];
	mog.Tx = length(boreholes);
end
if addRx
	if ~isfield(tmp.boreholes(mog.Rx),'fdata')
		tmp.boreholes(mog.Rx).fdata = [tmp.boreholes(mog.Rx).X tmp.boreholes(mog.Rx).Y tmp.boreholes(mog.Rx).Z; ...
				tmp.boreholes(mog.Rx).Xmax tmp.boreholes(mog.Rx).Ymax tmp.boreholes(mog.Rx).Zmax];
	end
	boreholes = [boreholes tmp.boreholes(mog.Rx)];
	mog.Rx = length(boreholes);
end

if isempty(mogs)
	mog.no_traces = 1:mog.data.ntrace;
else
	mog.no_traces = mogs(end).no_traces(end) + (1:mog.data.ntrace);
end
if ~isfield(mog,'fac_dt')
	mog.fac_dt = 1;
	mog.user_fac_dt = 0;
end
	
mogs = [mogs mog];
names_mog{end+1} = mog.name;

setappdata(handles.fig_bh_tomo_db,'boreholes',boreholes)
setappdata(handles.fig_bh_tomo_db,'mogs',mogs)
setappdata(handles.fig_bh_tomo_db,'names_mog',names_mog)
setappdata(handles.fig_bh_tomo_db,'air',air)
set(handles.edit_Tx_offset,'String',num2str(mogs(end).data.TxOffset))
set(handles.edit_Rx_offset,'String',num2str(mogs(end).data.RxOffset))
set(handles.edit_fac_dt,'String',num2str(mogs(end).fac_dt))
set(handles.checkbox_fac_dt,'Value',mogs(end).user_fac_dt)
set(handles.edit_f_mul_et,'String',num2str(mogs(end).f_et))

update_CDM_coord(handles, length(mogs))
update_tout(hObject, eventdata,handles)  %YH

h = getappdata(handles.fig_bh_tomo_db,'h');
h.saved = false;
setappdata(handles.fig_bh_tomo_db,'h',h)


function pushbutton_plotBH_Callback(hObject, eventdata, handles)
boreholes = getappdata(handles.fig_bh_tomo_db,'boreholes');
if isempty(boreholes), return, end

figure
for n=1:length(boreholes)
	plot3([boreholes(n).X boreholes(n).Xmax],...
		[boreholes(n).Y boreholes(n).Ymax],...
		[boreholes(n).Z boreholes(n).Zmax],'g','LineWidth',2)
	hold on
	plot3(boreholes(n).X, boreholes(n).Y, boreholes(n).Z_surf,'ro')
	text(boreholes(n).X, boreholes(n).Y, boreholes(n).Z_surf, boreholes(n).name)
end
grid on
hold off
xlabel('X [m]')
ylabel('Y [m]')
zlabel('Elevation [m]')
set(gca,'DataAspectRatio',[1 1 1])


function pushbutton_importerBH_Callback(hObject, eventdata, handles)
loc_str = getappdata(handles.fig_bh_tomo_db,'str');

[filename, pathname] = uigetfile('*.xyz', loc_str.s266);
if isequal(filename,0) || isequal(pathname,0)
	return
else
	name = filename(1:end-4);
	fdata = load([pathname,filename]);
end

boreholes = getappdata(handles.fig_bh_tomo_db,'boreholes');
no = length(boreholes)+1;

boreholes(no).name = char( name );
boreholes(no).X = fdata(1,1);
boreholes(no).Y = fdata(1,2);
boreholes(no).Z = fdata(1,3);
boreholes(no).Xmax = fdata(end,1);
boreholes(no).Ymax = fdata(end,2);
boreholes(no).Zmax = fdata(end,3);
boreholes(no).Z_surf = 0;
boreholes(no).scont = [];
boreholes(no).acont = [];
boreholes(no).fdata = fdata;
boreholes(no).Z_water = NaN;
boreholes(no).diam = 0;

set(handles.edit_forage_X,'String',num2str( boreholes(no).X ));
set(handles.edit_forage_Y,'String',num2str( boreholes(no).Y ));
set(handles.edit_forage_Z,'String',num2str( boreholes(no).Z ));
set(handles.edit_forage_Xmax,'String',num2str( boreholes(no).Xmax ));
set(handles.edit_forage_Ymax,'String',num2str( boreholes(no).Ymax ));
set(handles.edit_forage_Zmax,'String',num2str( boreholes(no).Zmax ));
set(handles.edit_forage_surface,'String',num2str( boreholes(no).Z_surf ));
set(handles.edit_forage_eau,'String',num2str( boreholes(no).Z_water ));
set(handles.edit_forage_diam,'String',num2str( boreholes(no).diam ));
str = get(handles.listbox_forages,'String');
str{no} = char( boreholes(no).name );
set(handles.listbox_forages,'String',str)
set(handles.listbox_forages,'Value',no)
set(handles.popupmenu_TxCDM,'String',str)
set(handles.popupmenu_RxCDM,'String',str)
h = getappdata(handles.fig_bh_tomo_db,'h');
h.saved = false;
setappdata(handles.fig_bh_tomo_db,'h',h)
setappdata(handles.fig_bh_tomo_db,'boreholes',boreholes);


function edit_fac_dt_Callback(hObject, eventdata, handles)
val = str2double(get(hObject,'String'));
if isnan(val)
    str = getappdata(handles.fig_bh_tomo_db,'str');
    errordlg(str.s54,str.s45,'modal')
    return
end
no = get(handles.listbox_CDM,'Value');
mogs = getappdata(handles.fig_bh_tomo_db,'mogs');

mogs(no).fac_dt = val;

setappdata(handles.fig_bh_tomo_db,'mogs',mogs)
h = getappdata(handles.fig_bh_tomo_db,'h');
h.saved = false;
setappdata(handles.fig_bh_tomo_db,'h',h)


function edit_fac_dt_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


function edit_Tx_offset_Callback(hObject, eventdata, handles)
no = get(handles.listbox_CDM,'Value');
mogs = getappdata(handles.fig_bh_tomo_db,'mogs');
offset = str2double(get(hObject,'String'));
if isnan(offset)
    str = getappdata(handles.fig_bh_tomo_db,'str');
    errordlg(str.s54,str.s45,'modal')
    return
end
mogs(no).data.TxOffset = offset;
setappdata(handles.fig_bh_tomo_db,'mogs',mogs)
if mogs(no).Tx ~= mogs(no).Rx
	update_CDM_coord(handles,no)
end
h = getappdata(handles.fig_bh_tomo_db,'h');
h.saved = false;
setappdata(handles.fig_bh_tomo_db,'h',h)


function edit_Tx_offset_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


function edit_Rx_offset_Callback(hObject, eventdata, handles)
no = get(handles.listbox_CDM,'Value');
mogs = getappdata(handles.fig_bh_tomo_db,'mogs');
offset = str2double(get(hObject,'String'));
if isnan(offset)
    str = getappdata(handles.fig_bh_tomo_db,'str');
    errordlg(str.s54,str.s45,'modal')
    return
end
mogs(no).data.RxOffset = offset;
setappdata(handles.fig_bh_tomo_db,'mogs',mogs)
if mogs(no).Tx ~= mogs(no).Rx
	update_CDM_coord(handles,no)
end
h = getappdata(handles.fig_bh_tomo_db,'h');
h.saved = false;
setappdata(handles.fig_bh_tomo_db,'h',h)


function edit_Rx_offset_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


function pushbutton_elaguerMOG_Callback(hObject, eventdata, handles)
no = get(handles.listbox_CDM,'Value');
mogs = getappdata(handles.fig_bh_tomo_db,'mogs');
par_decime = getappdata(handles.fig_bh_tomo_db,'par_decime');

[hh, in, par] = bh_tomo_decime(mogs(no), par_decime(no));
if ishandle(hh), delete(hh), end
mogs(no).in = in;
par_decime(no) = par;
setappdata(handles.fig_bh_tomo_db,'mogs',mogs)
setappdata(handles.fig_bh_tomo_db,'par_decime',par_decime)


function edit_fnom_Callback(hObject, eventdata, handles)
no = get(handles.listbox_CDM,'Value');
mogs = getappdata(handles.fig_bh_tomo_db,'mogs');
f = str2double(get(hObject,'String'));
mogs(no).data.rnomfreq = f;
setappdata(handles.fig_bh_tomo_db,'mogs',mogs)


function edit_fnom_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


function pushbutton_spectra_Callback(hObject, eventdata, handles)
no = get(handles.listbox_CDM,'Value');
mogs = getappdata(handles.fig_bh_tomo_db,'mogs');
bh_tomo_spectra(mogs(no))


function checkbox_fac_dt_Callback(hObject, eventdata, handles)
no = get(handles.listbox_CDM,'Value');
mogs = getappdata(handles.fig_bh_tomo_db,'mogs');
f = get(hObject,'Value');
mogs(no).user_fac_dt = f;
setappdata(handles.fig_bh_tomo_db,'mogs',mogs)


function edit_date_Callback(hObject, eventdata, handles)
no = get(handles.listbox_CDM,'Value');
mogs = getappdata(handles.fig_bh_tomo_db,'mogs');
d = get(hObject,'String');
mogs(no).date = d;
setappdata(handles.fig_bh_tomo_db,'mogs',mogs)


function edit_date_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% YH
function pushbutton_add_model3d_Callback(hObject, eventdata, handles)
model3d = getappdata(handles.fig_bh_tomo_db, 'model3d');
if ~isempty(model3d)
    n = length(model3d)+1;
else
    n = 1;
end
loc_str = getappdata(handles.fig_bh_tomo_db,'str');
prompt = loc_str.s295;
nameDefault = ['3dmodel',num2str(n)]; 
names = get(handles.listbox_model3d,'String');
for i =1:length(names)
    if strcmp(names(i),nameDefault)
        nameDefault = ['3dmodel',num2str(n+1)];
        break;
    end
end
title = 'Add 3d model';
nblines = 1;
answer = inputdlg(prompt,title,nblines,{nameDefault});
name = [];
if ~isempty(answer)
    name=answer{1};
end

if isempty(name)
	return;
else
    flag=0;
	for i=1:length(model3d)
		if strcmp(model3d(i).name, name)
			flag = 1;
			break;
		end
	end
	if flag==1
		rep=questdlg(['Do you want to replace ',name,'?']);
		if ~strcmp(rep,'Yes')
			return
        end
    end  
end
model3d(n).name = char( name );
model3d(n).mogs = [];
model3d(n).boreholes = [];
model3d(n).type = 'normal';

for i = 1:n
	str{i} = char( model3d(i).name );
end
set(handles.listbox_model3d,'String',str)
set(handles.listbox_model3d,'Value',n)
set(handles.listbox_CDM_model3d,'String','')
setappdata(handles.fig_bh_tomo_db,'model3d',model3d);
update_info(handles)
h = getappdata(handles.fig_bh_tomo_db,'h');
h.saved = false;
setappdata(handles.fig_bh_tomo_db,'h',h)


function pushbutton_remove_model3d_Callback(hObject, eventdata, handles)
model3d = getappdata(handles.fig_bh_tomo_db,'model3d');
no = get(handles.listbox_model3d,'Value');
ind=1:length(model3d);
ind = ind~=no;
model3d = model3d(ind);
liste = cell(0);
nn=1;
for n=1:length(model3d)
	liste{nn} = char( model3d(n).name );
	nn=nn+1;
end
set(handles.listbox_model3d,'String',liste,'Value',length(liste))
setappdata(handles.fig_bh_tomo_db,'model3d',model3d)
no=no-1;
if isempty(model3d)
    set(handles.listbox_CDM_model3d,'String','')
    return
elseif no==0
    no=1;
end
names_mog = getappdata(handles.fig_bh_tomo_db,'names_mog');
str = cell(1,length(model3d(no).mogs));
for n=1:length(model3d(no).mogs)
	str{n} = char( names_mog{model3d(no).mogs(n) } );
end
set(handles.listbox_CDM_model3d,'Value',no)
set(handles.listbox_CDM_model3d,'String',str)
h = getappdata(handles.fig_bh_tomo_db,'h');
h.saved = false;
setappdata(handles.fig_bh_tomo_db,'h',h)


function listbox_model3d_Callback(hObject, eventdata, handles)
model3d = getappdata(handles.fig_bh_tomo_db,'model3d');
names_mog = getappdata(handles.fig_bh_tomo_db,'names_mog');

no = get(handles.listbox_model3d,'Value'); 
str = cell(1,length(model3d(no).mogs));
for n=1:length(model3d(no).mogs)
	str{n} = char( names_mog{model3d(no).mogs(n) } );
end
set(handles.listbox_CDM_model3d,'Value',1)
set(handles.listbox_CDM_model3d,'String',str)


function listbox_model3d_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function tmp = prepare_grid3d_data(handles, no)
boreholes = getappdata(handles.fig_bh_tomo_db,'boreholes');
mogs = getappdata(handles.fig_bh_tomo_db,'mogs');
model3d = getappdata(handles.fig_bh_tomo_db,'model3d');
if isempty(model3d(no).mogs)
    msgbox('Please add mogs for this panel.','Message','modal');
    return;
end
tmp.boreholes = [];
for n=1:length(model3d(no).boreholes)
    tmp.boreholes = [tmp.boreholes boreholes(model3d(no).boreholes(n))];
end
tmp.Tx = [];
tmp.Rx = [];
tmp.TxCosDir = [];
tmp.RxCosDir = [];
tmp.Tx_Z_water = [];
tmp.Rx_Z_water = [];
tmp.in = [];
for n=1:length(model3d(no).mogs)
    tmp.in = [tmp.in; mogs(model3d(no).mogs(n)).in'];
    tmp.Tx = [tmp.Tx; [mogs(model3d(no).mogs(n)).data.Tx_x' ...
        mogs(model3d(no).mogs(n)).data.Tx_y' ...
        mogs(model3d(no).mogs(n)).data.Tx_z'] ];
    tmp.Rx = [tmp.Rx; [mogs(model3d(no).mogs(n)).data.Rx_x' ...
        mogs(model3d(no).mogs(n)).data.Rx_y' ...
        mogs(model3d(no).mogs(n)).data.Rx_z'] ];
    
    tmp.TxCosDir = [tmp.TxCosDir; mogs(model3d(no).mogs(n)).TxCosDir];
    tmp.RxCosDir = [tmp.RxCosDir; mogs(model3d(no).mogs(n)).RxCosDir];
    
    trou = boreholes(mogs(model3d(no).mogs(n)).Tx);
    
    if ~isnan(trou.Z_water) 
        x = interp1(trou.fdata(:,3), trou.fdata(:,1), trou.Z_water);
        y = interp1(trou.fdata(:,3), trou.fdata(:,2), trou.Z_water);
    else
        x = NaN;
        y = NaN;
    end
    tmp.Tx_Z_water = [tmp.Tx_Z_water; [x y trou.Z_water] ];
        
    trou = boreholes(mogs(model3d(no).mogs(n)).Rx);
    if ~isnan(trou.Z_water)
        x = interp1(trou.fdata(:,3), trou.fdata(:,1), trou.Z_water);
        y = interp1(trou.fdata(:,3), trou.fdata(:,2), trou.Z_water);
    else
        x = NaN;
        y = NaN;
    end
    tmp.Rx_Z_water = [tmp.Rx_Z_water; [x y trou.Z_water] ];
    
end
tmp.in = logical(tmp.in);


function pushbutton_create_grid_model3d_Callback(hObject, eventdata, handles)
no = get(handles.listbox_model3d,'Value');
tmp = prepare_grid3d_data(handles, no);
[hh,grid3d] = bh_tomo_grid3d(tmp);
if ishandle(hh), delete(hh), end
model3d = getappdata(handles.fig_bh_tomo_db,'model3d');

if ~isempty(grid3d.grx)
    model3d(no).grid3d = grid3d;
    setappdata(handles.fig_bh_tomo_db,'model3d', model3d)
end
h = getappdata(handles.fig_bh_tomo_db,'h');
h.saved = false;
setappdata(handles.fig_bh_tomo_db,'h',h)


function pushbutton_edit_grid_model3d_Callback(hObject, eventdata, handles)
model3d = getappdata(handles.fig_bh_tomo_db,'model3d');
no = get(handles.listbox_model3d,'value');
tmp = prepare_grid3d_data(handles, no);
if isfield(model3d(no), 'grid3d') && ~isempty(model3d(no).grid3d)  
    [hh,grid3d] = bh_tomo_grid3d(tmp,model3d(no).grid3d);
else
    [hh,grid3d] = bh_tomo_grid3d(tmp);
end
if ishandle(hh), delete(hh), end
if ~isempty(grid3d.grx)
    model3d(no).grid3d = grid3d;
    setappdata(handles.fig_bh_tomo_db,'model3d', model3d)
end
h = getappdata(handles.fig_bh_tomo_db,'h');
h.saved = false;
setappdata(handles.fig_bh_tomo_db,'h',h)


function pushbutton_add_CDM_model3d_Callback(hObject, eventdata, handles)
model3d = getappdata(handles.fig_bh_tomo_db,'model3d');
mogs = getappdata(handles.fig_bh_tomo_db,'mogs');
no = get(handles.listbox_model3d,'Value');

names_mog = getappdata(handles.fig_bh_tomo_db,'names_mog');
mogs = getappdata(handles.fig_bh_tomo_db,'mogs');

[no_mog,h] = ajouteMOG(names_mog);
delete(h)
if no_mog==0
    return
else
    flag=0;
	for i=1:length(model3d(no).mogs)
		if model3d(no).mogs(i)== no_mog
			flag = 1;
			break;
		end
	end
	if flag==1
		return;
    end  
end

nn = 1+length(model3d(no).mogs);
if nn==1
    model3d(no).boreholes(nn) = mogs(no_mog).Tx;
    model3d(no).boreholes(nn+1) = mogs(no_mog).Rx;
else
    ii=find(model3d(no).boreholes==mogs(no_mog).Tx,1);
    if isempty(ii) %#ok<ALIGN>
		% append to list
        model3d(no).boreholes = [model3d(no).boreholes mogs(no_mog).Tx];
	end
    ii=find(model3d(no).boreholes==mogs(no_mog).Rx,1);
    if isempty(ii) %#ok<ALIGN>
		% append to list
        model3d(no).boreholes = [model3d(no).boreholes mogs(no_mog).Rx];
	end
end
model3d(no).mogs(nn) = no_mog;
%
% update grid if it exists
%
str = cell(1,length(model3d(no).mogs));
for n=1:length(model3d(no).mogs)
	str{n} = char( names_mog{model3d(no).mogs(n) } );
end
set(handles.listbox_CDM_model3d,'String',str)
set(handles.listbox_CDM_model3d,'Value',nn)
setappdata(handles.fig_bh_tomo_db,'model3d',model3d);
h = getappdata(handles.fig_bh_tomo_db,'h');
h.saved = false;
setappdata(handles.fig_bh_tomo_db,'h',h)


function pushbutton_remove_CDM_model3d_Callback(hObject, eventdata, handles)
model3d = getappdata(handles.fig_bh_tomo_db,'model3d');
no = get(handles.listbox_model3d,'Value');
no_CDM = get(handles.listbox_CDM_model3d,'Value');
ind=1:length(model3d(no).mogs);
ind = ind~=no_CDM;
model3d(no).mogs = model3d(no).mogs(ind);
setappdata(handles.fig_bh_tomo_db,'model3d',model3d);
str = get(handles.listbox_CDM_model3d,'String');
nn=1;
str2=cell(0);
for n=1:length(ind)
	if ind(n)
		str2{nn} = str{n};
		nn=nn+1;
	end
end
no_CDM = no_CDM-1;
if no_CDM==0
	no_CDM=1;
end
set(handles.listbox_CDM_model3d,'Value',no_CDM,'String',str2);
h = getappdata(handles.fig_bh_tomo_db,'h');
h.saved = false;
setappdata(handles.fig_bh_tomo_db,'h',h)



function listbox_CDM_model3d_Callback(hObject, eventdata, handles)


function listbox_CDM_model3d_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


function fig_bh_tomo_db_CloseRequestFcn(hObject, eventdata, handles)
delete(hObject);
