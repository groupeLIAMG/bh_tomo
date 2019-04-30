function varargout = bh_tomo_phase_pick(varargin)
% BH_TOMO_PHASE_PICK M-file for bh_tomo_phase_pick.fig
%      BH_TOMO_PHASE_PICK, by itself, creates a new BH_TOMO_PHASE_PICK or raises the existing
%      singleton*.
%
%      H = BH_TOMO_PHASE_PICK returns the handle to a new BH_TOMO_PHASE_PICK or the handle to
%      the existing singleton*.
%
%      BH_TOMO_PHASE_PICK('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in BH_TOMO_PHASE_PICK.M with the given input arguments.
%
%      BH_TOMO_PHASE_PICK('Property','Value',...) creates a new BH_TOMO_PHASE_PICK or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before bh_tomo_phase_pick_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to bh_tomo_phase_pick_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help bh_tomo_phase_pick

% Last Modified by GUIDE v2.5 17-Dec-2012 13:51:14

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @bh_tomo_phase_pick_OpeningFcn, ...
                   'gui_OutputFcn',  @bh_tomo_phase_pick_OutputFcn, ...
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


% --- Executes just before bh_tomo_phase_pick is made visible.
function bh_tomo_phase_pick_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to bh_tomo_phase_pick (see VARARGIN)

% Choose default command line output for bh_tomo_phase_pick
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

h.db_file = get(handles.fig_bh_tomo_phase_pick,'UserData');

%h.str = get_str_locale();
setappdata(handles.fig_bh_tomo_phase_pick, 'h', h)

% UIWAIT makes bh_tomo_phase_pick wait for user response (see UIRESUME)
% uiwait(handles.fig_bh_tomo_phase_pick);


% --- Outputs from this function are returned to the command line.
function varargout = bh_tomo_phase_pick_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


function edit_SNR_threshold_Callback(hObject, eventdata, handles)


function edit_SNR_threshold_CreateFcn(hObject, eventdata, handles)

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


function checkbox_denoise_Callback(hObject, eventdata, handles)


function pushbutton_go_Callback(hObject, eventdata, handles)
mog = getappdata(handles.fig_bh_tomo_phase_pick, 'mog');
if isempty(mog)
	str = get_str_locale();
	warndlg(str.s52)
	return
end

par.trace_cut = 300;
par.L = str2double(get(handles.edit_L,'String'));
par.medfilt_lth = str2double(get(handles.edit_medfilt_lth,'String'));
par.TH = str2double(get(handles.edit_SNR_threshold,'String'));
par.denoise = logical(get(handles.checkbox_denoise,'Value'));
par.AIC_corr = false;
par.inv_polarite_cwt = logical(get(handles.checkbox_inv_polarite,'Value'));

pick = phase_picking(mog, par);

setappdata(handles.fig_bh_tomo_phase_pick, 'pick', pick)
h = getappdata(handles.fig_bh_tomo_phase_pick, 'h');
h.noTx = 1;
setappdata(handles.fig_bh_tomo_phase_pick, 'h', h)
update_figs(handles)
set(handles.pushbutton_sel, 'Enable','on')


function menu_file_Callback(hObject, eventdata, handles)


function menu_open_Callback(hObject, eventdata, handles)
h = getappdata(handles.fig_bh_tomo_phase_pick, 'h');

[h.no_mog,h.db_file,hh] = choisirMOG('UserData', h.db_file);
delete(hh)
if h.no_mog==0
	return
end
%set(handles.fig_bh_tomo_phase_pick, 'Name', ['bh_tomo_phase_pick - ', h.db_file])

load(h.db_file,'mogs')
mog = mogs(h.no_mog);

h.Tx = [mog.data.Tx_x; mog.data.Tx_y; mog.data.Tx_z]';
h.uTx = unique(h.Tx, 'rows');
cmax=max(max(abs(mog.traces)));
cmin= -cmax;
h.cminmax = [cmin cmax];

setappdata(handles.fig_bh_tomo_phase_pick, 'mog', mog)
setappdata(handles.fig_bh_tomo_phase_pick, 'h', h)
set(handles.pushbutton_prec, 'Enable','off')
set(handles.pushbutton_sel, 'Enable','off')

function menu_save_Callback(hObject, eventdata, handles)
h = getappdata(handles.fig_bh_tomo_phase_pick, 'h');
pick = getappdata(handles.fig_bh_tomo_phase_pick, 'pick');
load(h.db_file,'mogs')
mogs(h.no_mog).tt = pick.tt;
mogs(h.no_mog).et = pick.et;
mogs(h.no_mog).tt_fait = pick.tt_fait; %#ok<NASGU>
save(h.db_file,'mogs','-append')


function menu_quit_Callback(hObject, eventdata, handles)
delete(handles.fig_bh_tomo_phase_pick)

function edit_medfilt_lth_Callback(hObject, eventdata, handles)


function edit_medfilt_lth_CreateFcn(hObject, eventdata, handles)

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


function pushbutton_suiv_Callback(hObject, eventdata, handles)
h = getappdata(handles.fig_bh_tomo_phase_pick, 'h');
if h.noTx == size(h.uTx,1)-1
	set(hObject,'Enable','off')
else
	set(hObject,'Enable','on')
end
set(handles.pushbutton_prec,'Enable','on')
h.noTx = h.noTx+1;
setappdata(handles.fig_bh_tomo_phase_pick, 'h', h)
update_figs(handles)


function pushbutton_prec_Callback(hObject, eventdata, handles)
h = getappdata(handles.fig_bh_tomo_phase_pick, 'h');
if h.noTx == 2
	set(hObject,'Enable','off')
else
	set(hObject,'Enable','on')
end
set(handles.pushbutton_suiv,'Enable','on')
h.noTx = h.noTx-1;
setappdata(handles.fig_bh_tomo_phase_pick, 'h', h)
update_figs(handles)


function pushbutton_sel_Callback(hObject, eventdata, handles)
axes(handles.traces_contigues);
xlim = get(handles.traces_contigues,'XLim');
ylim = get(handles.traces_contigues,'YLim');
[x,y,b] = ginput(1);
while x>xlim(1) && x<xlim(2) && y>ylim(1) && y<ylim(2) && b==1
	update_trace_simple(handles,round(x))
	[x,y,b] = ginput(1);
end



function edit_L_Callback(hObject, eventdata, handles)


function edit_L_CreateFcn(hObject, eventdata, handles)

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


function update_figs(handles)
update_trace_simple(handles)
update_traces_contigues(handles)


function update_traces_contigues(handles)
h = getappdata(handles.fig_bh_tomo_phase_pick, 'h');
mog = getappdata(handles.fig_bh_tomo_phase_pick, 'mog');
pick = getappdata(handles.fig_bh_tomo_phase_pick, 'pick');

ind = all(h.Tx==repmat(h.uTx(h.noTx,:),mog.data.ntrace,1),2)';

colormap ramac_cmap(16)
no_traces = 1:mog.data.ntrace;
imagesc(no_traces(ind), mog.data.timestp, mog.traces(:,ind),'Parent',handles.traces_contigues)
caxis(handles.traces_contigues, h.cminmax)
update_t_lim(handles)
hold(handles.traces_contigues, 'on')
ind = ind & (pick.tt ~= -1);
plot(handles.traces_contigues, no_traces(ind),pick.pAIC(ind),'o')
plot(handles.traces_contigues, no_traces(ind),pick.tt(ind),'go')
plot(handles.traces_contigues, no_traces(ind),pick.tt(ind)-pick.et(ind),'ro')
plot(handles.traces_contigues, no_traces(ind),pick.tt(ind)+pick.et(ind),'ro')
hold(handles.traces_contigues, 'off')


function update_trace_simple(handles,varargin)
mog = getappdata(handles.fig_bh_tomo_phase_pick, 'mog');
if nargin == 1
	h = getappdata(handles.fig_bh_tomo_phase_pick, 'h');
	ind = all(h.Tx==repmat(h.uTx(h.noTx,:),mog.data.ntrace,1),2)';
	no_traces = 1:mog.data.ntrace;
	no_traces = no_traces(ind);
	no = no_traces(1);
else
	no = varargin{1};
end
pick = getappdata(handles.fig_bh_tomo_phase_pick, 'pick');
plot(handles.trace_simple, mog.data.timestp, mog.traces(:,no),'Color',[0.75 0.75 0.75])
hold(handles.trace_simple,'on')
ind = find(pick.no_traces==no);
if ~isempty(ind)
	m = size(pick.datad,1);
	plot(handles.trace_simple, mog.data.timestp(1:m),pick.datad(:,ind))
end
if pick.tt_fait(no) ~= -1
	alim = [str2double(get(handles.edit_Amin,'String')) str2double(get(handles.edit_Amax,'String')) ];
	plot(handles.trace_simple, [pick.pAIC(no) pick.pAIC(no)],alim,'b')
	plot(handles.trace_simple, [pick.tt(no) pick.tt(no)],alim,'g')
	plot(handles.trace_simple, [pick.tt(no)-pick.et(no) pick.tt(no)-pick.et(no)],alim,'r')
	plot(handles.trace_simple, [pick.tt(no)+pick.et(no) pick.tt(no)+pick.et(no)],alim,'r')
end
hold(handles.trace_simple,'off')
update_t_lim(handles)
update_a_lim(handles)
title(handles.trace_simple, ['Fa = ',num2str(pick.Fa(no)*1e-6),' SNR = ',num2str(pick.SNRR(no))])


function edit_tmin_Callback(hObject, eventdata, handles)
update_t_lim(handles)

function edit_tmin_CreateFcn(hObject, eventdata, handles)

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


function edit_tmax_Callback(hObject, eventdata, handles)
update_t_lim(handles)

function edit_tmax_CreateFcn(hObject, eventdata, handles)

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


function edit_Amin_Callback(hObject, eventdata, handles)
update_a_lim(handles)

function edit_Amin_CreateFcn(hObject, eventdata, handles)

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


function edit_Amax_Callback(hObject, eventdata, handles)
update_a_lim(handles)

function edit_Amax_CreateFcn(hObject, eventdata, handles)

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


%
% -------------------------------------------------------------------------
%
function update_t_lim(handles)
tlim = [str2double(get(handles.edit_tmin,'String')) str2double(get(handles.edit_tmax,'String')) ];
set(handles.trace_simple,'XLim',tlim)
set(handles.traces_contigues, 'YLim', tlim)

%
% -------------------------------------------------------------------------
%
function update_a_lim(handles)
alim = [str2double(get(handles.edit_Amin,'String')) str2double(get(handles.edit_Amax,'String')) ];
set(handles.trace_simple,'YLim',alim)


function checkbox_inv_polarite_Callback(hObject, eventdata, handles)


function pushbutton_show_cwt_Callback(hObject, eventdata, handles)
