function varargout = bh_tomo_pick(varargin)
% BH_TOMO_PICK M-file for bh_tomo_pick.fig
%      BH_TOMO_PICK, by itself, creates a new BH_TOMO_PICK or raises the existing
%      singleton*.
%
%      H = BH_TOMO_PICK returns the handle to a new BH_TOMO_PICK or the handle to
%      the existing singleton*.
%
%      BH_TOMO_PICK('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in BH_TOMO_PICK.M with the given input arguments.
%
%      BH_TOMO_PICK('Property','Value',...) creates a new BH_TOMO_PICK or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before bh_tomo_pick_OpeningFunction gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to bh_tomo_pick_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Copyright (C) 2008 Bernard Giroux
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
% Edit the above text to modify the response to help bh_tomo_pick

% Last Modified by GUIDE v2.5 05-Dec-2012 14:14:15

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @bh_tomo_pick_OpeningFcn, ...
                   'gui_OutputFcn',  @bh_tomo_pick_OutputFcn, ...
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


% --- Executes just before bh_tomo_pick is made visible.
function bh_tomo_pick_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to bh_tomo_pick (see VARARGIN)

% Choose default command line output for bh_tomo_pick
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

h.db_file = get(handles.fig_bh_tomo_pick,'UserData');
h.saved = true;
h.bin_no = 0;
h.theta_int = 5;
h.scale_fdom = 0;
h.fdom = [];
h.fmoy = 0;
h.useDetrend = get(handles.UseDetrendMenuItem,'Checked');
set(handles.text_bin_no,'String',num2str(h.bin_no))
set(handles.edit_theta_int,'String',num2str(h.theta_int))
set(handles.checkbox_scale_fdom, 'Value',h.scale_fdom)

h.str = get_str_locale();
setappdata(handles.fig_bh_tomo_pick, 'h', h)
set_String_locale(handles, h.str)



% UIWAIT makes bh_tomo_pick wait for user response (see UIRESUME)
% uiwait(handles.fig_bh_tomo_pick);



% --- Outputs from this function are returned to the command line.
function varargout = bh_tomo_pick_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


function FileMenu_Callback(hObject, eventdata, handles)


function OpenMenuItem_Callback(hObject, eventdata, handles)
h = getappdata(handles.fig_bh_tomo_pick, 'h');

[h.no_mog,h.db_file,hh] = choisirMOG('UserData', h.db_file);
delete(hh)
if h.no_mog==0
	return
end
set(handles.fig_bh_tomo_pick, 'Name', ['bh_tomo_pick - ', h.db_file])

classer = false;
load(h.db_file,'mogs','auto_pick')
mog = mogs(h.no_mog);
h_set = false;
if exist('auto_pick','var')
	if ~isempty(auto_pick) && length(auto_pick) >= h.no_mog
		tmp = auto_pick(h.no_mog);
		if ~isempty(tmp.no_mog)
            db_file = h.db_file;
			h = tmp;
			h.db_file = db_file;
            h_set = true;
			nbins = length(h.bins);
			ntraces = 0;
			for n=1:nbins
				ntraces = ntraces+length(h.bins{n}.no_trace);
			end
			set(handles.text_data_pc,'String',...
				[num2str(100*ntraces/mog.data.ntrace,'%3.1f'),' %'])
			set(handles.checkbox_scale_fdom, 'Value',h.scale_fdom)
			h.str = get_str_locale();
		end
	end
end
clear mogs auto_pick

if ~h_set

    h.useDetrend = get(handles.UseDetrendMenuItem,'Checked');
    if strcmp(h.useDetrend,'on')==1
        h.proc_data = detrend_rad(mog.traces);
    else
        h.proc_data = mog.traces;
    end
    
	absmax = repmat(max(abs(h.proc_data)), size(h.proc_data,1),1);
	h.proc_data = h.proc_data ./ absmax;
	clear absmax
	h.proc_data2 = h.proc_data;
    h.data_1cycle = [];
    h.data_1cycle2 = [];
		
	h.z_min_Tx = min([mog.data.Rx_z mog.data.Tx_z]);
	h.z_max_Tx = max([mog.data.Rx_z mog.data.Tx_z]);
    h.z_min_Rx = h.z_min_Tx;
    h.z_max_Rx = h.z_max_Tx;
    
	classer = true;
end

if ~h_set || isempty(h.snr)
    
    win_snr = str2double(get(handles.edit_window_snr,'String'));
    win_snr = round(win_snr/mog.data.timec);

    fac_f = 1.0;
    fac_t = 1.0;
    if strcmp(mog.data.tunits,'ns')
        % radar: freq nominale donnée en MHz
        fac_f = 10^6;
        fac_t = 10^-9;
    elseif strcmp(mog.data.tunits,'ms')
        % sismique: on assume que f dominante est en kHz
        disp('warning : assuming Tx nominal frequency in _kHz_')
        fac_f = 10^3;
        fac_t = 10^-3;
    else
        disp('warning : assuming Tx nominal frequency in _Hz_')
        disp('warning : assuming time step in _s_')
    end

    f0_s = mog.data.rnomfreq * fac_f;
    dt_s = mog.data.timec * mog.fac_dt * fac_t;

    h.snr = data_select(h.proc_data, f0_s, dt_s, win_snr);

end


set(handles.edit_t_max,'String',num2str(round(mog.data.timestp(length(mog.data.timestp)))))
set(handles.edit_t_min,'String',num2str(round(mog.data.timestp(1))))
set(handles.edit_A_max,'String','1')
set(handles.edit_A_min,'String','-1')
set(handles.text_win_units,'String',['[',mog.data.tunits,']'])
set(handles.edit_z_min_Tx,'String',num2str(h.z_min_Tx))
set(handles.edit_z_max_Tx,'String',num2str(h.z_max_Tx))
set(handles.edit_z_min_Rx,'String',num2str(h.z_min_Rx))
set(handles.edit_z_max_Rx,'String',num2str(h.z_max_Rx))
set(handles.UseDetrendMenuItem,'Checked', h.useDetrend)

setappdata(handles.fig_bh_tomo_pick, 'h', h)
setappdata(handles.fig_bh_tomo_pick, 'mog', mog)

if classer
	classer_traces(handles)

	if get(handles.checkbox_1cycle,'Value')==1 || get(handles.checkbox_scale_fdom,'Value')==1
		prepare_1cycle(handles);
	end
end
update_figs(handles)



%
% -------------------------------------------------------------------------
%
function SaveMenuItem_Callback(hObject, eventdata, handles)
sauve_mat(handles, false)



%
% -------------------------------------------------------------------------
%
function QuitMenuItem_Callback(hObject, eventdata, handles)
h=getappdata(handles.fig_bh_tomo_pick, 'h');
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
delete(handles.fig_bh_tomo_pick)


%
% -------------------------------------------------------------------------
%
function set_String_locale(handles, str)

%set(handles.,          'String', str.s);


%
% -------------------------------------------------------------------------
%
function slider_amp_Callback(hObject, eventdata, handles)
val = get(hObject,'Value');
set(handles.edit_A_min,'String',num2str(-val))
set(handles.edit_A_max,'String',num2str(val))
set(handles.axes_mean_tr, 'YLim', [-val val])



%
% -------------------------------------------------------------------------
%
function slider_amp_CreateFcn(hObject, eventdata, handles)
% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end


%
% -------------------------------------------------------------------------
%
function pushbutton_precedente_Callback(hObject, eventdata, handles)
precedent(handles)

function precedent(handles)
h = getappdata(handles.fig_bh_tomo_pick, 'h');
h.bin_no = h.bin_no-1;
if (h.bin_no==1)
    set(handles.pushbutton_precedente,'Enable','off')
end
set(handles.pushbutton_suivate,'Enable','on')
set(handles.text_bin_no,'String',num2str(h.bins{h.bin_no}.theta_mid))
set(handles.text_no_it,'String',num2str(h.bins{h.bin_no}.no_it_align))
setappdata(handles.fig_bh_tomo_pick,'h',h)
update_bin_tr(handles)
update_mean_tr(handles,-1,-1)




%
% -------------------------------------------------------------------------
%
function pushbutton_suivate_Callback(hObject, eventdata, handles)
suivant(handles)

function suivant(handles)
h = getappdata(handles.fig_bh_tomo_pick, 'h');
h.bin_no = h.bin_no+1;
if (h.bin_no==numel(h.bins))
    set(handles.pushbutton_suivate,'Enable','off')
end

set(handles.pushbutton_precedente,'Enable','on')
set(handles.text_bin_no,'String',num2str(h.bins{h.bin_no}.theta_mid))
set(handles.text_no_it,'String',num2str(h.bins{h.bin_no}.no_it_align))
setappdata(handles.fig_bh_tomo_pick,'h',h)
update_bin_tr(handles)
update_mean_tr(handles,-1,-1)


%
% -------------------------------------------------------------------------
%
function pushbutton_aligner_Callback(hObject, eventdata, handles)
aligner(handles)

function aligner(handles)
h = getappdata(handles.fig_bh_tomo_pick, 'h');
work = (h.proc_data(:, h.bins{h.bin_no}.no_trace))';
if get(handles.checkbox_1cycle,'Value')==1 || get(handles.checkbox_scale_fdom,'Value')==1
    pds = str2double(get(handles.edit_poids1cycle,'String'));
    tr2 = (h.data_1cycle(:, h.bins{h.bin_no}.no_trace))';
    work = (1-pds)*work + pds*tr2;
    h.bins{h.bin_no}.work1cycle = true;
else
    h.bins{h.bin_no}.work1cycle = false;
end
trace = get_ref_trace(handles, h.bin_no);
ne = length(trace);
shift_lag = zeros(1,size(work,1));

for nt=1:size(work,1)
%    [c,lags] = xcorr(trace, traces(nt,:));
    c = xcorr(trace, work(nt,:));
    [xcm, ixcm] = max(c);
    shift_lag(nt) = ixcm-ne;
%    work(nt,:) = circshift(work(nt,:),[0 shift_lag(nt)]);
end

%figure(1)
%hist(shift_lag,20)

%moy = mean(shift_lag)
%et = std(shift_lag)


if isempty(h.data_1cycle)
    for nt=1:size(work,1)
        work(nt,:) = circshift(work(nt,:),[0 shift_lag(nt)]);
    end
    h.proc_data(:, h.bins{h.bin_no}.no_trace) = work';
else
    tr2 = (h.data_1cycle(:, h.bins{h.bin_no}.no_trace))';
    for nt=1:size(tr2,1)
        tr2(nt,:) = circshift(tr2(nt,:),[0 shift_lag(nt)]);
    end
    h.data_1cycle(:, h.bins{h.bin_no}.no_trace) = tr2';

    work = (h.proc_data(:, h.bins{h.bin_no}.no_trace))';
    for nt=1:size(work,1)
        work(nt,:) = circshift(work(nt,:),[0 shift_lag(nt)]);
    end
    h.proc_data(:, h.bins{h.bin_no}.no_trace) = work';
end

h.bins{h.bin_no}.no_it_align = h.bins{h.bin_no}.no_it_align+1;
set(handles.text_no_it,'String',num2str(h.bins{h.bin_no}.no_it_align))
h.saved = false;
setappdata(handles.fig_bh_tomo_pick, 'h', h)

update_bin_tr(handles)
update_mean_tr(handles,-1,-1)


%
% -------------------------------------------------------------------------
%
function pushbutton_pick_mean_tr_Callback(hObject, eventdata, handles)
pick_mean_tr(handles)

function pick_mean_tr(handles)
h = getappdata(handles.fig_bh_tomo_pick, 'h');

xlim = get(handles.axes_mean_tr,'XLim');
ylim = get(handles.axes_mean_tr,'YLim');

axes(handles.axes_mean_tr)
[t1,a,button] = ginput(1);
pushOnce = true;
tt = -1;
et = -1;
while button ~= 2
	if gca ~= handles.axes_mean_tr || t1<xlim(1) || t1>xlim(2) || a<ylim(1) || a>ylim(2)
		return
	end
	if button==1
		tt = t1;
		update_mean_tr(handles, tt, et)
	elseif button==3
		et = abs(t1-tt);
		update_mean_tr(handles, tt, et)
	end
	[t1,a,button] = ginput(1);
	pushOnce=false;
end
% button == 2
if ~pushOnce
	h.fait(h.bin_no) = true;
	h.pick(h.bin_no) = tt;
	h.pick_et(h.bin_no) = et;

	h.saved = false;
	setappdata(handles.fig_bh_tomo_pick, 'h', h)
	update_done(handles)
	if h.bin_no<length(h.fait)
		pushbutton_suivate_Callback([], [], handles)
	end
end


%
% -------------------------------------------------------------------------
%
function edit_t_min_Callback(hObject, eventdata, handles)
tmin = str2double(get(hObject,'String'));
tmax = str2double(get(handles.edit_t_max,'String'));
set(handles.axes_mean_tr,'XLim', [tmin tmax])
set(handles.axes_bin_tr,'XLim', [tmin tmax])


%
% -------------------------------------------------------------------------
%
function edit_t_min_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


%
% -------------------------------------------------------------------------
%
function edit_t_max_Callback(hObject, eventdata, handles)
tmin = str2double(get(handles.edit_t_min,'String'));
tmax = str2double(get(hObject,'String'));
set(handles.axes_mean_tr,'XLim', [tmin tmax])
set(handles.axes_bin_tr,'XLim', [tmin tmax])


%
% -------------------------------------------------------------------------
%
function edit_t_max_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


%
% -------------------------------------------------------------------------
%
function edit_A_min_Callback(hObject, eventdata, handles)


%
% -------------------------------------------------------------------------
%
function edit_A_min_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


%
% -------------------------------------------------------------------------
%
function edit_A_max_Callback(hObject, eventdata, handles)


%
% -------------------------------------------------------------------------
%
function edit_A_max_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



%
% -------------------------------------------------------------------------
%
function classer_traces(handles)
h = getappdata(handles.fig_bh_tomo_pick, 'h');
h.bin_no = str2double(get(handles.text_bin_no,'String'));
h.theta_int = abs(str2double(get(handles.edit_theta_int,'String')));

h.proc_data = h.proc_data2;

% on assume que les forages sont quasi parallÃ¨les, i.e. que l'angle
% d'Ã©mission Ã  Tx est Ã  peu prÃ¨s Ã©gal Ã  l'angle d'incidence Ã  Rx

% on calcule l'angle du rai (droit) p/r Ã  l'horizontale

% h = getappdata(handles.fig_bh_tomo_pick, 'h');   %%%YH
mog = getappdata(handles.fig_bh_tomo_pick, 'mog');

lrai = sqrt((mog.data.Tx_x-mog.data.Rx_x).^2+...
	(mog.data.Tx_y-mog.data.Rx_y).^2+...
	(mog.data.Tx_z-mog.data.Rx_z).^2);
dz = mog.data.Rx_z-mog.data.Tx_z;
theta_rais = 180*asin(dz./lrai)/pi;

no_trace = 1:length(mog.tt);

theta_start = -h.theta_int/2;
tmp = theta_start;
while 1
    tmp = tmp - h.theta_int;
    if ( tmp > -90 )
        theta_start = tmp;
    else
        break;
    end
end
theta_end = h.theta_int/2;
tmp = theta_end;
while 1
    tmp = tmp + h.theta_int;
    if (tmp < 90)
        theta_end = tmp;
    else
        break;
    end
end

nbins = round( (theta_end-theta_start)/h.theta_int );

lower = theta_start;
no_bin = 0;
h.bins = {};
ntraces = 0;
for n=1:nbins
    upper = lower+h.theta_int;
    
    ind = theta_rais >= lower & theta_rais < upper;
	ind2 = mog.data.Rx_z <= h.z_max_Rx & mog.data.Rx_z >= h.z_min_Rx;
	ind3 = mog.data.Tx_z <= h.z_max_Tx & mog.data.Tx_z >= h.z_min_Tx;
	
	ind = ind & ind2 & ind3 & mog.in;
	
    if (sum(ind)>0)
        no_bin = no_bin+1;
        bin.no_trace = no_trace(ind);
        bin.theta_mid = 0.5*(lower+upper);
        bin.no_it_align = 0;
        h.bins{no_bin} = bin;
        h.bins{no_bin}.work1cycle = false;
        ntraces = ntraces + sum(ind);
    end    
    lower = upper;
end
h.fait = false(1,numel(h.bins));
h.pick = zeros(1,numel(h.bins))-1;
h.pick_et = zeros(1,numel(h.bins))-1;

h.bin_no = 1;
set(handles.text_bin_no,'String',num2str(h.bins{h.bin_no}.theta_mid))
set(handles.pushbutton_precedente,'Enable','off')
h.saved = false;
setappdata(handles.fig_bh_tomo_pick,'h',h)
set(handles.text_data_pc,'String',...
    [num2str(100*ntraces/sum(mog.in),'%3.1f'),' %'])





%
% -------------------------------------------------------------------------
%
function pushbutton_xcorr_Callback(hObject, eventdata, handles)
do_xcorr(handles)

function do_xcorr(handles)
h = getappdata(handles.fig_bh_tomo_pick, 'h');
mog = getappdata(handles.fig_bh_tomo_pick, 'mog');
seuil_fdom = str2double(get(handles.edit_seuil_fdom,'string'));


if isfield(h, 'snr')==false || isempty(h.snr)
	
	win_snr = str2double(get(handles.edit_window_snr,'String'));
	win_snr = round(win_snr/mog.data.timec);

	fac_f = 1.0;
	fac_t = 1.0;
	if strcmp(mog.data.tunits,'ns')
		% radar: freq nominale donnée en MHz
		fac_f = 10^6;
		fac_t = 10^-9;
	elseif strcmp(mog.data.tunits,'ms')
		% sismique: on assume que f dominante est en kHz
		disp('warning : assuming Tx nominal frequency in _kHz_')
		fac_f = 10^3;
		fac_t = 10^-3;
	else
		disp('warning : assuming Tx nominal frequency in _Hz_')
		disp('warning : assuming time step in _s_')
	end

	f0_s = mog.data.rnomfreq * fac_f;
	dt_s = mog.data.timec * mog.fac_dt * fac_t;

	h.snr = data_select(h.proc_data, f0_s, dt_s, win_snr);
end
if isfield(h, 'fdom')==false
	h.fdom = [];
end
if isempty(h.fdom)
	fmin = 1e-3*str2double(get(handles.edit_f_min,'String'));
	fmax = 1e-3*str2double(get(handles.edit_f_max,'String'));
	dt = mog.data.timec;
	f0 = 1e-3*mog.data.rnomfreq;

	Fs = 1/dt;
	HalfFs = Fs/2;
	Wp = 1.4*f0/HalfFs; Ws = 1.6*f0/HalfFs;
	Rp = 3; Rs = 40;
	[n,Wn] = cheb1ord(Wp,Ws,Rp,Rs);

	[b,a] = cheby1(n,0.5,Wn);

	if strcmp(h.useDetrend,'on')==1
		traces = detrend_rad(mog.traces);
	else
		traces = mog.traces;
	end
	absmax = repmat(max(abs(traces)), size(traces,1),1);
	traces = traces ./ absmax;

	nt = size(traces,2);
	h.fdom = zeros(1,nt)-1;
    hwb = waitbar(0,'Calcul des fréquences dominantes');
	inc_wb = 100;
    for n=1:nt
        trace = traces(:,n);
        if any(isfinite(trace))
            h.fdom(n) = get_fdom(filter(b,a,trace), dt, f0, fmin, fmax);
        else
            h.fdom(n) = -1;
        end
        if rem(n, inc_wb)==0, waitbar(n/nt,hwb), end
    end
	close(hwb)
end

show = get(handles.checkbox_show_fit,'Value');
if ~show, hwb = waitbar(0,h.str.s276); end
if show
    hshow = figure();
end
for nb = 1:numel(h.bins)
    if h.fait(nb) == false
        continue
    end
    traces = (h.proc_data2(:, h.bins{nb}.no_trace))';
	fdom = h.fdom(h.bins{nb}.no_trace);
	snr = h.snr(h.bins{nb}.no_trace);
	To4 = 1e3 * 0.25./fdom;  % période / 4
	To4(fdom==-1) = 0.5/(1e-3*mog.data.rnomfreq);
    
    if h.bins{nb}.work1cycle
        pds = str2double(get(handles.edit_poids1cycle,'String'));
        tr2 = (h.data_1cycle2(:, h.bins{nb}.no_trace))';
        traces = (1-pds)*traces + pds*tr2;
    end
    
    trace = get_ref_trace(handles, nb, h.bins{nb}.work1cycle);
    
%    itt = findnear(h.pick(nb), mog.data.timestp);
%    iet = findnear(h.pick(nb)-h.pick_et(nb), mog.data.timestp);
    
    ne = length(trace);
    shift_lag = zeros(1,size(traces,1));
    corr = zeros(1,size(traces,1));
	pick_et = h.pick_et(nb);
	b = pick_et ./ To4;
	m = (b-1)/-3;
    for nt=1:size(traces,1)
        %    [c,lags] = xcorr(trace, traces(nt,:));
        c = xcorr(trace, traces(nt,:),'coeff');
        [~, ixcm] = max(c);
        shift_lag(nt) = ixcm-ne;
		lsnr = log10(snr(nt));
		if lsnr>=3
			corr(nt) = 1;
		else
			corr(nt) = m(nt)*lsnr + b(nt);
		end
%        corr(nt) = xcm;
    end
    pick = h.pick(nb) - shift_lag*mog.data.timec;
    pick_et = pick_et./corr;
    mog.tt(h.bins{nb}.no_trace) = pick;
    mog.et(h.bins{nb}.no_trace) = pick_et;
    if h.scale_fdom==1
        % on remet à l'échelle initiale
        fac = h.fdom(h.bins{nb}.no_trace);
        ind1 = fac==-1;
        ind2 = h.snr(h.bins{nb}.no_trace) <= seuil_fdom;
        fac(ind1 | ind2') = h.fmoy;
        fac = h.fmoy./fac;
        mog.tt(h.bins{nb}.no_trace) = mog.tt(h.bins{nb}.no_trace) .* fac;
        mog.et(h.bins{nb}.no_trace) = mog.et(h.bins{nb}.no_trace) .* fac;
    end
    mog.tt_done(h.bins{nb}.no_trace) = true;
    
    if show
        figure(hshow)
        yt = 1:size(traces,1);
        imagesc(mog.data.timestp, yt, traces)
        colormap(gray)
        hold on
        plot(pick,yt,'g')
        plot(pick-pick_et,yt,'r--')
        plot(pick+pick_et,yt,'r--')
        hold off
        title(num2str(h.bins{nb}.theta_mid))
        pause
    else
        waitbar(nb/numel(h.bins), hwb)
    end
end
if ~show, close(hwb), end
h.saved = false;
setappdata(handles.fig_bh_tomo_pick, 'mog', mog)
setappdata(handles.fig_bh_tomo_pick, 'h', h)




%
% -------------------------------------------------------------------------
%
function sauve_mat(handles, msg)
if msg
	hh=msgbox('Sauvegarde intermediaire');
else
	hh=msgbox('Sauvegarde en cours');
end
h = getappdata(handles.fig_bh_tomo_pick, 'h');
mog = getappdata(handles.fig_bh_tomo_pick, 'mog');

load(h.db_file,'mogs','air','auto_pick')
mogs(h.no_mog) = mog; %#ok<NASGU>

if isempty( auto_pick ) %#ok<NODEF>
	clear auto_pick
end

auto_pick(h.no_mog) = h; %#ok<NASGU>

save(h.db_file,'mogs','air','auto_pick','-append')
if ishandle(hh)==1
	delete(hh)
end
h.saved = true;
setappdata(handles.fig_bh_tomo_pick, 'h', h)




%
% -------------------------------------------------------------------------
%
function update_figs(handles)
h = getappdata(handles.fig_bh_tomo_pick, 'h');
set(handles.text_bin_no,'String',num2str(h.bins{h.bin_no}.theta_mid))
set(handles.text_no_it,'String',num2str(h.bins{h.bin_no}.no_it_align))
update_done(handles)
update_bin_tr(handles)
update_mean_tr(handles,-1,-1)


%
% -------------------------------------------------------------------------
%
function update_done(handles)
h = getappdata(handles.fig_bh_tomo_pick, 'h');
imagesc(h.fait, 'Parent',handles.axes_done)
set(handles.axes_done,'CLim',[0 1])
set(handles.axes_done,'YTick',[])
ind = 5:5:numel(h.bins);
theta = zeros(numel(h.bins),1);
for n=1:length(theta)
    theta(n) = h.bins{n}.theta_mid;
end
set(handles.axes_done,'XTick',ind);
set(handles.axes_done,'XTickLabel',num2str(theta(ind)));
%colormap(handles.axes_done,[1 0 0; 1 0 0; 0 1 0])

%if (sum(h.fait)==numel(h.fait))
%    set(handles.pushbutton_xcorr,'Enable','on')
%end


%
% -------------------------------------------------------------------------
%
function update_bin_tr(handles)
h = getappdata(handles.fig_bh_tomo_pick, 'h');
mog = getappdata(handles.fig_bh_tomo_pick, 'mog');

traces = (h.proc_data(:, h.bins{h.bin_no}.no_trace))';
if ( get(handles.checkbox_1cycle,'Value')==1 || get(handles.checkbox_scale_fdom,'Value')==1 ) ...
        && get(handles.checkbox_showOrig,'Value')==0 && ~isempty(h.data_1cycle)
    pds = str2double(get(handles.edit_poids1cycle,'String'));
    tr2 = (h.data_1cycle(:, h.bins{h.bin_no}.no_trace))';
    traces = (1-pds)*traces + pds*tr2;
end

yt = 1:size(traces,1);
imagesc(mog.data.timestp, yt, traces, 'Parent', handles.axes_bin_tr)
colormap(handles.axes_bin_tr,gray)
tmin = str2double(get(handles.edit_t_min, 'String'));
tmax = str2double(get(handles.edit_t_max, 'String'));
set(handles.axes_bin_tr,'XLim', [tmin tmax])



%
% -------------------------------------------------------------------------
%
function update_mean_tr(handles, tt, et)
mog = getappdata(handles.fig_bh_tomo_pick, 'mog');
h = getappdata(handles.fig_bh_tomo_pick, 'h');

if ( h.bins{h.bin_no}.no_it_align == 0 )
    traces = (h.proc_data(:, h.bins{h.bin_no}.no_trace))';
    %snr = zeros(1,size(traces,1));
    %win_snr = str2double(get(handles.edit_window_snr,'String'));
    %win_snr = round(win_snr/mog.data.timec);
    %for n=1:length(snr)
    %    snr(n) = 1/rmsv(traces(n,1:win_snr));
    %end
    [m, i] = max(h.snr(h.bins{h.bin_no}.no_trace));
    if ( get(handles.checkbox_1cycle,'Value')==1 || ...
            get(handles.checkbox_scale_fdom,'Value')==1 ) && ...
        ~isempty(h.data_1cycle) 
        traces = (h.data_1cycle(:, h.bins{h.bin_no}.no_trace))';
    end
    trace = traces(i,:);
else
    traces = (h.proc_data(:, h.bins{h.bin_no}.no_trace))';
    if ( get(handles.checkbox_1cycle,'Value')==1 || get(handles.checkbox_scale_fdom,'Value')==1 )...
            && get(handles.checkbox_showOrig,'Value')==0
        pds = str2double(get(handles.edit_poids1cycle,'String'));
        tr2 = (h.data_1cycle(:, h.bins{h.bin_no}.no_trace))';
        traces = (1-pds)*traces + pds*tr2;
    end
    trace = mean(traces);
end

plot(handles.axes_mean_tr, mog.data.timestp, trace)
if tt==-1
	h = getappdata(handles.fig_bh_tomo_pick, 'h');
	tt = h.pick(h.bin_no);
	et = h.pick_et(h.bin_no);
end
if tt>-1
	hold on
	plot(handles.axes_mean_tr, [tt tt], [-1 1],'g')
	if et>-1
		plot(handles.axes_mean_tr, [tt-et tt-et], [-1 1],'r')
		plot(handles.axes_mean_tr, [tt+et tt+et], [-1 1],'r')
	end
	hold off
end
tmin = str2double(get(handles.edit_t_min, 'String'));
tmax = str2double(get(handles.edit_t_max, 'String'));
Amin = str2double(get(handles.edit_A_min, 'String'));
Amax = str2double(get(handles.edit_A_max, 'String'));
set(handles.axes_mean_tr, 'XLim',[tmin tmax],'YLim', [Amin Amax])




%
% -------------------------------------------------------------------------
%
function trace = get_ref_trace(handles, bin_no, varargin)
h = getappdata(handles.fig_bh_tomo_pick, 'h');
%mog = getappdata(handles.fig_bh_tomo_pick, 'mog');
traces = (h.proc_data(:, h.bins{bin_no}.no_trace))';
%win_snr = str2double(get(handles.edit_window_snr,'String'));
%win_snr = round(win_snr/mog.data.timec);

if nargin==3
    work1cycle = varargin{1};
else
    work1cycle = get(handles.checkbox_1cycle,'Value')==1 || get(handles.checkbox_scale_fdom,'Value')==1;
end

if ( h.bins{bin_no}.no_it_align == 0 )
%    snr = zeros(1,size(traces,1));
%    for n=1:length(snr)
%        snr(n) = 1/rmsv(traces(n,1:win_snr));
%    end
    [m, i] = max(h.snr(h.bins{bin_no}.no_trace));
    if work1cycle
        pds = str2double(get(handles.edit_poids1cycle,'String'));
        tr2 = (h.data_1cycle(:, h.bins{bin_no}.no_trace))';
        traces = (1-pds)*traces + pds*tr2;
    end
    trace = traces(i,:);
else
    if work1cycle
        pds = str2double(get(handles.edit_poids1cycle,'String'));
        tr2 = (h.data_1cycle(:, h.bins{bin_no}.no_trace))';
        traces = (1-pds)*traces + pds*tr2;
    end
    trace = mean(traces);
end



%
% -------------------------------------------------------------------------
%
function edit_window_snr_Callback(hObject, eventdata, handles)
h = getappdata(handles.fig_bh_tomo_pick, 'h');
mog = getappdata(handles.fig_bh_tomo_pick, 'mog');

win_snr = str2double(get(handles.edit_window_snr,'String'));
win_snr = round(win_snr/mog.data.timec);

fac_f = 1.0;
fac_t = 1.0;
if strcmp(mog.data.tunits,'ns')
    % radar: freq nominale donnée en MHz
    fac_f = 10^6;
    fac_t = 10^-9;
elseif strcmp(mog.data.tunits,'ms')
    % sismique: on assume que f dominante est en kHz
    disp('warning : assuming Tx nominal frequency in _kHz_')
    fac_f = 10^3;
    fac_t = 10^-3;
else
    disp('warning : assuming Tx nominal frequency in _Hz_')
    disp('warning : assuming time step in _s_')
end

f0_s = mog.data.rnomfreq * fac_f;
dt_s = mog.data.timec * mog.fac_dt * fac_t;

h.snr = data_select(h.proc_data, f0_s, dt_s, win_snr);
h.saved = false;
setappdata(handles.fig_bh_tomo_pick, 'h', h)


%
% -------------------------------------------------------------------------
%
function edit_window_snr_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


%
% -------------------------------------------------------------------------
%
function edit_z_min_Rx_Callback(hObject, eventdata, handles)
h = getappdata(handles.fig_bh_tomo_pick, 'h');
h.z_min_Rx = str2double(get(hObject,'String'));
setappdata(handles.fig_bh_tomo_pick, 'h', h)
classer_traces(handles)
update_figs(handles)


%
% -------------------------------------------------------------------------
%
function edit_z_min_Rx_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


%
% -------------------------------------------------------------------------
%
function edit_z_max_Rx_Callback(hObject, eventdata, handles)
h = getappdata(handles.fig_bh_tomo_pick, 'h');
h.z_max_Rx = str2double(get(hObject,'String'));
setappdata(handles.fig_bh_tomo_pick, 'h', h)
classer_traces(handles)
update_figs(handles)


%
% -------------------------------------------------------------------------
%
function edit_z_max_Rx_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


%
% -------------------------------------------------------------------------
%
function edit_z_min_Tx_Callback(hObject, eventdata, handles)
h = getappdata(handles.fig_bh_tomo_pick, 'h');
h.z_min_Tx = str2double(get(hObject,'String'));
setappdata(handles.fig_bh_tomo_pick, 'h', h)
classer_traces(handles)
update_figs(handles)



%
% -------------------------------------------------------------------------
%
function edit_z_min_Tx_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



%
% -------------------------------------------------------------------------
%
function edit_z_max_Tx_Callback(hObject, eventdata, handles)
h = getappdata(handles.fig_bh_tomo_pick, 'h');
h.z_max_Tx = str2double(get(hObject,'String'));
setappdata(handles.fig_bh_tomo_pick, 'h', h)
classer_traces(handles)
update_figs(handles)



%
% -------------------------------------------------------------------------
%
function edit_z_max_Tx_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function edit_theta_int_Callback(hObject, eventdata, handles)
classer_traces(handles)
update_figs(handles)

%--------------------------------------------------------------
function edit_theta_int_CreateFcn(hObject, eventdata, handles)

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



%
% -------------------------------------------------------------------------
%
function pushbutton_reinit_Callback(hObject, eventdata, handles)
h = getappdata(handles.fig_bh_tomo_pick, 'h');
h.fait(h.bin_no) = false;
h.pick(h.bin_no) = -1;
h.pick_et(h.bin_no) = -1;
h.bins{h.bin_no}.no_it_align = 0;
h.proc_data(:, h.bins{h.bin_no}.no_trace) = h.proc_data2(:, h.bins{h.bin_no}.no_trace);
if get(handles.checkbox_1cycle,'Value')==1 || get(handles.checkbox_scale_fdom,'Value')==1
    h.data_1cycle(:, h.bins{h.bin_no}.no_trace) = h.data_1cycle2(:, h.bins{h.bin_no}.no_trace);
end
h.saved = false;
setappdata(handles.fig_bh_tomo_pick, 'h', h)
update_figs(handles)



%
% -------------------------------------------------------------------------
%
function checkbox_show_fit_Callback(hObject, eventdata, handles)
update_figs(handles)



%
% -------------------------------------------------------------------------
%
function checkbox_1cycle_Callback(hObject, eventdata, handles)
if get(hObject,'Value')==1
    set(handles.checkbox_showOrig,'Enable','on');
else
    set(handles.checkbox_showOrig,'Value',0);
    set(handles.checkbox_showOrig,'Enable','off');
end
update_figs(handles)


%
% -------------------------------------------------------------------------
%
function checkbox_scale_fdom_Callback(hObject, eventdata, handles)
if get(hObject,'Value')==1
    set(handles.checkbox_showOrig,'Enable','on');
else
    set(handles.checkbox_showOrig,'Value',0);
    set(handles.checkbox_showOrig,'Enable','off');
end
update_figs(handles)


%
% -------------------------------------------------------------------------
%
function prepare_1cycle(handles)
h = getappdata(handles.fig_bh_tomo_pick, 'h');
mog = getappdata(handles.fig_bh_tomo_pick, 'mog');

fmin = 1e-3*str2double(get(handles.edit_f_min,'String'));
fmax = 1e-3*str2double(get(handles.edit_f_max,'String'));
seuil_xcorr = str2double(get(handles.edit_seuil_1cycle,'String'));
seuil = str2double(get(handles.edit_seuil,'string'));
seuil_fdom = str2double(get(handles.edit_seuil_fdom,'string'));


dt = mog.data.timec;
f0 = 1e-3*mog.data.rnomfreq;

Fs = 1/dt;
HalfFs = Fs/2;
Wp = 1.4*f0/HalfFs; Ws = 1.6*f0/HalfFs;
Rp = 3; Rs = 40;
[n,Wn] = cheb1ord(Wp,Ws,Rp,Rs);

[b,a] = cheby1(n,0.5,Wn);

if strcmp(h.useDetrend,'on')==1
    h.proc_data = detrend_rad(mog.traces);
else
    h.proc_data = mog.traces;
end
absmax = repmat(max(abs(h.proc_data)), size(h.proc_data,1),1);
h.proc_data = h.proc_data ./ absmax;
clear absmax
h.proc_data2 = h.proc_data;
traces = h.proc_data;


%debug=0;
nt = size(traces,2);

h.fdom = zeros(1,nt)-1;
test=0;
inc_wb = round(nt/100);

if get(handles.checkbox_1cycle,'Value')==1
    hwb = waitbar(0,'Isolation of the first cycle');
    for n=1:nt
        trace = traces(:,n);
%        snr(n) = 1/rmsv(trace(1:win_snr));
        if h.snr(n) > seuil
%			if h.snr(n) < 50
				[traces(:,n), h.fdom(n)] = get1cycle(trace, b, a, dt, f0, seuil_xcorr, fmin, fmax);
%			else
%				[traces(:,n), h.fdom(n)] = get1cycle(trace, [], [], dt, f0, seuil_xcorr, fmin, fmax);
%			end
            if test==1 && rem(n,100)==0
                figure(1)
                plot(mog.data.timestp, trace, 'b', ...
                    mog.data.timestp, filtfilt(b,a,trace), 'r', ...
                    mog.data.timestp, traces(:,n),'g', ...
                    mog.data.timestp([idebut idebut]),[-1 1],'k')
                title(['fdom = ',num2str(1000*fdom),', snr = ',num2str(h.snr(n))])
                pause
            end
        else
            h.fdom(n) = get_fdom(filter(b,a,trace), dt, f0, fmin, fmax);
        end
        if rem(n, inc_wb)==0, waitbar(n/nt,hwb), end
    end
    close(hwb)
else
    hwb = waitbar(0,'Calculation of dominant frequencies');
    for n=1:nt
        trace = traces(:,n);
        %snr(n) = 1/rmsv(trace(1:win_snr));
%		if h.snr(n) < 50
			h.fdom(n) = get_fdom(filter(b,a,trace), dt, f0, fmin, fmax);
%		else
%			h.fdom(n) = get_fdom(trace, dt, f0, fmin, fmax);
%		end
        if rem(n, inc_wb)==0, waitbar(n/nt,hwb), end
    end
    close(hwb)
end

if get(handles.checkbox_scale_fdom,'Value')==1
    h.scale_fdom = 1;
	hwb = waitbar(0,'Scaling traces');
	h.fmoy = mean(h.fdom(h.fdom~=-1));
	for n=1:nt
		if h.fdom(n) ~= -1 && h.snr(n) > seuil_fdom
            if test==1
                figure(10)
                plot(mog.data.timestp, traces(:,n),'g')
                hold on
            end
			timestp = mog.data.timestp * h.fdom(n)/h.fmoy;
			traces(:,n) = interp1(timestp, traces(:,n), mog.data.timestp, [], 0);
            if test==1
                plot(mog.data.timestp, traces(:,n),'b')
                title(['ratio = ',num2str(h.fdom(n)/h.fmoy),', fmoy = ',num2str(h.fmoy),', fdom = ',num2str(h.fdom(n))])
                hold off
                pause
            end
			h.proc_data(:,n) = interp1(timestp, h.proc_data(:,n), mog.data.timestp, [], 0);
		end
		if rem(n, inc_wb)==0, waitbar(n/nt,hwb), end
	end	
	close(hwb)
	h.proc_data2 = h.proc_data;
else
    h.scale_fdom = 0;
end
h.data_1cycle = traces;
h.data_1cycle2 = traces;
setappdata(handles.fig_bh_tomo_pick, 'h', h)



%
% -------------------------------------------------------------------------
%
function checkbox_showOrig_Callback(hObject, eventdata, handles)
update_figs(handles)



%
% -------------------------------------------------------------------------
%
function edit_seuil_Callback(hObject, eventdata, handles)
%prepare_1cycle(handles)
%update_figs(handles)


%
% -------------------------------------------------------------------------
%
function edit_seuil_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end




%
% -------------------------------------------------------------------------
%
function [indmin, indmax, indzer] = extr(x)
%EXTR  finds extrema and zero-crossings
%
% [indmin, indmax, indzer] = EXTR(x)
%
% inputs : - x : analyzed signal
%
% outputs : - indmin = indices of minima
%           - indmax = indices of maxima
%           - indzer = indices of zero-crossings
%
% Adapté de:
% G. Rilling, last modification: July 2002
% gabriel.rilling@ens-lyon.fr

m = length(x);

if nargout > 2
    x1=x(1:m-1);
    x2=x(2:m);
    indzer = find(x1.*x2<0);
    
    if any(x == 0)
        iz = find( x==0 );
        indz = []; %#ok<NASGU>
        if any(diff(iz)==1)
            zer = x == 0;
            dz = diff([0 zer 0]);
            debz = find(dz == 1);
            finz = find(dz == -1)-1;
            indz = round((debz+finz)/2);
        else
            indz = iz;
        end
        indzer = sort([indzer indz]);
    end
end
  
d = diff(x);

n = length(d);
d1 = d(1:n-1);
d2 = d(2:n);
indmin = find(d1.*d2<0 & d1<0)+1;
indmax = find(d1.*d2<0 & d1>0)+1;

if any(d==0)
  
  imax = [];
  imin = [];
  
  bad = (d==0);
  dd = diff([0 bad 0]);
  debs = find(dd == 1);
  fins = find(dd == -1);
  if debs(1) == 1
    if length(debs) > 1
      debs = debs(2:end);
      fins = fins(2:end);
    else
      debs = [];
      fins = [];
    end
  end
  if isempty(debs)==false
    if fins(end) == m
      if length(debs) > 1
        debs = debs(1:(end-1));
        fins = fins(1:(end-1));

      else
        debs = [];
        fins = [];
      end      
    end
  end
  lc = length(debs);
  if lc > 0
    for k = 1:lc
      if d(debs(k)-1) > 0
        if d(fins(k)) < 0
          imax = [imax round((fins(k)+debs(k))/2)]; %#ok<AGROW>
        end
      else
        if d(fins(k)) > 0
          imin = [imin round((fins(k)+debs(k))/2)]; %#ok<AGROW>
        end
      end
    end
  end
  
  if isempty(imax)==false
    indmax = sort([indmax imax]);
  end

  if isempty(imin)==false
    indmin = sort([indmin imin]);
  end
  
end  



%
% -------------------------------------------------------------------------
%
function edit_poids1cycle_Callback(hObject, eventdata, handles)
update_figs(handles)



%
% -------------------------------------------------------------------------
%
function edit_poids1cycle_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


%
% -------------------------------------------------------------------------
%
function pushbutton_show1trace_Callback(hObject, eventdata, handles)
h = getappdata(handles.fig_bh_tomo_pick, 'h');
mog = getappdata(handles.fig_bh_tomo_pick, 'mog');
xlim = get(handles.axes_bin_tr,'XLim');
ylim = get(handles.axes_bin_tr,'YLim');

axes(handles.axes_bin_tr)
[t1,a] = ginput(1);
if gca ~= handles.axes_bin_tr || t1<xlim(1) || t1>xlim(2) || a<ylim(1) || a>ylim(2)
    return
end

no = h.bins{h.bin_no}.no_trace(round(a));
%win_snr = str2double(get(handles.edit_window_snr,'String'));
%win_snr = round(win_snr/mog.data.timec);
%snr = 1/rmsv(h.proc_data(1:win_snr,no));

figure
plot(mog.data.timestp,h.proc_data(:,no))
if get(handles.checkbox_1cycle,'Value')==1 || get(handles.checkbox_scale_fdom,'Value')==1
    hold on
    pds = str2double(get(handles.edit_poids1cycle,'String'));
    trace = (1-pds)*h.proc_data(:,no) + pds*h.data_1cycle(:,no);
    plot(mog.data.timestp,trace,'g')
    hold off
end
xlabel('Time')
title(['Trace no ',num2str(mog.no_traces(no)),', Elevation - Tx: ',num2str(mog.data.Tx_z(no)),', Rx: ',...
    num2str(mog.data.Rx_z(no)),', SNR: ',num2str(h.snr(no))])


%
% -------------------------------------------------------------------------
%
function pushbutton_eliminer1trace_Callback(hObject, eventdata, handles)
h = getappdata(handles.fig_bh_tomo_pick, 'h');
mog = getappdata(handles.fig_bh_tomo_pick, 'mog');
xlim = get(handles.axes_bin_tr,'XLim');
ylim = get(handles.axes_bin_tr,'YLim');

axes(handles.axes_bin_tr)
[t1,a,b] = ginput(1);
while b==1
    if gca ~= handles.axes_bin_tr || t1<xlim(1) || t1>xlim(2) || a<ylim(1) || a>ylim(2)
        return
    end

    no = h.bins{h.bin_no}.no_trace(round(a));

    ind = h.bins{h.bin_no}.no_trace == no;
    h.bins{h.bin_no}.no_trace = h.bins{h.bin_no}.no_trace(~ind);

    setappdata(handles.fig_bh_tomo_pick, 'h', h)

    update_figs(handles)
    nbins = length(h.bins);
    ntraces = 0;
    for n=1:nbins
        ntraces = ntraces+length(h.bins{n}.no_trace);
    end
    set(handles.text_data_pc,'String',...
        [num2str(100*ntraces/mog.data.ntrace,'%3.1f'),' %'])

    axes(handles.axes_bin_tr)
    [t1,a,b] = ginput(1);
end



function edit_seuil_1cycle_Callback(hObject, eventdata, handles)


function edit_seuil_1cycle_CreateFcn(hObject, eventdata, handles)

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


function edit_f_min_Callback(hObject, eventdata, handles)


function edit_f_min_CreateFcn(hObject, eventdata, handles)

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


function edit_f_max_Callback(hObject, eventdata, handles)


function edit_f_max_CreateFcn(hObject, eventdata, handles)

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


function edit_seuil_fdom_Callback(hObject, eventdata, handles)


function edit_seuil_fdom_CreateFcn(hObject, eventdata, handles)

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


function pushbutton_prepare1cycle_Callback(hObject, eventdata, handles)
h = getappdata(handles.fig_bh_tomo_pick, 'h');
reponse = true;
if ~isempty(h.data_1cycle)
    b=questdlg('Previously processed data will be overwritten.  Proceed?');
    if strcmp(b,'Yes')==0
        reponse = false;
    end
end
if reponse==true
    prepare_1cycle(handles);
end


function EditMenu_Callback(hObject, eventdata, handles)


function ReinitMenuItem_Callback(hObject, eventdata, handles)
b=questdlg('About to reinitialize everything.  Proceed?');
if strcmp(b,'Yes')==0
    return
end
h = getappdata(handles.fig_bh_tomo_pick, 'h');

load(h.db_file,'mogs')
mog = mogs(h.no_mog);

if strcmp(h.useDetrend,'on')==1
    h.proc_data = detrend_rad(mog.traces);
else
    h.proc_data = mog.traces;
end
absmax = repmat(max(abs(h.proc_data)), size(h.proc_data,1),1);
h.proc_data = h.proc_data ./ absmax;
clear absmax
h.proc_data2 = h.proc_data;
h.data_1cycle = [];
h.data_1cycle2 = [];

h.z_min_Tx = min([mog.data.Rx_z mog.data.Tx_z]);
h.z_max_Tx = max([mog.data.Rx_z mog.data.Tx_z]);
h.z_min_Rx = h.z_min_Tx;
h.z_max_Rx = h.z_max_Tx;

set(handles.edit_t_max,'String',num2str(round(mog.data.timestp(length(mog.data.timestp)))))
set(handles.edit_t_min,'String',num2str(round(mog.data.timestp(1))))
set(handles.edit_A_max,'String','1')
set(handles.edit_A_min,'String','-1')
set(handles.text_win_units,'String',['[',mog.data.tunits,']'])
set(handles.edit_z_min_Tx,'String',num2str(h.z_min_Tx))
set(handles.edit_z_max_Tx,'String',num2str(h.z_max_Tx))
set(handles.edit_z_min_Rx,'String',num2str(h.z_min_Rx))
set(handles.edit_z_max_Rx,'String',num2str(h.z_max_Rx))
set(handles.checkbox_1cycle,'Value',0)
set(handles.checkbox_scale_fdom,'Value',0)

setappdata(handles.fig_bh_tomo_pick, 'h', h)
setappdata(handles.fig_bh_tomo_pick, 'mog', mog)

classer_traces(handles)

update_figs(handles)


function UseDetrendMenuItem_Callback(hObject, eventdata, handles)
b=questdlg('About to reinitialize everything.  Proceed?');
if strcmp(b,'Yes')==0
    return
end

h = getappdata(handles.fig_bh_tomo_pick, 'h');
mog = getappdata(handles.fig_bh_tomo_pick, 'mog');
if strcmp(h.useDetrend,'on')==1
    h.useDetrend = 'off';
    set(handles.UseDetrendMenuItem,'Checked','off');
else
    h.useDetrend = 'on';
    set(handles.UseDetrendMenuItem,'Checked','on');
end

if get(handles.checkbox_1cycle,'Value')==1 || get(handles.checkbox_scale_fdom,'Value')==1
    prepare_1cycle(handles)
else
    if strcmp(h.useDetrend,'on')==1
        h.proc_data = detrend_rad(mog.traces);
    else
        h.proc_data = mog.traces;
    end
    absmax = repmat(max(abs(h.proc_data)), size(h.proc_data,1),1);
    h.proc_data = h.proc_data ./ absmax;
    clear absmax
    h.proc_data2 = h.proc_data;
    h.data_1cycle = [];
    h.data_1cycle2 = [];
    setappdata(handles.fig_bh_tomo_pick, 'h', h)

end


function ActionMenu_Callback(hObject, eventdata, handles)


function pickMeanTraceMenuItem_Callback(hObject, eventdata, handles)
pick_mean_tr(handles)


function alignTracesMenuItem_Callback(hObject, eventdata, handles)
aligner(handles)


function precedentMenuItem_Callback(hObject, eventdata, handles)
precedent(handles)

function suivantMenuItem_Callback(hObject, eventdata, handles)
suivant(handles)


function showMeanTraceMenuItem_Callback(hObject, eventdata, handles)
mog = getappdata(handles.fig_bh_tomo_pick, 'mog');
h = getappdata(handles.fig_bh_tomo_pick, 'h');


theta = zeros(numel(h.bins),1);
for n=1:length(theta)
    theta(n) = h.bins{n}.theta_mid;
end

tracesM = zeros(length(theta), mog.data.nptsptrc);

for n=1:length(theta)
    if ( h.bins{n}.no_it_align == 0 )
        traces = (h.proc_data(:, h.bins{n}.no_trace))';
        [~, i] = max(h.snr(h.bins{n}.no_trace));
        if ( get(handles.checkbox_1cycle,'Value')==1 || ...
                get(handles.checkbox_scale_fdom,'Value')==1 ) && ...
                ~isempty(h.data_1cycle)
            traces = (h.data_1cycle(:, h.bins{n}.no_trace))';
        end
        tracesM(n,:) = traces(i,:);
    else
        traces = (h.proc_data(:, h.bins{n}.no_trace))';
        if ( get(handles.checkbox_1cycle,'Value')==1 || get(handles.checkbox_scale_fdom,'Value')==1 )...
                && get(handles.checkbox_showOrig,'Value')==0
            pds = str2double(get(handles.edit_poids1cycle,'String'));
            tr2 = (h.data_1cycle(:, h.bins{n}.no_trace))';
            traces = (1-pds)*traces + pds*tr2;
        end
        tracesM(n,:) = mean(traces);
    end
end

tt = h.pick;
%et = h.pick_et;

ttmin = min(tt(tt~=-1));
dt = mog.data.timec;

t = mog.data.timestp-ttmin;
for n=1:length(theta)
    AMax = max(abs(tracesM(n,:)));
    tracesM(n,:) = 5*tracesM(n,:)./AMax + theta(n);
    if ( tt(n)~=-1 )
        npts = round((tt(n)-ttmin)/dt);
        tracesM(n,:) = circshift(tracesM(n,:)', -npts)';
    end
end

figure
plot(tracesM(1,:), t);
hold on
for n=2:length(theta)
    plot(tracesM(n,:), t);
end
xlim = get(gca,'XLim');
plot(xlim, [0 0], 'r--')
hold off
ylim = get(gca,'YLim');
ylim(1) = -10;
ylim(2) = 80;
set(gca,'YDir','reverse','YLim',ylim);

    
    
    
