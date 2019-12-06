function varargout = bh_tomo_amp(varargin)
% BH_TOMO_AMP M-file for bh_tomo_amp.fig
%      BH_TOMO_AMP, by itself, creates a new BH_TOMO_AMP or raises the existing
%      singleton*.
%
%      H = BH_TOMO_AMP returns the handle to a new BH_TOMO_AMP or the handle to
%      the existing singleton*.
%
%      BH_TOMO_AMP('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in BH_TOMO_AMP.M with the given input arguments.
%
%      BH_TOMO_AMP('Property','Value',...) creates a new BH_TOMO_AMP or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before bh_tomo_amp_OpeningFunction gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to bh_tomo_amp_OpeningFcn via varargin.
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

% Edit the above text to modify the response to help bh_tomo_amp

% Last Modified by GUIDE v2.5 30-Apr-2019 14:49:38

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @bh_tomo_amp_OpeningFcn, ...
                   'gui_OutputFcn',  @bh_tomo_amp_OutputFcn, ...
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


% --- Executes just before bh_tomo_amp is made visible.
function bh_tomo_amp_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to bh_tomo_amp (see VARARGIN)

% Choose default command line output for bh_tomo_amp
handles.output = hObject;
handles.init_data = 0;

% Update handles structure
guidata(hObject, handles);

h.db_file = get(handles.fig_bh_amp,'UserData');
h.init_data = 0;
h.ffac = 1e6;
h.fUnits = 'MHz';
h.Fmin_centroid = 1*h.ffac;
h.Fmax_centroid = 175*h.ffac;
h.st_mean = false;
h.fen_App = 30;
h.lrai = [];
h.updateFcentroid_S = true;
h.updateFcentroid_fft = true;
h.thetaCut = 45;
h.A0 = 0;
h.z_surf = -0.01;
h.zmax_rais = [];
h.name_Ldc = {};
h.str = get_str_locale();
set(handles.f_min_centroide,'String',num2str(h.Fmin_centroid/h.ffac))
set(handles.f_max_centroide,'String',num2str(h.Fmax_centroid/h.ffac))
set(handles.fen_App,'String',num2str(h.fen_App))
set(handles.tCut,'String',num2str(h.thetaCut))
set(handles.z_surface,'String',num2str(h.z_surf))
setappdata(handles.fig_bh_amp, 'h', h)
set_String_locale(handles, h.str)

% UIWAIT makes bh_tomo_amp wait for user response (see UIRESUME)
% uiwait(handles.fig_bh_amp);


% --- Outputs from this function are returned to the command line.
function varargout = bh_tomo_amp_OutputFcn(hObject, eventdata, handles)
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;

% --- Executes during object creation, after setting all properties.
function no_trace_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



%
% --------------------------------------------------------------------
%
function no_trace_Callback(hObject, eventdata, handles)

h = getappdata(handles.fig_bh_amp, 'h');
no_trace = str2double(get(hObject,'string'));
if isnan(no_trace)
    errordlg(h.str.s54, h.str.s45,'modal')
    set(hObject,'string',num2str(h.no_trace))
	return
end
mog = getappdata(handles.fig_bh_amp, 'mog');
if ~mog.in(no_trace)
    errordlg('Trace was excluded after pruning',h.str.s45,'modal')
    set(hObject,'string',num2str(h.no_trace))
    return 
end
h.no_trace = no_trace;
setappdata(handles.fig_bh_amp, 'h', h)
update_tout(hObject, eventdata, handles)


%
% --------------------------------------------------------------------
%
function trace_suivante_Callback(hObject, eventdata, handles)
h = getappdata(handles.fig_bh_amp, 'h');
mog = getappdata(handles.fig_bh_amp, 'mog');
no = h.no_trace+1;
while no <= size(h.proc_data,2)
    if mog.in(no)
        break
    end
    no = no+1;
end
if no <= size(h.proc_data,2) && mog.in(no)
    h.no_trace = no;
    setappdata(handles.fig_bh_amp, 'h', h)
    update_tout(hObject, eventdata, handles)
end


%
% --------------------------------------------------------------------
%
function trace_precedente_Callback(hObject, eventdata, handles)
h = getappdata(handles.fig_bh_amp, 'h');
mog = getappdata(handles.fig_bh_amp, 'mog');
no = h.no_trace-1;
while no >= 1
    if mog.in(no)
        break
    end
    no = no-1;
end
if no >= 1 && mog.in(no)
    h.no_trace = no;
    setappdata(handles.fig_bh_amp, 'h', h)
    update_tout(hObject, eventdata, handles)
end


%
% --------------------------------------------------------------------
%
function prochaine_trace_Callback(hObject, eventdata, handles)
h = getappdata(handles.fig_bh_amp, 'h');
mog = getappdata(handles.fig_bh_amp, 'mog');
if get(handles.traite_tt_data,'Value')==1
    for n=1:size(h.proc_data,2)
        if h.data_pick_min(n) ~= -1 && ~mog.amp_done(n) && mog.in(n) %#ok<AND2>
            h.no_trace = n;
            break
        end
    end
else
    for n=1:size(h.proc_data,2)
        if ~mog.amp_done(n) && mog.in(n)
            h.no_trace = n;
            break
        end
    end
end
setappdata(handles.fig_bh_amp, 'h', h)
update_tout(hObject, eventdata, handles)

%
% --------------------------------------------------------------------
%
function fichier_Callback(hObject, eventdata, handles)

%
% --------------------------------------------------------------------
%
function ouvrir_Callback(hObject, eventdata, handles)
h = getappdata(handles.fig_bh_amp, 'h');

[h.no_mog,h.db_file,hh] = choisirMOG('UserData', h.db_file);
delete(hh)
if h.no_mog==0
	return
end
load(h.db_file,'mogs','boreholes')
mog = mogs(h.no_mog);

h.proc_data = mog.traces;
h.Fs = 1e9/mog.data.timec; % Hz
if strcmp(mog.data.tunits,'ms')
	h.Fs = 1e3/mog.data.timec; % Hz
end
h.data_dewow = [];
f = h.Fs * (0:round(mog.data.nptsptrc/2)-1)/mog.data.nptsptrc;
ind = f<=h.Fmax_centroid;
h.iFmax_centroid = length(f(ind));
ind = f<=h.Fmin_centroid;
h.iFmin_centroid = length(f(ind));
h.data_pick_min = mog.tt-mog.et;
%if ~isempty(mog.t0_merged)
%    h.data_pick_min = h.data_pick_min + mog.t0_merged;
%end
h.data_pick_min(mog.tt==-1) = -1;
ind = mog.amp_tmin == -1;
ind2 = mog.tt ~= -1;
mog.amp_tmin(ind&ind2) = mog.tt(ind&ind2)-mog.et(ind&ind2);
%if ~isempty(mog.t0_merged)
%    mog.amp_tmin(ind&ind2) = mog.amp_tmin(ind&ind2) + mog.t0_merged(ind&ind2);
%end

set(handles.edit_omega, 'String', mog.data.antennas( regexp(mog.data.antennas,'\d') ))

h.z_surf = min( [boreholes(mog.Tx).Z boreholes(mog.Rx).Z] ) - 0.001;
h.no_trace = 1;
h.init_data = 1;
set(handles.t_min,'String','0');
set(handles.t_max,'String', num2str(mog.data.timestp(length(mog.data.timestp))));
set(handles.z_surface,'String',num2str(h.z_surf))
setappdata(handles.fig_bh_amp, 'h', h)
setappdata(handles.fig_bh_amp, 'mog', mog)
get_rais_droits(handles)
set(handles.pushbutton_rc,'String',h.str.s154)
prochaine_trace_Callback(hObject, eventdata, handles)
%update_tout(hObject, eventdata, handles)

%
% --------------------------------------------------------------------
%
function fichier_inversion_Callback(hObject, eventdata, handles)
SaveMenuItem_Callback(hObject, eventdata, handles);

%
% --------------------------------------------------------------------
%
function type_analyse_Callback(hObject, eventdata, handles)
if (get(hObject,'Value')==1 || get(hObject,'Value')==2 ) && ...
		strcmp(get(handles.tSshowFMenuItem,'Checked'),'on')
	set(handles.pushbutton_fit_spectre,'Enable','on')
else
	set(handles.pushbutton_fit_spectre,'Enable','off')
end
update_spectre(handles)
update_positions_info(handles)

%
% --------------------------------------------------------------------
%
function quitter_Callback(hObject, eventdata, handles)
delete(handles.fig_bh_amp)

%
% --------------------------------------------------------------------
%
function update_tout(hObject, eventdata, handles)
h = getappdata(handles.fig_bh_amp, 'h');
set(handles.no_trace,'String',num2str(h.no_trace))
update_trace_simple(hObject, eventdata, handles)
update_spectre(handles)
update_t_lim(hObject, eventdata, handles)
calcule_App(hObject, eventdata, handles)
update_positions_info(handles)

%
% --------------------------------------------------------------------
%
function update_trace_simple(hObject, eventdata, handles)
h = getappdata(handles.fig_bh_amp, 'h');
mog = getappdata(handles.fig_bh_amp, 'mog');
axes(handles.trace)

if get(handles.filtre_ondelette,'Value')==1
  m = mean(mog.traces(:,h.no_trace));
  plot(mog.data.timestp, mog.traces(:,h.no_trace)-m, 'color',[0.8 0.8 0.8])
  hold on
  plot(mog.data.timestp, h.proc_data(:,h.no_trace))
else
  plot(mog.data.timestp, h.proc_data(:,h.no_trace))
  hold on
end
  
ylim=get(handles.trace,'YLim');
update_t_lim(hObject, eventdata, handles)
if mog.tt_done(h.no_trace) & h.data_pick_min(h.no_trace) ~= -1 %#ok<AND2>
    x = [h.data_pick_min(h.no_trace) ...
        h.data_pick_min(h.no_trace)+h.fen_App ...
        h.data_pick_min(h.no_trace)+h.fen_App ...
		h.data_pick_min(h.no_trace)];
    y = kron(ylim, [1 1]);
    p = patch(x, y, [0.85 0.325 0.098]);
    p.EdgeColor = 'none';
    p.FaceAlpha = 0.4;
end

grid on
xlabel([h.str.s22,' [',mog.data.tunits,']'])
ylabel(h.str.s21)
if length(mog.amp_tmin)>=h.no_trace
  plot([mog.amp_tmin(h.no_trace) mog.amp_tmin(h.no_trace)],ylim,'color',[0.466 0.674 0.118]);
end
if length(mog.amp_tmax)>=h.no_trace
  plot([mog.amp_tmax(h.no_trace) mog.amp_tmax(h.no_trace)],ylim,'color',[0.466 0.674 0.118]);
end
hold off

%
% --------------------------------------------------------------------
%
function update_t_lim(hObject, eventdata, handles)
tlim = [str2double(get(handles.t_min,'String')) str2double(get(handles.t_max,'String')) ];
set(handles.trace,'XLim',tlim)

%
% --------------------------------------------------------------------
%
function update_positions_info(handles)
h = getappdata(handles.fig_bh_amp, 'h');
mog = getappdata(handles.fig_bh_amp, 'mog');
Tx_x = mog.data.Tx_x(h.no_trace);
Tx_y = mog.data.Tx_y(h.no_trace);
Tx_z = mog.data.Tx_z(h.no_trace);
Rx_x = mog.data.Rx_x(h.no_trace);
Rx_y = mog.data.Rx_y(h.no_trace);
Rx_z = mog.data.Rx_z(h.no_trace);
no_db = mog.no_traces(h.no_trace);
texte{1} = '';
texte{2} = ['MOG: ',mog.name];
texte{3} = '';
texte{4} = ['Positions Tx -- Rx (db # ',num2str(no_db),')'];
texte{5} = '';
texte{6} = '       x        y        z';
texte{7} = ['Tx:   ',num2str(Tx_x,'%5.2f'),'   ',num2str(Tx_y, '%5.2f'),'   ',num2str(Tx_z,'%5.2f')];
texte{8} = ['Rx:   ',num2str(Rx_x,'%5.2f'),'   ',num2str(Rx_y, '%5.2f'),'   ',num2str(Rx_z,'%5.2f')];
texte{9} = [num2str(sum(mog.in)), ' traces, (', num2str(size(h.proc_data,2)),' total)'];%[num2str(size(h.proc_data,2)),' traces'];
texte{10} = '';
if get(handles.type_analyse,'Value')==1 || get(handles.type_analyse,'Value')==2
    texte{11} = [h.str.s87,': ',num2str(round(mog.fcentroid(h.no_trace)/h.ffac)),' MHz'];
    texte{12} = [h.str.s88,': ',num2str(round(mog.scentroid(h.no_trace)/(h.ffac^2))),' MHz'];
elseif get(handles.type_analyse,'Value')==3
	texte{11} = [h.str.s83,': ',num2str(round(mog.App(h.no_trace)))];
else
	texte{11} = [h.str.s243,': ',num2str(round(mog.App(h.no_trace)))];
end
set(handles.info,'String',texte)


%
% --------------------------------------------------------------------
%
function t_min_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



%
% --------------------------------------------------------------------
%
function t_min_Callback(hObject, eventdata, handles)
update_t_lim(hObject, eventdata, handles)


%
% --------------------------------------------------------------------
%
function t_max_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


%
% --------------------------------------------------------------------
%
function t_max_Callback(hObject, eventdata, handles)
update_t_lim(hObject, eventdata, handles)


%
% --------------------------------------------------------------------
%
function sauve_mat(handles)
hh=msgbox('Saving Data');
h = getappdata(handles.fig_bh_amp, 'h');
mog = getappdata(handles.fig_bh_amp, 'mog');
load(h.db_file,'mogs','air')
mogs(h.no_mog) = mog; %#ok<NASGU>
save(h.db_file,'mogs','-append')
if ishandle(hh)==1
	delete(hh)
end



%
% --------------------------------------------------------------------
%
function pointe(hObject, eventdata, handles)
get_tmin_tmax(hObject, eventdata, handles)
update_tout(hObject, eventdata, handles)


%
% --------------------------------------------------------------------
%
function get_tmin_tmax(hObject, eventdata, handles)
h = getappdata(handles.fig_bh_amp, 'h');
mog = getappdata(handles.fig_bh_amp, 'mog');
h.updateFcentroid_S = true;
h.updateFcentroid_fft = true;

button = 1; %#ok<NASGU>
%tcmp = 0;
%acmp = 1;
pointe_actif=1;
xlim = get(handles.trace,'XLim');
ylim = get(handles.trace,'YLim');
n_trace_traite = 0;
tmin = mog.amp_tmin(h.no_trace);
tmax = mog.amp_tmax(h.no_trace);
update_tout(hObject, eventdata, handles)
while pointe_actif
  [t1,a,button] = ginput(1);
  if gca ~= handles.trace || t1<xlim(1) || t1>xlim(2) || a<ylim(1) || a>ylim(2)
		%pointe_actif=0;
		setappdata(handles.fig_bh_amp, 'h', h)
		return
  end
  t=t1;
  if button==1 
		if tmin==-1
			tmin = t1;
			tmax = t1;
			update_images_while_picking(hObject, eventdata, handles, tmin, tmax)
			[t2,a,button] = ginput(1);
			if gca ~= handles.trace || t2<xlim(1) || t2>xlim(2) || a<ylim(1) || a>ylim(2)
				%pointe_actif=0;
				setappdata(handles.fig_bh_amp, 'h', h)
				return
			end
			if t2>t1, tmax = t2; else tmin=t2; end
			t=t2;
		else
			if t1>tmin
				tmax=t1;
			else
				tmax=tmin;
				tmin=t1;
			end
		end
  end
  while button==1
		if t<=tmin
			tmin=t;
		elseif t>=tmax
			%
			%ii = findnear(t,mog.data.timestp);
			%m=
			tmax=t;
		elseif t-tmin < tmax-t
			tmin=t;
		else
			tmax=t;
		end
		update_images_while_picking(hObject, eventdata, handles, tmin, tmax)
		[t,a,button] = ginput(1);
		if gca ~= handles.trace || t<xlim(1) || t>xlim(2) || a<ylim(1) || a>ylim(2)
			%pointe_actif=0;
			setappdata(handles.fig_bh_amp, 'h', h)
			return
		end
  end
  if button==2
		mog.amp_done(h.no_trace) = true;
		mog.amp_tmin(h.no_trace) = tmin;
		mog.amp_tmax(h.no_trace) = tmax;
  elseif button==3
		mog.amp_done(h.no_trace) = true;
		mog.amp_tmin(h.no_trace) = -1;
		mog.amp_tmax(h.no_trace) = -1;
  end
  n_trace_traite = n_trace_traite+1;
  h.no_trace = h.no_trace+1;
  if h.no_trace<size(h.proc_data,2)
		
		if get(handles.traite_tt_data,'Value')==0
			while mog.amp_done(h.no_trace)
				h.no_trace = h.no_trace+1;
				if h.no_trace>length(mog.amp_done)
					break
				end
			end
		elseif sum(mog.amp_done | h.data_pick_min==-1)~=length(mog.amp_done)
			while ( mog.amp_done(h.no_trace) == true | ...
							h.data_pick_min(h.no_trace) == -1 ) %#ok<OR2>
				h.no_trace = h.no_trace+1;
				if h.no_trace>length(mog.amp_done)
					break
				end
			end
		else
			while ( h.data_pick_min(h.no_trace) == -1 )
				h.no_trace = h.no_trace+1;
				if h.no_trace>length(mog.amp_done)
					break
				end
			end
		end
  end
  if h.no_trace>size(h.proc_data,2)
		h.no_trace = h.no_trace-1;
		pointe_actif=0;
  end
  set(handles.no_trace,'String',num2str(h.no_trace))
  setappdata(handles.fig_bh_amp, 'h', h)
  setappdata(handles.fig_bh_amp, 'mog', mog)
  update_tout(hObject, eventdata, handles)
  tmin = mog.amp_tmin(h.no_trace);
  tmax = mog.amp_tmax(h.no_trace);
  if tmax == -1, tmax = tmin+h.fen_App; end
  ylim = get(handles.trace,'YLim');
end


%
% -------------------------------------------------------------------------
%
function pushbutton_reset_trace_Callback(hObject, eventdata, handles)
h = getappdata(handles.fig_bh_amp, 'h');
mog = getappdata(handles.fig_bh_amp, 'mog');
h.updateFcentroid_S = true;
h.updateFcentroid_fft = true;
mog.amp_done(h.no_trace) = false;
mog.amp_tmin(h.no_trace) = -1;
mog.amp_tmax(h.no_trace) = -1;
setappdata(handles.fig_bh_amp, 'h', h)
setappdata(handles.fig_bh_amp, 'mog', mog)
update_tout(hObject, eventdata, handles)


%
% --------------------------------------------------------------------
%
function pushbutton_manual_pick_Callback(hObject, eventdata, handles)
pointe(hObject, eventdata, handles)

%
% --------------------------------------------------------------------
%
function update_spectre(handles)
axes(handles.spectre)
cla
h = getappdata(handles.fig_bh_amp, 'h');
mog = getappdata(handles.fig_bh_amp, 'mog');
if (mog.amp_done(h.no_trace))
  if mog.amp_tmax(h.no_trace) > mog.amp_tmin(h.no_trace)
	plot_spectre(handles, ...
				 mog.amp_tmin(h.no_trace), ...
				 mog.amp_tmax(h.no_trace))
  end
end

%
% --------------------------------------------------------------------
%
function update_images_while_picking(hObject, eventdata, handles, tmin, tmax)

%h = getappdata(handles.fig_bh_amp, 'h');
update_trace_simple(hObject, eventdata, handles)
axes(handles.trace)
ylim = get(gca,'YLim');
hold on
plot([tmin tmin],ylim,'r');
plot([tmax tmax],ylim,'r');
hold off

if tmax>tmin
  plot_spectre(handles, tmin, tmax)
  update_positions_info(handles)
end

%
% --------------------------------------------------------------------
%
function plot_spectre(handles, tmin, tmax)
if get(handles.type_analyse,'Value')==1
	plot_fft(handles, tmin, tmax)
elseif get(handles.type_analyse,'Value')==2
	plot_S(handles, tmin, tmax)
end

%
% --------------------------------------------------------------------
%
function calcule_App(hObject, eventdata, handles)
h = getappdata(handles.fig_bh_amp, 'h');
mog = getappdata(handles.fig_bh_amp, 'mog');

imin = findnear(mog.amp_tmin(h.no_trace), mog.data.timestp);
if numel(imin)>1, imin=imin(1); end
imax = findnear(mog.amp_tmax(h.no_trace), mog.data.timestp);
if numel(imax)>1, imax=imax(1); end

if imax<=imin
    return
end

data = h.proc_data(:,h.no_trace);
mog.App(h.no_trace) = max(data(imin:imax)) - min(data(imin:imax));
setappdata(handles.fig_bh_amp, 'mog', mog)

%
% --------------------------------------------------------------------
%
function plot_fft(handles, tmin, tmax)
h = getappdata(handles.fig_bh_amp, 'h');
mog = getappdata(handles.fig_bh_amp, 'mog');

imin = findnear(tmin, mog.data.timestp);
imax = findnear(tmax, mog.data.timestp);
if numel(imin)>1, imin=imin(1); end
if numel(imax)>1, imax=imax(1); end

data = h.proc_data(:,h.no_trace);
%fenetre = hanning(mog.data.nptsptrc);

%data = h.proc_data(imin:imax,h.no_trace);
fenetre = zeros(mog.data.nptsptrc,1);
fenetre(imin:imax) = hanning(imax-imin+1);

data = data.*fenetre;
%figure(1)
%plot(data);
S = abs(fft(data));
N = length(S);
S = S(1:round(N/2))';
f = h.Fs * (0:length(S)-1)/N;
ind = ( f<=h.Fmax_centroid & f>=h.Fmin_centroid );
f = f(ind);
S = S(ind);
mog.fcentroid(h.no_trace) = sum(f.*S)/sum(S);
mog.scentroid(h.no_trace) = sum(((f-mog.fcentroid(h.no_trace)).^2).*S)/sum(S);

axes(handles.spectre)
plot(f/h.ffac,S)
ylim=get(gca,'YLim');
hold on
plot([mog.fcentroid(h.no_trace) mog.fcentroid(h.no_trace)]./h.ffac,ylim,'g-.')
hold off
ylabel(h.str.s21)
xlabel([h.str.s77,' [',h.fUnits,']'])

setappdata(handles.fig_bh_amp, 'mog', mog)


%
% --------------------------------------------------------------------
%
function plot_S(handles, tmin, tmax)
h = getappdata(handles.fig_bh_amp, 'h');
mog = getappdata(handles.fig_bh_amp, 'mog');

imin = findnear(tmin, mog.data.timestp);
imax = findnear(tmax, mog.data.timestp);
if numel(imin)>1, imin=imin(1); end
if numel(imax)>1, imax=imax(1); end

data = h.proc_data(:,h.no_trace);
[Pxx,t,f] = st_bg(data, h.iFmin_centroid, h.iFmax_centroid, 1/h.Fs);
Pxx = abs(Pxx);
t=t*1e9;
if strcmp(mog.data.tunits,'ms')
	t = t*1e3;
end
if h.st_mean
  S = Pxx(:,imin:imax);
  F=repmat(f',1,1+imax-imin);
  mog.fcentroid(h.no_trace) = mean(sum(F.*S)./sum(S));
  mog.scentroid(h.no_trace) = mean(sum(((F-mog.fcentroid(h.no_trace)).^2).*S)./sum(S));
else
%   [Pmax,ind] = max(Pxx);
%   ind=1:length(ind);
%   [Pmax,ind] = max(Pmax(ind>=imin & ind<=imax));
%   S = Pxx(:,ind+imin-1)';
%   tS = t(ind+imin-1);
  [S, tS] = getS(handles, data, imin, imax, Pxx, t);
  mog.fcentroid(h.no_trace) = sum(f.*S)/sum(S);
  mog.scentroid(h.no_trace) = sum(((f-mog.fcentroid(h.no_trace)).^2).*S)/sum(S);
end

axes(handles.spectre)
if strcmp(get(handles.tSshowTFMenuItem,'Checked'),'on')
	contourf(t,f/h.ffac,Pxx)
	xlim=get(gca,'XLim');
	ylim=get(gca,'YLim');
	hold on
	plot([tmin tmin],ylim,'r')
	plot([tmax tmax],ylim,'r')
	plot(xlim,[mog.fcentroid(h.no_trace) mog.fcentroid(h.no_trace)]./h.ffac,'g-.')
	if ~h.st_mean
		plot([tS tS],ylim,'c--')
	end
	hold off
	set(gca,'Xlim',[mog.data.timestp(1) mog.data.timestp(length(mog.data.timestp))]);
	xlabel([h.str.s22,' [',mog.data.tunits,']'])
	ylabel([h.str.s77,' [',h.fUnits,']'])
	colormap(flipud(cmr))
else
	plot(f/h.ffac,S)
	ylim=get(gca,'YLim');
	hold on
	plot([mog.fcentroid(h.no_trace) mog.fcentroid(h.no_trace)]/h.ffac,ylim,'g-.')
	hold off
	ylabel(h.str.s21)
	xlabel([h.str.s77,' [',h.fUnits,']'])
end
if ~h.st_mean
	axes(handles.trace)
	hold on
	ylim=get(gca,'YLim');
	plot([tS tS],ylim,'c--')
	hold off
end
setappdata(handles.fig_bh_amp, 'mog', mog)


%
% --------------------------------------------------------------------
%
function SaveMenuItem_Callback(hObject, eventdata, handles)
if get(handles.checkbox_hybride,'Value')==0
	update_Amp(handles)
	ajusteA0(hObject, eventdata, handles)
else
	ajusteHybride(hObject, eventdata, handles)
end
sauve_mat(handles)

%
% --------------------------------------------------------------------
%
function HelpMenuItem_Callback(hObject, eventdata, handles)


%
% --------------------------------------------------------------------
%
function [tau,data,m,et]=prepare_spectral(hObject, eventdata, handles)
update_fcentroid(handles, false);
h = getappdata(handles.fig_bh_amp, 'h');
mog = getappdata(handles.fig_bh_amp, 'mog');

ind = getInd(handles);
sigma_s2 = median(mog.scentroid(ind));
%sigma_s2 = 10;

f_Tx = h.ffac*str2double(get(handles.edit_omega,'String'));
data = f_Tx * mog.fcentroid(ind) / sigma_s2;
lrai = (h.lrai(ind))';

if get(handles.ponderationSB,'Value')==1
  w = calculSB(handles);
  et = 1./w;
  m = lscov([lrai ones(size(lrai))], data', w);
else
  m = [lrai ones(size(lrai))]\data';
  et = zeros(size(data));
end

A0 = m(2);
tau = A0 - data;

%
% --------------------------------------------------------------------
%
function [tau,Acorr,m,et]=prepare_App(hObject, eventdata, handles)
update_Amp(handles);
h = getappdata(handles.fig_bh_amp, 'h');
mog = getappdata(handles.fig_bh_amp, 'mog');

ind = getInd(handles);

Acorr = log( mog.App(ind).*h.lrai(ind) ./ (diagRayon(h.thetaTx(ind),handles)...
	.* diagRayon(h.thetaRx(ind),handles)) );

if get(handles.ponderationSB,'Value')==1
  w = calculSB(handles);
  et = 1./w;
  lrai = h.lrai(ind);
  m = lscov([lrai' ones(size(lrai'))], Acorr', w);
else
  m = [h.lrai(ind)' ones(size(h.lrai(ind)'))]\Acorr';
  et = zeros(size(Acorr));
end

A0 = m(2);
tau = A0 - Acorr;


%
% --------------------------------------------------------------------
%
function [tau,data,Acorr,mF,mA]=prepare_hybrid(handles)
h = getappdata(handles.fig_bh_amp, 'h');

if strcmp(get(handles.hybFftMenuItem,'Checked'),'on') && h.updateFcentroid_fft == true;
	update_fcentroid_fft(handles, flag)
elseif h.updateFcentroid_S == true;
	update_fcentroid_S(handles, flag)
end
if strcmp(get(handles.hybPtPMenuItem,'Checked'),'on')
	update_App(handles)
else
	update_Amax(handles)
end

mog = getappdata(handles.fig_bh_amp, 'mog');

ind = getInd(handles);
f_Tx = h.ffac*str2double(get(handles.edit_omega,'String'));
sigma_s2 = median(mog.scentroid(ind));
data = f_Tx * mog.fcentroid(ind) / sigma_s2;
lrai = (h.lrai(ind))';

mF = [lrai ones(size(lrai))]\data';
tau = mF(2) - data;
Acorr = log( mog.App(ind).*h.lrai(ind) ./ (diagRayon(h.thetaTx(ind),handles)...
	.* diagRayon(h.thetaRx(ind),handles)) );
mA = [h.lrai(ind)' ones(size(h.lrai(ind)'))]\Acorr';
tau = [tau mA(2) - Acorr];


%
% --------------------------------------------------------------------
%
function fl = rashkovskij(f0, K)
fl = 0.70710678118655*f0*sqrt(1+1/K);

%
% --------------------------------------------------------------------
%
function r = spectre_ricker(w, w0)
r = 2/(sqrt(pi)*w0^3) * w.^2 .* exp( -w.^2/w0^2 );

%
% --------------------------------------------------------------------
%
function dr = spectre_d_ricker(w, w0)
dr = 4/sqrt(pi) * (w/w0^3 - w.^3/w0^5) .* exp( -w.^2/w0^2 );


%
% --------------------------------------------------------------------
%
function stats_Callback(hObject, eventdata, handles)
if get(handles.checkbox_hybride,'Value')==1
	stats_hybrid(handles);
elseif get(handles.type_analyse,'Value')==1 || get(handles.type_analyse,'Value')==2
    stats_fcentroid(hObject, eventdata, handles);
else
    stats_App(hObject, eventdata, handles);
end

%
% --------------------------------------------------------------------
%
function stats_App(hObject, eventdata, handles)

[tau,Acorr,m,et]=prepare_App(hObject, eventdata, handles);

h = getappdata(handles.fig_bh_amp, 'h');
mog = getappdata(handles.fig_bh_amp, 'mog');

ind = getInd(handles);

h_stat = figure;
set(h_stat, 'Position',[349 69 565 671]);

subplot(321)
plot(h.lrai(ind),mog.App(ind),'o')
xlabel(h.str.s86)
ylabel(h.str.s83)

if get(handles.ponderationSB,'Value')==1
	subplot(322)
	plot(h.lrai(ind), 1./et, 'o')
	xlabel(h.str.s86)
	ylabel(h.str.s84)
end

subplot(323)
plot(h.lrai(ind), Acorr, 'o')
xlim=get(gca,'XLim');
hold on
plot(xlim,[xlim' [1 1]']*m,'g-')
hold off
xlabel(h.str.s86)
ylabel(h.str.s85)
title(['A_0 = ',num2str(m(2))])

subplot(324)
plot(h.lrai(ind),tau./h.lrai(ind),'o')
xlabel(h.str.s86)
ylabel('\alpha_{app}')

subplot(325)
plot(h.lrai(ind),tau,'o')
xlabel(h.str.s86)
ylabel('\tau')

subplot(326)
plot(h.thetaRais(ind),tau./h.lrai(ind),'o')
xlabel(h.str.s85)
ylabel('\alpha_{app}')

suptitle([num2str(sum(ind)),' traces'])

%
% --------------------------------------------------------------------
%
function stats_fcentroid(hObject, eventdata, handles)
%[tau,data,m,et]=prepare_spectral(hObject, eventdata, handles);
[tau,data,m]=prepare_spectral(hObject, eventdata, handles);

h = getappdata(handles.fig_bh_amp, 'h');
mog = getappdata(handles.fig_bh_amp, 'mog');
ind = getInd(handles);

h_stat = figure;
set(h_stat, 'Position',[349 69 565 671]);

subplot(321)
plot(h.lrai(ind), mog.fcentroid(ind)/h.ffac,'o')
xlabel(h.str.s86)
ylabel([h.str.s80,' [',h.fUnits,']'])

subplot(322)
plot(h.lrai(ind),mog.scentroid(ind)/(h.ffac^2),'o')
xlabel(h.str.s86)
ylabel([h.str.s81,' [',h.fUnits,'^2]'])

subplot(323)
plot(h.lrai(ind), data,'o')
xlim=get(gca,'XLim');
hold on
plot(xlim,[xlim' [1 1]']*m,'g-')
xlabel(h.str.s86)
ylabel([h.str.s82,' [',h.fUnits,'^{-1}]'])
sigma_s2 = median(mog.scentroid(ind));
f_Tx = h.ffac*str2double(get(handles.edit_omega,'string'));
title(['f_0 = ',num2str(m(2)*sigma_s2/(h.ffac*f_Tx))])

subplot(324)
plot(h.lrai(ind), tau./h.lrai(ind),'o')
xlabel(h.str.s86)
ylabel('\alpha_{app}')

subplot(325)
plot(h.lrai(ind), tau,'o')
xlabel(h.str.s86)
ylabel('\tau')

%
% --------------------------------------------------------------------
%
function pointe_traces_tt_Callback(hObject, eventdata, handles)


%
% --------------------------------------------------------------------
%
function f_max_centroide_Callback(hObject, eventdata, handles)
h = getappdata(handles.fig_bh_amp, 'h');
mog = getappdata(handles.fig_bh_amp, 'mog');
h.Fmax_centroid = h.ffac*str2double(get(hObject,'String'));
f = h.Fs * (0:round(mog.data.nptsptrc/2)-1)/mog.data.nptsptrc;
ind = f<=h.Fmax_centroid;
h.iFmax_centroid = length(f(ind));
setappdata(handles.fig_bh_amp, 'h', h)
if ~isempty(mog.amp_tmin)
    update_spectre(handles)
    update_positions_info(handles)
end


%
% --------------------------------------------------------------------
%
function f_max_centroide_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


%
% --------------------------------------------------------------------
%
function traite_tt_data_Callback(hObject, eventdata, handles)

%
% --------------------------------------------------------------------
%
function ponderationSB_Callback(hObject, eventdata, handles)

%
% --------------------------------------------------------------------
%
function SB = calculSB(handles)
%h = getappdata(handles.fig_bh_amp, 'h');
mog = getappdata(handles.fig_bh_amp, 'mog');

ind = getInd(handles);
no = 1:mog.data.ntrace;
no=no(ind);

fac_f = 1.0;
fac_t = 1.0;
if strcmp(mog.data.tunits,'ns')
    % radar: freq nominale donn�e en MHz
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

freq = mog.data.rnomfreq * fac_f;
dt = mog.data.timec * fac_t;

traces = mog.traces(:, no);
% bruit pris sur derni�re 20 ns
win_snr = round(20/mog.data.timec);
SB = data_select(traces, freq, dt, win_snr);


%
% --------------------------------------------------------------------
%
function update_fcentroid(handles, flag)
h = getappdata(handles.fig_bh_amp, 'h');

if get(handles.type_analyse,'Value')==1 && h.updateFcentroid_fft == true;
	update_fcentroid_fft(handles, flag)
elseif (get(handles.type_analyse,'Value')==2 || get(handles.checkbox_hybride,'Value')==1) && h.updateFcentroid_S == true;
	update_fcentroid_S(handles, flag)
end

%
% --------------------------------------------------------------------
%
function update_Amp(handles)
if get(handles.type_analyse,'Value')==3
	update_App(handles)
elseif get(handles.type_analyse,'Value')==4
	update_Amax(handles)
end

%
% --------------------------------------------------------------------
%
function update_App(handles)
h = getappdata(handles.fig_bh_amp, 'h');
mog = getappdata(handles.fig_bh_amp, 'mog');

ind = getInd(handles);
no = 1:mog.data.ntrace;
no=no(ind);

for n=1:length(no)
    data = h.proc_data(:,no(n));
    
	if mog.amp_tmax(no(n)) ~= -1
      % tmin must be picked
      imin = findnear(mog.amp_tmin(no(n)), mog.data.timestp);
	  imax = findnear(mog.amp_tmax(no(n)), mog.data.timestp);
	  if numel(imax)>1, imax=imax(1); end
    else
      imin = findnear(h.data_pick_min(no(n)), mog.data.timestp);
      if numel(imin)>1, imin=imin(1); end
	  imax = imin + round(h.fen_App/mog.data.timec);
	end
    mog.App(no(n)) = max(data(imin:imax)) - min(data(imin:imax));
end
setappdata(handles.fig_bh_amp, 'mog', mog)


%
% --------------------------------------------------------------------
%
function update_Amax(handles)
h = getappdata(handles.fig_bh_amp, 'h');
mog = getappdata(handles.fig_bh_amp, 'mog');

ind = getInd(handles);
no = 1:mog.data.ntrace;
no=no(ind);

for n=1:length(no)
    data = h.proc_data(:,no(n));
    imin = findnear(h.data_pick_min(no(n)), mog.data.timestp);
    if numel(imin)>1, imin=imin(1); end
	if mog.amp_tmax(no(n)) ~= -1
	  imax = findnear(mog.amp_tmax(no(n)), mog.data.timestp);
	  if numel(imax)>1, imax=imax(1); end
	else
	  imax = imin + round(h.fen_App/mog.data.timec);
	end
    mog.App(no(n)) = max(abs(data(imin:imax)));
end
setappdata(handles.fig_bh_amp, 'mog', mog)


%
% --------------------------------------------------------------------
%
function update_fcentroid_fft(handles, flag)
h = getappdata(handles.fig_bh_amp, 'h');
mog = getappdata(handles.fig_bh_amp, 'mog');

ind = getInd(handles);
no = 1:mog.data.ntrace;
no=no(ind);
if ~flag
    hh=waitbar(0, h.str.s78);
end
for n=1:length(no)
    if mog.fcentroid(no(n)) ~= 0 && flag
		disp(n)
        continue
    end
    if ~flag
        waitbar(n/length(no));
    end
    data = h.proc_data(:,no(n));
    imin = findnear(mog.amp_tmin(no(n)), mog.data.timestp);
    imax = findnear(mog.amp_tmax(no(n)), mog.data.timestp);
    if numel(imin)>1, imin=imin(1); end
    if numel(imax)>1, imax=imax(1); end

	fenetre = zeros(mog.data.nptsptrc,1);
	fenetre(imin:imax) = hanning(imax-imin+1);
	
	S = abs(fft(data.*fenetre));
	N = length(S);
	S = S(1:round(N/2))';
	f = h.Fs * (0:length(S)-1)/N;
%    ind = f<=h.Fmax_centroid;
	ind = ( f<=h.Fmax_centroid & f>=h.Fmin_centroid );
    f = f(ind);
    S = S(ind);
	mog.fcentroid(no(n)) = sum(f.*S)/sum(S);
	mog.scentroid(no(n)) = sum(((f-mog.fcentroid(no(n))).^2).*S)/sum(S);
end
if ~flag
    close(hh);
end
h.updateFcentroid_fft = false;
setappdata(handles.fig_bh_amp, 'h', h)
setappdata(handles.fig_bh_amp, 'mog', mog)


%
% --------------------------------------------------------------------
%
function update_fcentroid_S(handles, flag)
h = getappdata(handles.fig_bh_amp, 'h');
mog = getappdata(handles.fig_bh_amp, 'mog');

ind = getInd(handles);
no = 1:mog.data.ntrace;
no = no(ind);
if ~flag
    hh=waitbar(0, h.str.s78);
end
for n=1:length(no)
    if mog.fcentroid(no(n)) ~= 0 && flag
		disp(n)
        continue
    end
    if ~flag
        waitbar(n/length(no));
    end
    data = h.proc_data(:,no(n));
	imin = findnear(mog.amp_tmin(no(n)), mog.data.timestp);
    imax = findnear(mog.amp_tmax(no(n)), mog.data.timestp);
    if numel(imin)>1, imin=imin(1); end
    if numel(imax)>1, imax=imax(1); end
	[Pxx,t,f] = st_bg(data, h.iFmin_centroid, h.iFmax_centroid, 1/h.Fs);
	Pxx = abs(Pxx);
	t=t*1000;
	if h.st_mean
	  S = Pxx(:,imin:imax);
	  F=repmat(f',1,1+imax-imin);
	  mog.fcentroid(no(n)) = mean(sum(F.*S)./sum(S));
	  mog.scentroid(no(n)) = mean(sum(((F-mog.fcentroid(no(n))).^2).*S)./sum(S));
	else
% 	  [Pmax,ind] = max(Pxx);
% 	  ind=1:length(ind);
% 	  [Pmax,ind] = max(Pmax(ind>=imin & ind<=imax));
% 	  S = Pxx(:,ind+imin-1)';

      %[S, tS] = getS(handles, data, imin, imax, Pxx, t);
	  S = getS(handles, data, imin, imax, Pxx, t);

	  mog.fcentroid(no(n)) = sum(f.*S)/sum(S);
	  mog.scentroid(no(n)) = sum(((f-mog.fcentroid(no(n))).^2).*S)/sum(S);
	end
end
if ~flag
    close(hh);
end
h.updateFcentroid_S = false;
setappdata(handles.fig_bh_amp, 'h', h)
setappdata(handles.fig_bh_amp, 'mog', mog)

%
% --------------------------------------------------------------------
%
function edit_omega_Callback(hObject, eventdata, handles)


%
% --------------------------------------------------------------------
%
function edit_omega_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

%
% --------------------------------------------------------------------
%
function Amax = getAmax(handles)
mog = getappdata(handles.fig_bh_amp, 'mog');

ind = getInd(handles);
no = 1:mog.data.ntrace;
no=no(ind);

Amax = -1*ones(size(no));

for n=1:length(no)
    imin = findnear(mog.amp_tmin(no(n)), mog.data.timestp);
    imax = findnear(mog.amp_tmax(no(n)), mog.data.timestp);
    if numel(imin)>1, imin=imin(1); end
    if numel(imax)>1, imax=imax(1); end
    Amax(no(n)) = max(h.proc_data(imin:imax,no(n)));
end


%
% --------------------------------------------------------------------
%
function ind = getInd(handles)
h = getappdata(handles.fig_bh_amp, 'h');
mog = getappdata(handles.fig_bh_amp, 'mog');
ind2 = (mog.amp_tmin ~= -1);
ind = mog.amp_done & ind2;
ind = ind & (~isnan(h.lrai)) & mog.amp_tmin<mog.amp_tmax;

ind2 = ( h.thetaTx > (h.thetaCut*pi/180) & ...
		 h.thetaTx < ((180-h.thetaCut)*pi/180) );

ind = ind(:) & ind2(:);
ind2 = ( h.thetaRx > h.thetaCut*pi/180 & ...
		 h.thetaRx < (180-h.thetaCut)*pi/180 );
ind = ind(:) & ind2(:);

ind = ind & mog.in(:);

if ~isempty(h.zmax_rais)
    ind2 = h.zmax_rais < h.z_surf;
    ind = ind(:) & ind2(:);
end

%
% --------------------------------------------------------------------
%
function fen_App_Callback(hObject, eventdata, handles)
h = getappdata(handles.fig_bh_amp, 'h');
h.fen_App = str2double(get(hObject,'String'));
setappdata(handles.fig_bh_amp, 'h', h)
update_trace_simple(hObject, eventdata, handles)
update_positions_info(hObject)


%
% --------------------------------------------------------------------
%
function fen_App_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


%
% --------------------------------------------------------------------
%
function filtre_ondelette_Callback(hObject, eventdata, handles)
h = getappdata(handles.fig_bh_amp, 'h');
mog = getappdata(handles.fig_bh_amp, 'mog');
if ~h.init_data
    warndlg(h.str.s52, h.str.s53)
    set(hObject,'Value',0);
    return
end
if get(hObject,'Value')==0
    h.proc_data = mog.traces;
else
  if isempty( mog.fw )
    mog.fw = filtrage_wavelet2(mog.traces);
  end
  h.proc_data = detrend_rad(mog.fw);
end
setappdata(handles.fig_bh_amp,'h',h)
setappdata(handles.fig_bh_amp,'mog',mog)
update_tout(hObject, eventdata, handles)


%
% --------------------------------------------------------------------
%
function get_rais_droits(handles)
h = getappdata(handles.fig_bh_amp, 'h');
mog = getappdata(handles.fig_bh_amp, 'mog');
h.lrai = sqrt((mog.data.Tx_x-mog.data.Rx_x).^2+...
	(mog.data.Tx_y-mog.data.Rx_y).^2+...
	(mog.data.Tx_z-mog.data.Rx_z).^2);
%dz = mog.data.Rx_z-mog.data.Tx_z;
%h.thetaTx = pi/2-asin(dz./h.lrai); % assume boreholes verticaux

a = (mog.data.Rx_x-mog.data.Tx_x)./h.lrai;
b = (mog.data.Rx_y-mog.data.Tx_y)./h.lrai;
c = (mog.data.Rx_z-mog.data.Tx_z)./h.lrai;
h.thetaTx = acos(dot([a' b' c'], mog.TxCosDir,2))';
h.thetaRx = h.thetaTx;
h.thetaRais = h.thetaTx;
setappdata(handles.fig_bh_amp,'h',h)


%
% --------------------------------------------------------------------
%
function get_rais_courbes(handles,tomo)
h = getappdata(handles.fig_bh_amp, 'h');
mog = getappdata(handles.fig_bh_amp, 'mog');
ind1 = false(size(mog.no_traces));
ind2 = [];
for n=1:length(ind1)
    ii = find( mog.no_traces(n)==tomo.no_trace );
    if ~isempty(ii)
        ind1(n) = true;
        ind2 = [ind2 ii]; %#ok<AGROW>
    end
end
if isempty(ind2)
    uiwait(warndlg(h.str.s188))
    get_rais_droits(handles)
    return
end

h.lrai = nan(1,mog.data.ntrace);
h.thetaTx = h.lrai;
h.thetaRx = h.lrai;
h.zmax_rais = h.lrai;
%h.lrai(ind1) = sum(tomo.L(ind2,:)')';
h.lrai(ind1) = sum(tomo.L(ind2,:), 2)';

dim = size(tomo.rays{1}, 2);
nn=1;
for n=1:length(ind1)
    if ind1(n)
        % Tx
        n2 = ind2(nn);
        nn=nn+1;
        dz = tomo.rays{n2}(3,dim)-tomo.rays{n2}(1,dim);
        hyp = sqrt(sum( (tomo.rays{n2}(3,:)-tomo.rays{n2}(1,:)).^2 ));
        h.thetaTx(n) = pi/2-asin(dz/hyp);
        % Rx
        nseg = size(tomo.rays{n2},1);
        dz = tomo.rays{n2}(nseg,dim)-tomo.rays{n2}(nseg-2,dim);
        hyp = sqrt(sum( (tomo.rays{n2}(nseg,:)-tomo.rays{n2}(nseg-2,:)).^2 ));
        h.thetaRx(n) = pi/2-asin(dz/hyp);
        % z max
        h.zmax_rais(n) = max( tomo.rays{n2}(:,dim) );
    end
end

setappdata(handles.fig_bh_amp,'h',h)

%
% --------------------------------------------------------------------
%
function type_rais_SelectionChangeFcn(hObject, eventdata, handles)

%
% --------------------------------------------------------------------
%
function f_min_centroide_Callback(hObject, eventdata, handles)
h = getappdata(handles.fig_bh_amp, 'h');
mog = getappdata(handles.fig_bh_amp, 'mog');
h.Fmin_centroid = h.ffac*str2double(get(hObject,'String'));
f = h.Fs * (0:round(mog.data.nptsptrc/2)-1)/mog.data.nptsptrc;
ind = f<=h.Fmin_centroid;
h.iFmin_centroid = length(f(ind));
setappdata(handles.fig_bh_amp, 'h', h)
if ~isempty(mog.amp_tmin)
    update_spectre(handles)
    update_positions_info(handles)
end

%
% --------------------------------------------------------------------
%
function f_min_centroide_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



%
% --------------------------------------------------------------------
%
function tCut_Callback(hObject, eventdata, handles)
h = getappdata(handles.fig_bh_amp, 'h');
h.thetaCut = str2double(get(hObject,'String'));
setappdata(handles.fig_bh_amp, 'h', h)



%
% --------------------------------------------------------------------
%
function tCut_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end




%
% --------------------------------------------------------------------
%
function ajusteA0(hObject, eventdata, handles)
h = getappdata(handles.fig_bh_amp, 'h');
mog = getappdata(handles.fig_bh_amp, 'mog');
ind = getInd(handles);
h2.lrai = h.lrai(ind);
if get(handles.type_analyse,'Value')==1 || get(handles.type_analyse,'Value')==2
    [h2.tau,h2.data,h2.m,h2.et]=prepare_spectral(hObject, eventdata, handles);
	h2.rawdata = mog.fcentroid(ind);
	h2.var = mog.scentroid(ind);
	h2.label{1} = [h.str.s80,' [',h.fUnits,']'];  %YH add h.
	h2.label{2} = [h.str.s81,' [',h.fUnits,'^2]']; %YH add h.
	h2.label{3} = [h.str.s82,' [',h.fUnits,'^{-1}]']; %YH add h.
	h2.type = 'fce';
	h2.f_Tx = h.ffac*str2double(get(handles.edit_omega,'String'));
elseif get(handles.type_analyse,'Value')==3 || get(handles.type_analyse,'Value')==4
    [h2.tau,h2.data,h2.m,h2.et]=prepare_App(hObject, eventdata, handles);
	h2.rawdata = mog.App(ind);
	h2.var = 1./h2.et;
	h2.label{1} = h.str.s83;
	h2.label{2} = h.str.s84;
	h2.label{3} = h.str.s85;
	h2.type = 'amp';
end
h2.label{4} = '\alpha_{app}';
h2.label{5} = '\tau';
h2.label{6} = h.str.s86;

[A0,hh]=fitA0(h2);
close(hh)

if get(handles.type_analyse,'Value')==1 || get(handles.type_analyse,'Value')==2
	mog.tauFce = -1*ones(size(mog.amp_done));
	mog.tauFce_et = -1*ones(size(mog.amp_done));
	mog.tauFce(ind) = A0 - h2.data;
	mog.tau_params.f0 = A0;
	mog.tauFce_et(ind) = h2.et;
elseif get(handles.type_analyse,'Value')==3 || get(handles.type_analyse,'Value')==4
	mog.tauApp = -1*ones(size(mog.amp_done));
	mog.tauApp_et = -1*ones(size(mog.amp_done));
	mog.tauApp(ind) = A0 - h2.data;
	mog.tau_params.A0 = A0;
	mog.tauApp_et(ind) = h2.et;
end

mog.amp_name_Ldc = h.name_Ldc;
setappdata(handles.fig_bh_amp, 'mog', mog)


%
% --------------------------------------------------------------------
%
function ajusteHybride(hObject, eventdata, handles)
if strcmp(get(handles.hybStMenuItem,'Checked'), 'on')
	update_fcentroid_S(handles, false);
else
	update_fcentroid_fft(handles, false);
end
if strcmp(get(handles.hybPtPMenuItem,'Checked'), 'on')
	update_App(handles);
else
	update_Amax(handles);
end
h = getappdata(handles.fig_bh_amp, 'h');
mog = getappdata(handles.fig_bh_amp, 'mog');

ind = getInd(handles);
h2.f = h.ffac*str2double(get(handles.edit_omega,'String'));
h2.fcentroid = mog.fcentroid(ind);
h2.scentroid = mog.scentroid(ind);
h2.lrai = (h.lrai(ind))';
h2.Acorr = log( mog.App(ind).*h.lrai(ind) ./ (diagRayon(h.thetaTx(ind),handles)...
	.* diagRayon(h.thetaRx(ind),handles)) );
h2.name_Ldc = h.name_Ldc;
h2.db_file = h.db_file;
h2.no_traces = mog.no_traces(ind);
h2.Tx = [mog.data.Tx_x(ind)' mog.data.Tx_y(ind)' mog.data.Tx_z(ind)'];
h2.ffac = h.ffac;
h2.fUnits = h.fUnits;

[tau,hh,params] = fitHybride(h2, mog.tau_params);
close(hh)
if tau == 0, return, end
if params.poids_App == 1
	mog.tauApp = -1*ones(size(mog.amp_done));
	mog.tauApp_et = -1*ones(size(mog.amp_done));
	mog.tauApp(ind) = tau;
elseif params.poids_fc == 1
	mog.tauFce = -1*ones(size(mog.amp_done));
	mog.tauFce_et = -1*ones(size(mog.amp_done));
	mog.tauFce(ind) = tau;
else
	mog.tauHyb = -1*ones(size(mog.amp_done));
	mog.tauHyb_et = -1*ones(size(mog.amp_done));
	mog.tauHyb(ind) = tau;
end
mog.tau_params = params;
mog.amp_name_Ldc = h.name_Ldc;
setappdata(handles.fig_bh_amp, 'mog', mog)


%
% --------------------------------------------------------------------
%
function z_surface_Callback(hObject, eventdata, handles)
h = getappdata(handles.fig_bh_amp, 'h');
h.z_surf = str2double(get(hObject,'String'));
setappdata(handles.fig_bh_amp, 'h', h)


%
% --------------------------------------------------------------------
%
function z_surface_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


%
% --------------------------------------------------------------------
%
function update_Data(handles)
if get(handles.type_analyse,'Value')==1 || get(handles.type_analyse,'Value')==2
	update_fcentroid(handles, false);
else
	update_Amp(handles);
end



%
% -------------------------------------------------------------------------
%
function set_String_locale(handles, str)

h = getappdata(handles.fig_bh_amp, 'h');
set(handles.type_analyse,          'String', str.s63)
set(handles.pushbutton_manual_pick,'String', str.s67)
set(handles.no_trace_label,        'String', str.s24)
set(handles.trace_suivante,        'String', str.s16)
set(handles.trace_precedente,      'String', str.s15)
set(handles.prochaine_trace,       'String', str.s17)
set(handles.pushbutton_reset_trace,'String', str.s18)
set(handles.stats,                 'String', str.s11)
set(handles.ponderationSB,         'String', str.s68)
set(handles.traite_tt_data,        'String', str.s69)
set(handles.uipanel1,              'Title',  str.s13)
set(handles.omega_label,           'String', [str.s70,' [',h.fUnits,']'])
set(handles.Fmin_centroid_label,   'String', str.s71)
set(handles.Fmax_centroid_label,   'String', str.s72)
set(handles.fenetre_App,           'String', str.s73)
set(handles.filtre_ondelette,      'String', str.s19)
set(handles.fichier,               'Label',  str.s25)
set(handles.ouvrir,                'Label',  str.s26)
set(handles.SaveMenuItem,          'Label',  str.s29)
set(handles.quitter,               'Label',  str.s31)
set(handles.HelpMenuItem,          'Label',  str.s32)
set(handles.pushbutton_rc,         'String', str.s154)
set(handles.checkbox_utiliser_fen, 'String', str.s155)
set(handles.checkbox_hybride,      'String', str.s65)
set(handles.showHelpMenuItem,      'Label',  str.s79)
set(handles.diagrammeMenuItem,     'Label',  str.s167)
set(handles.tSMenuItem,            'Label',  str.s178)
set(handles.tSshowFMenuItem,       'Label',  str.s238)
set(handles.tSshowTFMenuItem,      'Label',  str.s239)
set(handles.HybrideMenuItem,       'Label',  str.s240)
set(handles.hybPtPMenuItem,        'Label',  str.s241)
set(handles.hybAbsMenuItem,        'Label',  str.s242)
set(handles.hybStMenuItem,         'Label',  str.s244)
set(handles.hybFftMenuItem,        'Label',  str.s245)
set(handles.stMaxEMenuItem,        'Label',  str.s246)
set(handles.st1MenuItem,           'Label',  str.s247)
set(handles.stMaxAMenuItem,        'Label',  str.s248)
set(handles.stFbMenuItem,          'Label',  str.s249)
set(handles.maxwellMenuItem,       'Label',  str.s254)
set(handles.rayleighMenuItem,      'Label',  str.s255)
set(handles.normaleMenuItem,       'Label',  str.s256)
set(handles.logNormalMenuItem,     'Label',  str.s257)
set(handles.fitSMenuItem,          'Label',  str.s258)
set(handles.pushbutton_fit_spectre,'String', str.s264)
set(handles.fUnitsMenuItem,        'Label',  str.s267)

%set(handles.,      'String', str.s)


%
% -------------------------------------------------------------------------
%
function checkbox_utiliser_fen_Callback(hObject, eventdata, handles)


%
% -------------------------------------------------------------------------
%
function pushbutton_rc_Callback(hObject, eventdata, handles)
h = getappdata(handles.fig_bh_amp, 'h');
[tomo, h.name_Ldc, hh] = get_tomo_Ldc(h.db_file);
delete(hh)
setappdata(handles.fig_bh_amp, 'h', h)
if ~isempty(tomo)
    set(hObject,'String',[h.name_Ldc{1},'-',h.name_Ldc{2}])
    get_rais_courbes(handles, tomo)
else
    set(hObject,'String','Load curved rays')
    get_rais_droits(handles);
end

%
% -------------------------------------------------------------------------
%
function pushbutton_rc_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


%
% -------------------------------------------------------------------------
%
function checkbox_hybride_Callback(hObject, eventdata, handles)


%
% -------------------------------------------------------------------------
%
function stats_hybrid(handles)
[tau,data,Acorr,mF,mA]=prepare_hybrid(handles);

h = getappdata(handles.fig_bh_amp, 'h');
mog = getappdata(handles.fig_bh_amp, 'mog');

ind = getInd(handles);
h_stat = figure;
set(h_stat, 'Position',[349 69 565 671]);

subplot(321)
plot(h.lrai(ind), Acorr, 'go')
xlim=get(gca,'XLim');
hold on
plot(xlim,[xlim' [1 1]']*mA,'-')
hold off
xlabel(h.str.s86)
ylabel(h.str.s85)
title(['A_0 = ',num2str(mA(2))])

subplot(322)
plot(h.lrai(ind), data,'o')
xlim=get(gca,'XLim');
hold on
plot(xlim,[xlim' [1 1]']*mF,'g-')
xlabel(h.str.s86)
ylabel([h.str.s82,' [',h.fUnits,'^{-1}]'])
sigma_s2 = median(mog.scentroid(ind));
f_Tx = h.ffac*str2double(get(handles.edit_omega,'string'));
title(['f_0 = ',num2str(mF(2)*sigma_s2/(h.ffac*f_Tx))])

subplot(323)
plot(h.lrai(ind) ,tau(1:numel(tau)/2)./h.lrai(ind),'o')
hold on
plot(h.lrai(ind) ,tau(1+numel(tau)/2:end)./h.lrai(ind),'go')
hold off
xlabel(h.str.s86)
ylabel('\alpha_{app}')

subplot(324)
plot(h.lrai(ind),tau(1:numel(tau)/2),'o')
hold on
plot(h.lrai(ind) ,tau(1+numel(tau)/2:end),'go')
hold off
xlabel(h.str.s86)
ylabel('\tau')


function OptionsMenu_Callback(hObject, eventdata, handles)


function diagrammeMenuItem_Callback(hObject, eventdata, handles)


function sinDiagMenuItem_Callback(hObject, eventdata, handles)
if strcmp(get(hObject,'Checked'),'on')
	set(hObject,'Checked','off')
	set(handles.sin2DiagMenuItem,'Checked','on')
else
	set(hObject,'Checked','on')
	set(handles.sin2DiagMenuItem,'Checked','off')
end


function sin2DiagMenuItem_Callback(hObject, eventdata, handles)
if strcmp(get(hObject,'Checked'),'on')
	set(hObject,'Checked','off')
	set(handles.sinDiagMenuItem,'Checked','on')
else
	set(hObject,'Checked','on')
	set(handles.sinDiagMenuItem,'Checked','off')
end


function corr = diagRayon(theta,handles)
if strcmp(get(handles.sinDiagMenuItem,'Checked'),'on')
	corr = sin(theta);
elseif strcmp(get(handles.sin2DiagMenuItem,'Checked'),'on')
	corr = sin(theta).^2;
end


function showHelpMenuItem_Callback(hObject, eventdata, handles)
h = getappdata(handles.fig_bh_amp, 'h');
helpdlg(h.str.s75, h.str.s32)


function tSMenuItem_Callback(hObject, eventdata, handles)


function tSshowFMenuItem_Callback(hObject, eventdata, handles)
if strcmp(get(hObject,'Checked'),'on')
	set(hObject,'Checked','off')
	set(handles.tSshowTFMenuItem,'Checked','on')
	set(handles.pushbutton_fit_spectre,'Enable','off')
else
	set(hObject,'Checked','on')
	set(handles.tSshowTFMenuItem,'Checked','off')
	if get(handles.type_analyse,'Value')==1 || get(handles.type_analyse,'Value')==2
		set(handles.pushbutton_fit_spectre,'Enable','on')
	end
end
guidata(hObject, handles);
update_spectre(handles)
update_positions_info(handles)


function tSshowTFMenuItem_Callback(hObject, eventdata, handles)
if strcmp(get(hObject,'Checked'),'on')
	set(hObject,'Checked','off')
	set(handles.tSshowFMenuItem,'Checked','on')
	if get(handles.type_analyse,'Value')==1 || get(handles.type_analyse,'Value')==2
		set(handles.pushbutton_fit_spectre,'Enable','on')
	end
else
	set(hObject,'Checked','on')
	set(handles.tSshowFMenuItem,'Checked','off')
	set(handles.pushbutton_fit_spectre,'Enable','off')
end
set(handles.pushbutton_fit_spectre,'Enable','off')
guidata(hObject, handles);
update_spectre(handles)
update_positions_info(handles)


function HybrideMenuItem_Callback(hObject, eventdata, handles)


function hybPtPMenuItem_Callback(hObject, eventdata, handles)
if strcmp(get(hObject,'Checked'),'on')
	set(hObject,'Checked','off')
	set(handles.hybAbsMenuItem,'Checked','on')
else
	set(hObject,'Checked','on')
	set(handles.hybAbsMenuItem,'Checked','off')
end


function hybAbsMenuItem_Callback(hObject, eventdata, handles)
if strcmp(get(hObject,'Checked'),'on')
	set(hObject,'Checked','off')
	set(handles.hybPtPMenuItem,'Checked','on')
else
	set(hObject,'Checked','on')
	set(handles.hybPtPMenuItem,'Checked','off')
end


function hybStMenuItem_Callback(hObject, eventdata, handles)
if strcmp(get(hObject,'Checked'),'on')
	set(hObject,'Checked','off')
	set(handles.hybFftMenuItem,'Checked','on')
else
	set(hObject,'Checked','on')
	set(handles.hybFftMenuItem,'Checked','off')
end


function hybFftMenuItem_Callback(hObject, eventdata, handles)
if strcmp(get(hObject,'Checked'),'on')
	set(hObject,'Checked','off')
	set(handles.hybStMenuItem,'Checked','on')
else
	set(hObject,'Checked','on')
	set(handles.hybStMenuItem,'Checked','off')
end


function stMaxEMenuItem_Callback(hObject, eventdata, handles)
if strcmp(get(hObject,'Checked'),'on')
	set(hObject,'Checked','off')
	set(handles.st1MenuItem,'Checked','on')
	set(handles.stMaxAMenuItem,'Checked','off')
	set(handles.stFbMenuItem,'Checked','off')
else
	set(hObject,'Checked','on')
	set(handles.st1MenuItem,'Checked','off')
	set(handles.stFbMenuItem,'Checked','off')
	set(handles.stMaxAMenuItem,'Checked','off')
end
update_spectre(handles)
update_positions_info(handles)



function st1MenuItem_Callback(hObject, eventdata, handles)
if strcmp(get(hObject,'Checked'),'on')
	set(hObject,'Checked','off')
	set(handles.stMaxEMenuItem,'Checked','on')
	set(handles.stFbMenuItem,'Checked','off')
	set(handles.stMaxAMenuItem,'Checked','off')
else
	set(hObject,'Checked','on')
	set(handles.stMaxEMenuItem,'Checked','off')
	set(handles.stFbMenuItem,'Checked','off')
	set(handles.stMaxAMenuItem,'Checked','off')
end
update_spectre(handles)
update_positions_info(handles)



function stMaxAMenuItem_Callback(hObject, eventdata, handles)
if strcmp(get(hObject,'Checked'),'on')
	set(hObject,'Checked','off')
	set(handles.stMaxEMenuItem,'Checked','off')
	set(handles.st1MenuItem,'Checked','on')
	set(handles.stFbMenuItem,'Checked','off')
else
	set(hObject,'Checked','on')
	set(handles.stMaxEMenuItem,'Checked','off')
	set(handles.st1MenuItem,'Checked','off')
	set(handles.stFbMenuItem,'Checked','off')
end
update_spectre(handles)
update_positions_info(handles)



function stFbMenuItem_Callback(hObject, eventdata, handles)
if strcmp(get(hObject,'Checked'),'on')
	set(hObject,'Checked','off')
	set(handles.stMaxEMenuItem,'Checked','off')
	set(handles.st1MenuItem,'Checked','on')
	set(handles.stMaxAMenuItem,'Checked','off')
else
	set(hObject,'Checked','on')
	set(handles.stMaxEMenuItem,'Checked','off')
	set(handles.st1MenuItem,'Checked','off')
	set(handles.stMaxAMenuItem,'Checked','off')
end
update_spectre(handles)
update_positions_info(handles)




function [S, tS] = getS(handles, trace, imin, imax, Pxx, t)
if strcmp(get(handles.stFbMenuItem,'Checked'),'on')
	S = Pxx(:,imin)';
	tS = t(imin);
elseif strcmp(get(handles.stMaxEMenuItem,'Checked'),'on')
	[Pmax,ind] = max(Pxx);
	ind=1:length(ind);
	[Pmax,ind] = max(Pmax(ind>=imin & ind<=imax));
	S = Pxx(:,ind+imin-1)';
	tS = t(ind+imin-1);
elseif strcmp(get(handles.st1MenuItem,'Checked'),'on')
	[Amin, iAmin] = min(trace(imin:imax));
	[Amax, iAmax] = max(trace(imin:imax));
	if iAmin(1)<iAmax(1)
		S = Pxx(:,imin+iAmin(1)-1)';
		tS = t(imin+iAmin(1)-1);
	else
		S = Pxx(:,imin+iAmax(1)-1)';
		tS = t(imin+iAmax(1)-1);
	end
elseif strcmp(get(handles.stMaxAMenuItem,'Checked'),'on')
	[Amax, iAmax] = max(abs(trace(imin:imax)));
	S = Pxx(:,imin+iAmax(1)-1)';
	tS = t(imin+iAmax(1)-1);
end

function pushbutton_fit_spectre_Callback(hObject, eventdata, handles)
%%%%%%%%%%%%%%YH
h = getappdata(handles.fig_bh_amp, 'h');   %add to get data
hf=get(handles.spectre,'Children');   % h--->hf
if isempty(hf)  %add
    msgbox('Please pick up a trace.','Message','modal');
    return;
end
for n=1:length(hf)   
	if get(hf(n),'Color')==[0 0 1] %#ok<BDSCA> % on cherche le spectre trac? en bleu
		break
	end
end
f = get(hf(n),'XData');   
S = get(hf(n),'YData');   
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if strcmp(get(handles.maxwellMenuItem,'Checked'),'on')
	type = 1;
	coul = 'r';
elseif strcmp(get(handles.rayleighMenuItem,'Checked'),'on')
	type = 2;
	coul = 'c';
elseif strcmp(get(handles.gammaMenuItem,'Checked'),'on')
	type = 3;
	coul = 'm';
else
	return;
end
 [S2, param] = fitS(f*h.ffac, S, type);
param.a
axes(handles.spectre)
hold on
plot(f,S2,coul)
hold off



function maxwellMenuItem_Callback(hObject, eventdata, handles)
set(hObject,'Checked','on')
set(handles.normaleMenuItem,'Checked','off')
set(handles.rayleighMenuItem,'Checked','off')
set(handles.gammaMenuItem,'Checked','off')
set(handles.logNormalMenuItem,'Checked','off')


function normaleMenuItem_Callback(hObject, eventdata, handles)
set(hObject,'Checked','on')
set(handles.maxwellMenuItem,'Checked','off')
set(handles.rayleighMenuItem,'Checked','off')
set(handles.gammaMenuItem,'Checked','off')
set(handles.logNormalMenuItem,'Checked','off')


function rayleighMenuItem_Callback(hObject, eventdata, handles)
set(hObject,'Checked','on')
set(handles.maxwellMenuItem,'Checked','off')
set(handles.gammaMenuItem,'Checked','off')
set(handles.normaleMenuItem,'Checked','off')
set(handles.logNormalMenuItem,'Checked','off')


function logNormalMenuItem_Callback(hObject, eventdata, handles)
set(hObject,'Checked','on')
set(handles.maxwellMenuItem,'Checked','off')
set(handles.rayleighMenuItem,'Checked','off')
set(handles.gammaMenuItem,'Checked','off')
set(handles.normaleMenuItem,'Checked','off')


function gammaMenuItem_Callback(hObject, eventdata, handles)
set(hObject,'Checked','on')
set(handles.maxwellMenuItem,'Checked','off')
set(handles.rayleighMenuItem,'Checked','off')
set(handles.normaleMenuItem,'Checked','off')
set(handles.logNormalMenuItem,'Checked','off')


function fitSMenuItem_Callback(hObject, eventdata, handles)


function fUnitsMenuItem_Callback(hObject, eventdata, handles)


function mhzMenuItem_Callback(hObject, eventdata, handles)
set(hObject,'Checked','on')
set(handles.khzMenuItem,'Checked','off')
set(handles.hzMenuItem,'Checked','off')
h = getappdata(handles.fig_bh_amp, 'h');
h.fUnits = 'MHz';
h.ffac = 1e6;
setappdata(handles.fig_bh_amp, 'h', h)
set(handles.omega_label,           'String', [h.str.s70,' [',h.fUnits,']'])
update_tout(hObject, eventdata, handles)



function khzMenuItem_Callback(hObject, eventdata, handles)
set(hObject,'Checked','on')
set(handles.mhzMenuItem,'Checked','off')
set(handles.hzMenuItem,'Checked','off')
h = getappdata(handles.fig_bh_amp, 'h');
h.fUnits = 'kHz';
h.ffac = 1e3;
setappdata(handles.fig_bh_amp, 'h', h)
set(handles.omega_label,           'String', [h.str.s70,' [',h.fUnits,']'])
update_tout(hObject, eventdata, handles)



function hzMenuItem_Callback(hObject, eventdata, handles)
set(hObject,'Checked','on')
set(handles.khzMenuItem,'Checked','off')
set(handles.mhzMenuItem,'Checked','off')
h = getappdata(handles.fig_bh_amp, 'h');
h.fUnits = 'Hz';
h.ffac = 1;
setappdata(handles.fig_bh_amp, 'h', h)
set(handles.omega_label,           'String', [h.str.s70,' [',h.fUnits,']'])
update_tout(hObject, eventdata, handles)


function pushbutton_autopick_Callback(hObject, eventdata, handles)

if get(handles.type_analyse,'Value') ~= 3
    warndlg('Autopicker working only for peak-to-peak amplitude')
    return
end
mog = getappdata(handles.fig_bh_amp, 'mog');

% reset variables
mog.amp_tmin(:) = -1;
mog.amp_tmax(:) = -1;
mog.App(:) = 0;
mog.amp_done(:) = false;

dt = mog.data.timestp(2) - mog.data.timestp(1);
for n=1:mog.data.ntrace
    if mog.tt_done(n)
        [Ama, tma] = findpeaks(mog.traces(:,n), mog.data.timestp);
        [Ami, tmi] = findpeaks(-mog.traces(:,n), mog.data.timestp);
        [~, i] = min(abs(tma-mog.tt(n)));
        tmax = tma(i);
        Amax = Ama(i);
        ind = tmi>tmax;  % we search for Amax then Amin which comes _after_
        tmi = tmi(ind);
        [~, i] = min(abs(tmi-mog.tt(n)));
        tmin = tmi(i);
        Amin = -Ami(i);
        
        mog.amp_tmin(n) = tmax-dt; % start of window is just before Amax
		mog.amp_tmax(n) = tmin+dt; % end of window is just after Amin
        mog.App(n) = Amax - Amin;
        mog.amp_done(n) = 1;
    end
end
setappdata(handles.fig_bh_amp, 'mog', mog)
update_tout(hObject, eventdata, handles)
  
  
