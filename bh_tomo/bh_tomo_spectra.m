function varargout = bh_tomo_spectra(varargin)
% BH_TOMO_SPECTRA M-file for bh_tomo_spectra.fig
%      BH_TOMO_SPECTRA, by itself, creates a new BH_TOMO_SPECTRA or raises the existing
%      singleton*.
%
%      H = BH_TOMO_SPECTRA returns the handle to a new BH_TOMO_SPECTRA or the handle to
%      the existing singleton*.
%
%      BH_TOMO_SPECTRA('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in BH_TOMO_SPECTRA.M with the given input arguments.
%
%      BH_TOMO_SPECTRA('Property','Value',...) creates a new BH_TOMO_SPECTRA or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before bh_tomo_spectra_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to bh_tomo_spectra_OpeningFcn via varargin.
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

% Edit the above text to modify the response to help bh_tomo_spectra

% Last Modified by GUIDE v2.5 14-Oct-2013 19:54:07

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @bh_tomo_spectra_OpeningFcn, ...
                   'gui_OutputFcn',  @bh_tomo_spectra_OutputFcn, ...
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


% --- Executes just before bh_tomo_spectra is made visible.
function bh_tomo_spectra_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to bh_tomo_spectra (see VARARGIN)

% Choose default command line output for bh_tomo_spectra
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);


h.mog = varargin{1};
h.Tx = unique(h.mog.data.Tx_z);
nos = cell(length(h.Tx),1);
for n = 1:length(h.Tx)
	nos{n} = num2str(n);
end
set(handles.popupmenu_groupe, 'String', nos)

h.rdata = detrend_rad(h.mog.data.rdata);
Amax = kron(max(abs(h.rdata)), ones(size(h.rdata,1),1));
h.rdata = h.rdata./Amax;

setappdata(handles.fig_bh_tomo_spectra, 'h', h)
update_figs(handles)

% UIWAIT makes bh_tomo_spectra wait for user response (see UIRESUME)
% uiwait(handles.fig_bh_tomo_spectra);


% --- Outputs from this function are returned to the command line.
function varargout = bh_tomo_spectra_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


function popupmenu_groupe_Callback(hObject, eventdata, handles)
update_figs(handles)

function popupmenu_groupe_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


function update_figs(handles)
h = getappdata(handles.fig_bh_tomo_spectra, 'h');

n = get(handles.popupmenu_groupe,'Value');
fmax = str2double(get(handles.edit_fmax,'String'));
SNRmax = str2double(get(handles.edit_max_snr,'String'));

dt = h.mog.data.timec * h.mog.fac_dt;

ind = find( h.Tx(n) == h.mog.data.Tx_z );
traces = h.rdata(:,ind);
%nfft = 2^nextpow2(size(traces,1));
nfft = 2^(1+nextpow2(size(traces,1)));

fac_f = 1.0;
fac_t = 1.0;
if strcmp(h.mog.data.tunits,'ns')
    % radar: freq nominale donnée en MHz
    fac_f = 10^6;
    fac_t = 10^-9;
elseif strcmp(h.mog.data.tunits,'ms')
    % sismique: on assume que f dominante est en kHz
    disp('warning : assuming Tx nominal frequency in _kHz_')
    fac_f = 10^3;
    fac_t = 10^-3;
else
    disp('warning : assuming Tx nominal frequency in _Hz_')
    disp('warning : assuming time step in _s_')
end

f0 = h.mog.data.rnomfreq * fac_f;
dt = dt * fac_t;

Fs = 1/dt;


% bruit pris sur dernière 20 ns
win_snr = round(20/h.mog.data.timec);
snr = data_select(traces, f0, dt, win_snr);

if get(handles.checkbox_lowpass,'Value') == 1

	HalfFs = Fs/2;
	Wp = 1.4*f0/HalfFs;
	Ws = 1.6*f0/HalfFs;
	Rp = 3;
	Rs = 40;
	[nc,Wn] = cheb1ord(Wp,Ws,Rp,Rs);

	[b,a] = cheby1(nc,0.5,Wn);
	
	for nt=1:size(traces,2)
%		if snr(nt) < 50
			traces(:,nt) = filtfilt(b,a,traces(:,nt));
%		end
	end
end

if get(handles.checkbox_fdom,'Value') == 1
	[tmp,freq] = pburg(traces(:,1),2,nfft,Fs);
	fdom = zeros(size(traces,2), 1);
	Pxx = zeros(length(tmp), size(traces,2));
	
	min_fdom = str2double(get(handles.edit_min_fdom,'String')) * fac_f;
	max_fdom = str2double(get(handles.edit_max_fdom,'String')) * fac_f;
	
	for nt=1:size(traces,2)
		[fdom(nt), Pxx(:,nt)] = get_fdom(traces(:,nt), dt, f0, min_fdom, max_fdom);
	end
else
	method = get(handles.popupmenu_method,'Value');
	if method == 1
		[tmp, freq] = pwelch(traces(:,1),[],[],[],Fs);
		Pxx = zeros(length(tmp), size(traces,2));
		Pxx(:,1) = tmp;
		for nt=2:size(traces,2)
			Pxx(:,nt) = pwelch(traces(:,nt),[],[],[],Fs);
		end
	else
		ordre = method;
		[tmp,freq] = pburg(traces(:,1),ordre,nfft,Fs);
		Pxx = zeros(length(tmp), size(traces,2));
		Pxx(:,1) = tmp;
		for nt=2:size(traces,2)
			Pxx(:,nt) = pburg(traces(:,nt),ordre,nfft,Fs);
		end
	end
end

z = h.mog.data.Rx_z(ind);

timestp = h.mog.data.timestp * h.mog.fac_dt;

axes(handles.axes1)
imagesc(timestp, z, traces')
xlabel(['Time [',h.mog.data.tunits,']'],'FontSize',12)
ylabel(['Rx elevation [',h.mog.data.cunits,']'],'FontSize',12)
title('Normalized amplitudes','FontSize',14)
set(handles.axes1,'YDir','normal')
	
axes(handles.axes2)
imagesc(freq/fac_f, z, log10(Pxx'))
xlabel('Frequency [MHz]','FontSize',12)
title('log_{10} Power spectra','FontSize',14)
set(handles.axes2, 'XLim', [0 fmax],'YDir','normal')
if get(handles.checkbox_fdom,'Value') == 1
	hold(handles.axes2,'on')
	ind = fdom~=-1;
	plot(handles.axes2, fdom(ind)/fac_f,z(ind),'ko')
	hold(handles.axes2,'off')
end

ylim = get(handles.axes1,'YLim');
plot(handles.axes3, snr, z)
%plot(handles.axes3, log10(snr), z)
%xlabel(handles.axes3,'log_{10}(snr)')
title(handles.axes3,'Signal-to-Noise ratio','FontSize',14)
set(handles.axes3, 'YLim',ylim,'YDir','reverse','XLim',[0 SNRmax])

set(handles.text_suptitle, 'String', ['Tx elevation: ',num2str(h.Tx(n)),' m'])


function edit_fmax_Callback(hObject, eventdata, handles)
fmax = str2double(get(hObject,'String'));
set(handles.axes2, 'XLim', [0 fmax])


function edit_fmax_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


function checkbox_lowpass_Callback(hObject, eventdata, handles)
update_figs(handles)


function popupmenu_method_Callback(hObject, eventdata, handles)
update_figs(handles)

function popupmenu_method_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


function checkbox_fdom_Callback(hObject, eventdata, handles)
update_figs(handles)

function edit_min_fdom_Callback(hObject, eventdata, handles)
update_figs(handles)

function edit_min_fdom_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


function edit_max_fdom_Callback(hObject, eventdata, handles)
update_figs(handles)

function edit_max_fdom_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


function edit_max_snr_Callback(hObject, eventdata, handles)
SNRmax = str2double(get(handles.edit_max_snr,'String'));
set(handles.axes3, 'YLim',ylim,'YDir','reverse','XLim',[0 SNRmax])


function edit_max_snr_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
