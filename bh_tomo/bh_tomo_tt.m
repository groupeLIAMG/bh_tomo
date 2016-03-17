function varargout = bh_tomo_tt(varargin)
% BH_TOMO_TT M-file for bh_tomo_tt.fig
%      BH_TOMO_TT, by itself, creates a new BH_TOMO_TT or raises the existing
%      singleton*.
%
%      H = BH_TOMO_TT returns the handle to a new BH_TOMO_TT or the handle to
%      the existing singleton*.
%
%      BH_TOMO_TT('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in BH_TOMO_TT.M with the given input arguments.
%
%      BH_TOMO_TT('Property','Value',...) creates a new BH_TOMO_TT or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before bh_tomo_tt_OpeningFunction gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to bh_tomo_tt_OpeningFcn via varargin.
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

% Edit the above text to modify the response to help bh_tomo_tt

% Last Modified by GUIDE v2.5 16-Mar-2016 11:31:29

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @bh_tomo_tt_OpeningFcn, ...
                   'gui_OutputFcn',  @bh_tomo_tt_OutputFcn, ...
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

%
% -------------------------------------------------------------------------
%
function bh_tomo_tt_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to bh_tomo_tt (see VARARGIN)

% Choose default command line output for bh_tomo_tt
handles.output = hObject;
handles.init_data = 0;
% Update handles structure
guidata(hObject, handles);

h.db_file = get(handles.fig_bh_tomo_tt,'UserData');
h.init_data = false;
h.init_av = false;
h.init_ap = false;
h.saved = true;
h.str = get_str_locale();
setappdata(handles.fig_bh_tomo_tt, 'h', h)
set_String_locale(handles, h.str)

if exist('swt','file') ~=2
    set(handles.filtre_ondelette,'Enable','off')
end


% UIWAIT makes bh_tomo_tt wait for user response (see UIRESUME)
% uiwait(handles.fig_bh_tomo_tt);

%
% --- Outputs from this function are returned to the command line.
%
function varargout = bh_tomo_tt_OutputFcn(hObject, eventdata, handles)
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;

%
% -------------------------------------------------------------------------
%
function FileMenu_Callback(hObject, eventdata, handles)


%
% -------------------------------------------------------------------------
%
function OpenMenuItem_Callback(hObject, eventdata, handles)
h = getappdata(handles.fig_bh_tomo_tt, 'h');

[h.no_mog,h.db_file,hh] = choisirMOG('UserData', h.db_file);
delete(hh)
if h.no_mog==0
	return
end
set(handles.fig_bh_tomo_tt, 'Name', ['bh_tomo_tt - ',h.db_file])

load(h.db_file,'mogs','air')
mog = mogs(h.no_mog);
h.init_data = true;
clear mogs
if isempty(mog.av)
	av = [];
	set(handles.radiobutton_av,'Enable','off')
	h.init_av = false;
else
	av = air(mog.av);
	set(handles.radiobutton_av,'String',[h.str.s02,': ',av.name])
	set(handles.radiobutton_av,'Enable','on')
	h.init_av = true;
end
if isempty(mog.ap)
	ap = [];
	set(handles.radiobutton_ap,'Enable','off')
	h.init_ap = false;
else
	ap = air(mog.ap);
	set(handles.radiobutton_ap,'String',[h.str.s03,': ',ap.name])
	set(handles.radiobutton_ap,'Enable','on')
	h.init_ap = true;
end

set(handles.radiobutton_data,'String',[h.str.s23,': ',mog.name])

%h.proc_data = detrend_rad(mog.data.rdata);
%cmin=-max(max(abs(h.proc_data)));
%cmax= -cmin;
%h.cminmax = [cmin cmax];

h.data_vrp = false;
if ~isempty(findstr(mog.data.csurvmod, 'VRP'))
  h.data_vrp = true;
end

set(handles.t_max,'String',num2str(round(mog.data.timestp(length(mog.data.timestp)))))

h.pickTx = mog.ttTx;
h.doneTx = mog.ttTx_done;
setappdata(handles.fig_bh_tomo_tt, 'h', h)
setappdata(handles.fig_bh_tomo_tt, 'mog', mog)
setappdata(handles.fig_bh_tomo_tt, 'av', av)
setappdata(handles.fig_bh_tomo_tt, 'ap', ap)

if ~isempty(mog.data.tdata)
    set(handles.checkbox_pickTx,'Enable','on')
else
    set(handles.checkbox_pickTx,'Enable','off')
end

set(handles.filtre_ondelette,'Value',0)
set(handles.radiobutton_data,'Value',1)
radiobutton_data_Callback(hObject, eventdata, handles)



%
% -------------------------------------------------------------------------
%
function PrintMenuItem_Callback(hObject, eventdata, handles)
printdlg(handles.fig_bh_tomo_tt)


%
% -------------------------------------------------------------------------
%
function CloseMenuItem_Callback(hObject, eventdata, handles)
%selection = questdlg(['Close ' get(handles.fig_bh_tomo_tt,'Name') '?'],...
%                     ['Close ' get(handles.fig_bh_tomo_tt,'Name') '...'],...
%                     'Yes','No','Yes');
%if strcmp(selection,'No')
%    return;
%end
h=getappdata(handles.fig_bh_tomo_tt, 'h');
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
delete(handles.fig_bh_tomo_tt)


%
% -------------------------------------------------------------------------
%
function action_CreateFcn(hObject, eventdata, handles)
if ispc
    set(hObject,'BackgroundColor','white');
else
    set(hObject,'BackgroundColor',get(0,'defaultUicontrolBackgroundColor'));
end


%
% -------------------------------------------------------------------------
%
function action_Callback(hObject, eventdata, handles)

%
% -------------------------------------------------------------------------
%
function trace_simple_ButtonDownFcn(hObject, eventdata, handles)


%
% -------------------------------------------------------------------------
%
function rb_data_Callback(hObject, eventdata, handles)

%
% -------------------------------------------------------------------------
%
function rb_av_Callback(hObject, eventdata, handles)

%
% -------------------------------------------------------------------------
%
function rb_ap_Callback(hObject, eventdata, handles)

%
% -------------------------------------------------------------------------
%
function mutual_exclude(off)
set(off,'Value',0)

%
% -------------------------------------------------------------------------
%
function trace_traite_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


%
% -------------------------------------------------------------------------
%
function trace_traite_Callback(hObject, eventdata, handles)
h = getappdata(handles.fig_bh_tomo_tt, 'h');
mog = getappdata(handles.fig_bh_tomo_tt, 'mog');
av = getappdata(handles.fig_bh_tomo_tt, 'av');
ap = getappdata(handles.fig_bh_tomo_tt, 'ap');
no_trace = str2double(get(hObject,'string'));
if isnan(no_trace)
    errordlg(h.str.s54,h.str.s45,'modal')
	set(hObject,'string',num2str(h.no_trace))
	return
end
nos = 1:mog.data.ntrace;
nos = nos(mog.in);
if isempty( find(no_trace==nos,1) )
	errordlg('Trace was excluded after pruning',h.str.s45,'modal')
	set(hObject,'string',num2str(h.no_trace))
	return
elseif (get(handles.radiobutton_data,'Value')==1) && no_trace>mog.data.ntrace
    errordlg([h.str.s55,num2str(mog.data.ntrace)],h.str.s45,'modal')
	set(hObject,'string',num2str(h.no_trace))
	return
elseif (get(handles.radiobutton_av,'Value')==1) && no_trace > av.data.ntrace
    errordlg([h.str.s55,num2str(av.data.ntrace)],h.str.s45,'modal')
	set(hObject,'string',num2str(h.no_trace))
	return
elseif (get(handles.radiobutton_ap,'Value')==1) && no_trace > ap.data.ntrace
    errordlg([h.str.s55,num2str(ap.data.ntrace)],h.str.s45,'modal')
	set(hObject,'string',num2str(h.no_trace))
	return
end
h.no_trace = no_trace;
setappdata(handles.fig_bh_tomo_tt, 'h', h)
update_tout(hObject, eventdata, handles)

%
% -------------------------------------------------------------------------
%
function update_tout(hObject, eventdata, handles)
h = getappdata(handles.fig_bh_tomo_tt, 'h');
set(handles.trace_traite,'String',num2str(h.no_trace))
update_trace_simple(hObject, eventdata, handles)
update_traces_contigues(hObject, eventdata, handles)
update_positions_info(hObject, eventdata, handles)

%
% -------------------------------------------------------------------------
%
function update_trace_simple(hObject, eventdata, handles)
h = getappdata(handles.fig_bh_tomo_tt, 'h');
mog = getappdata(handles.fig_bh_tomo_tt, 'mog');
dc = median(h.proc_data(1:50,h.no_trace));
if get(handles.checkbox_pickTx,'Value')==0
    if get(handles.filtre_ondelette,'Value')==1
        plot(handles.trace_simple, mog.data.timestp, detrend_rad(mog.data.rdata(:,h.no_trace))-dc, 'color',[0.8 0.8 0.8])
        hold(handles.trace_simple, 'on')
        plot(handles.trace_simple, mog.data.timestp, h.proc_data(:,h.no_trace)-dc)
    else
        plot(handles.trace_simple, mog.data.timestp, h.proc_data(:,h.no_trace)-dc)
        hold(handles.trace_simple, 'on')
    end
    if ~isempty(mog.data.tdata)
        m1 = max(abs(h.proc_data(:,h.no_trace)-dc));
        m2 = max(abs(mog.data.tdata(:,h.no_trace)));
        plot(handles.trace_simple, mog.data.timestp,mog.data.tdata(:,h.no_trace)*m1/m2,'color',[0 0.6 0],'LineStyle','-.')
        if get(handles.filtre_ondelette,'Value')==1
           % legend(handles.trace_simple,'detrended Rx','Rx','Tx')
        else
            %legend(handles.trace_simple,'Rx','Tx')
        end
        
    end

    grid(handles.trace_simple, 'on')
    xlabel(handles.trace_simple, [h.str.s22,' [',mog.data.tunits,']'])
    ylabel(handles.trace_simple, h.str.s21)
    update_t_lim(hObject, eventdata, handles)
    if get(handles.dynALim,'Value') == 1.0
        al = 1.05*max(abs(h.proc_data(:,h.no_trace)-dc));
        set(handles.a_min,'String',num2str(-al))
        set(handles.a_max,'String',num2str(al))
    end
    update_a_lim(hObject, eventdata, handles)

    ylim = get(handles.trace_simple,'YLim');
    if length(h.pick)>=h.no_trace
        plot(handles.trace_simple, [h.pick(h.no_trace) h.pick(h.no_trace)],ylim,'g');
    end
    if length(h.pick_et)>=h.no_trace && h.pick_et(h.no_trace) ~= -1
        tmp = h.pick(h.no_trace)-h.pick_et(h.no_trace);
        plot(handles.trace_simple, [tmp tmp],ylim,'r');
        tmp = h.pick(h.no_trace)+h.pick_et(h.no_trace);
        plot(handles.trace_simple, [tmp tmp],ylim,'r');
    end
    
    if  ~isempty(h.pickTx)
        if length(h.pickTx)>=h.no_trace
            plot(handles.trace_simple, [h.pickTx(h.no_trace) h.pickTx(h.no_trace)],ylim,'k-');
        end
    end
else
    if get(handles.filtre_ondelette,'Value')==1
        plot(handles.trace_simple, mog.data.timestp, detrend_rad(mog.data.rdata(:,h.no_trace))-dc, 'color',[0.8 0.8 0.8],'LineStyle','-.')
        hold(handles.trace_simple, 'on')
        plot(handles.trace_simple, mog.data.timestp, h.proc_data(:,h.no_trace)-dc,'-.')
    else
        plot(handles.trace_simple, mog.data.timestp, h.proc_data(:,h.no_trace)-dc,'-.')
        hold(handles.trace_simple, 'on')
    end
    m1 = max(abs(h.proc_data(:,h.no_trace)-dc));
    m2 = max(abs(mog.data.tdata(:,h.no_trace)));
    plot(handles.trace_simple, mog.data.timestp,mog.data.tdata(:,h.no_trace)*m1/m2,'color',[0 0.6 0])
    if get(handles.filtre_ondelette,'Value')==1
       % legend(handles.trace_simple,'detrended Rx','Rx','Tx')
    else
       % legend(handles.trace_simple,'Rx','Tx')
	end

	grid(handles.trace_simple, 'on')
    xlabel(handles.trace_simple, [h.str.s22,' [',mog.data.tunits,']'])
    ylabel(handles.trace_simple, h.str.s21)
    update_t_lim(hObject, eventdata, handles)
    if get(handles.dynALim,'Value') == 1.0
        al = 1.05*max(abs(h.proc_data(:,h.no_trace)-dc));
        set(handles.a_min,'String',num2str(-al))
        set(handles.a_max,'String',num2str(al))
    end
    update_a_lim(hObject, eventdata, handles)

end
hold(handles.trace_simple, 'off')
%
% -------------------------------------------------------------------------
%
function update_traces_contigues(hObject, eventdata, handles)
h = getappdata(handles.fig_bh_tomo_tt, 'h');
mog = getappdata(handles.fig_bh_tomo_tt, 'mog');
av = getappdata(handles.fig_bh_tomo_tt, 'av');
ap = getappdata(handles.fig_bh_tomo_tt, 'ap');
%axes(handles.traces_contigues)
%b1=hsv(64);
%b2=gray(64);
%colormap([b2(1:32,:);b1;gray(32);b1;b2(33:64,:)])
colormap ramac_cmap(16)
tmin_axe=str2double(get(handles.t_min,'String'));
tmax_axe=str2double(get(handles.t_max,'String'));
if (get(handles.radiobutton_data,'Value')==1)
  %    ind = (mog.data.Tx_z(h.no_trace) == mog.data.Tx_z);
  if ~h.data_vrp
    ind = false(size(mog.data.Tx_z));
	ind(h.no_trace)=true;
	nn=h.no_trace+1;
	while nn<=mog.data.ntrace && ...
			  mog.data.Tx_z(nn) == mog.data.Tx_z(h.no_trace)
	  ind(nn)=true;
	  nn=nn+1;
	end
	nn=h.no_trace-1;
	while nn>0 && mog.data.Tx_z(nn) == mog.data.Tx_z(h.no_trace)
	  ind(nn)=true;
	  nn=nn-1;
	end
  else
	ind = false(size(mog.data.Tx_z));
	ind(h.no_trace)=true;
	nn=h.no_trace+1;
	while  nn<=mog.data.ntrace && ...
			   mog.data.Tx_x(nn) == mog.data.Tx_x(h.no_trace)
	  ind(nn)=true;
	  nn=nn+1;
	end
	nn=h.no_trace-1;
	while nn>0 && mog.data.Tx_x(nn) == mog.data.Tx_x(h.no_trace)
	  ind(nn)=true;
	  nn=nn-1;
	end
  end
%  ind = ind & mog.in;
  no_traces = 1:mog.data.ntrace;
  ind2 = ind&h.done;
  out = ~mog.in;
elseif (get(handles.radiobutton_av,'Value')==1)
  ind = 1:av.data.ntrace;
  no_traces = ind;
  ind2 = av.tt_done;
  out = false(size(ind));
elseif (get(handles.radiobutton_ap,'Value')==1)
  ind = 1:ap.data.ntrace;
  no_traces = ind;
  ind2 = ap.tt_done;
  out = false(size(ind));
end

imagesc(no_traces(ind), mog.data.timestp, h.proc_data(:,ind),'Parent',handles.traces_contigues)
caxis(handles.traces_contigues, h.cminmax)
ylabel(handles.traces_contigues, [h.str.s22,' [',mog.data.tunits,']'])
xlabel(handles.traces_contigues, h.str.s24)
update_t_lim(hObject, eventdata, handles)
hold(handles.traces_contigues, 'on')
vapp = calculVapp(handles);
if vapp~=0 && get(handles.radiobutton_data,'Value')==1 && ...
		get(handles.montre_vapp,'Value')==1
  l = sqrt((mog.data.Tx_x(ind)-mog.data.Rx_x(ind)).^2+...
		   (mog.data.Tx_y(ind)-mog.data.Rx_y(ind)).^2+...
		   (mog.data.Tx_z(ind)-mog.data.Rx_z(ind)).^2);
  tvapp=l./vapp;
  plot(handles.traces_contigues, no_traces(ind),tvapp,'y:')
end
%tmp = (tmax_axe-1)*ones(size(no_traces(ind2)));
if any(ind & out)
	plot(handles.traces_contigues, no_traces(ind & out), tmax_axe-1,'ro','MarkerFaceColor','r');
end
plot(handles.traces_contigues, h.no_trace,tmin_axe+1,'ys','MarkerFaceColor','y')
plot(handles.traces_contigues, no_traces(ind2),tmax_axe-ones(size(no_traces(ind2))),'gs','MarkerFaceColor','g')
plot(handles.traces_contigues, no_traces(ind),h.pick(ind),'go')
ind = ind & (h.pick_et ~= -1);
plot(handles.traces_contigues, no_traces(ind),h.pick(ind)-h.pick_et(ind),'ro')
plot(handles.traces_contigues, no_traces(ind),h.pick(ind)+h.pick_et(ind),'ro')

hold(handles.traces_contigues, 'off')

%
% -------------------------------------------------------------------------
%
function update_t_lim(hObject, eventdata, handles)
tlim = [str2double(get(handles.t_min,'String')) str2double(get(handles.t_max,'String')) ];
set(handles.trace_simple,'XLim',tlim)
set(handles.traces_contigues, 'YLim', tlim)

%
% -------------------------------------------------------------------------
%
function update_a_lim(hObject, eventdata, handles)
alim = [str2double(get(handles.a_min,'String')) str2double(get(handles.a_max,'String')) ];
set(handles.trace_simple,'YLim',alim)


%
% -------------------------------------------------------------------------
%
function t_min_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

%
% -------------------------------------------------------------------------
%
function t_min_Callback(hObject, eventdata, handles)
update_t_lim(hObject, eventdata, handles)

%
% -------------------------------------------------------------------------
%
function t_max_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


%
% -------------------------------------------------------------------------
%
function t_max_Callback(hObject, eventdata, handles)
update_t_lim(hObject, eventdata, handles)

%
% -------------------------------------------------------------------------
%
function a_min_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

%
% -------------------------------------------------------------------------
%
function a_min_Callback(hObject, eventdata, handles)
update_a_lim(hObject, eventdata, handles)

%
% -------------------------------------------------------------------------
%
function a_max_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

%
% -------------------------------------------------------------------------
%
function a_max_Callback(hObject, eventdata, handles)
update_a_lim(hObject, eventdata, handles)


%
% -------------------------------------------------------------------------
%
function a_zoom_CreateFcn(hObject, eventdata, handles)
usewhitebg = 1;
if usewhitebg
    set(hObject,'BackgroundColor',[.9 .9 .9]);
else
    set(hObject,'BackgroundColor',get(0,'defaultUicontrolBackgroundColor'));
end


%
% -------------------------------------------------------------------------
%
function a_zoom_Callback(hObject, eventdata, handles)
val = 10^get(hObject,'Value');
set(handles.a_min,'String',-val);
set(handles.a_max,'String', val);
update_trace_simple(hObject, eventdata, handles)

%
% -------------------------------------------------------------------------
%
function update_positions_info(hObject, eventdata, handles)
h = getappdata(handles.fig_bh_tomo_tt, 'h');
mog = getappdata(handles.fig_bh_tomo_tt, 'mog');
av = getappdata(handles.fig_bh_tomo_tt, 'av');
ap = getappdata(handles.fig_bh_tomo_tt, 'ap');
if (get(handles.radiobutton_data,'Value')==1)
    Tx_x = mog.data.Tx_x(h.no_trace);
    Tx_y = mog.data.Tx_y(h.no_trace);
    Tx_z = mog.data.Tx_z(h.no_trace);
    Rx_x = mog.data.Rx_x(h.no_trace);
    Rx_y = mog.data.Rx_y(h.no_trace);
    Rx_z = mog.data.Rx_z(h.no_trace);
    Tx_d = mog.Tx_z_orig(h.no_trace);
	in = mog.in;
elseif (get(handles.radiobutton_av,'Value')==1)
    Tx_x = av.data.Tx_x(h.no_trace);
    Tx_y = av.data.Tx_y(h.no_trace);
    Tx_z = av.data.Tx_z(h.no_trace);
    Rx_x = av.data.Rx_x(h.no_trace);
    Rx_y = av.data.Rx_y(h.no_trace);
    Rx_z = av.data.Rx_z(h.no_trace);
    Tx_d = av.d_TxRx(h.no_trace);
	in = av.in;
elseif (get(handles.radiobutton_ap,'Value')==1)
    Tx_x = ap.data.Tx_x(h.no_trace);
    Tx_y = ap.data.Tx_y(h.no_trace);
    Tx_z = ap.data.Tx_z(h.no_trace);
    Rx_x = ap.data.Rx_x(h.no_trace);
    Rx_y = ap.data.Rx_y(h.no_trace);
    Rx_z = ap.data.Rx_z(h.no_trace);
    Tx_d = ap.d_TxRx(h.no_trace);
	in = ap.in;
end
texte{1} = '';
texte{2} = 'Positions Tx -- Rx';
texte{3} = '';
texte{4} = '       x        y        z';
texte{5} = ['Tx:   ',num2str(Tx_x,'%5.2f'),'   ',num2str(Tx_y, '%5.2f'),'   ',num2str(Tx_z,'%5.2f'),' (',num2str(Tx_d),')'];
texte{6} = ['Rx:   ',num2str(Rx_x,'%5.2f'),'   ',num2str(Rx_y, '%5.2f'),'   ',num2str(Rx_z,'%5.2f')];
texte{7} = [num2str(sum(in)), ' traces, (', num2str(size(h.proc_data,2)),' total)'];
texte{8} = '';
texte{9} = [num2str(round(100*sum(h.done(in))/sum(in))),' % done'];
set(handles.positionsTxRx,'String',texte)

%
% -------------------------------------------------------------------------
%
function pointe(hObject, eventdata, handles)

if get(handles.checkbox_pickTx,'Value') == 1
    pointe_simple(hObject, eventdata, handles)
elseif get(handles.action,'Value') == 1
  pointe_ecart_type(hObject, eventdata, handles)
elseif get(handles.action,'Value') == 2
  pointe_simple(hObject, eventdata, handles)
end
update_tout(hObject, eventdata, handles)
set(handles.pushbutton_active_pointe,'Value',0)

if get(handles.checkbox_interm_save,'Value')==1
	sauve_mat(handles,true)
end

%
% -------------------------------------------------------------------------
%
function pointe_contigues(hObject, eventdata, handles)

if get(handles.action,'Value') == 1
	if get(handles.radiobutton_tt,'Value') == 1
		pointe_contigues_simple(hObject, eventdata, handles)
    elseif get(handles.radiobutton_et,'Value') == 1
		pointe_contigues_ecart_type(hObject, eventdata, handles)
    else
        selectionne_contigues(hObject, eventdata, handles)
	end
elseif get(handles.action,'Value') == 2
  pointe_contigues_simple(hObject, eventdata, handles)
end
update_tout(hObject, eventdata, handles)
set(handles.pushbutton_active_pointe,'Value',0)

if get(handles.checkbox_interm_save,'Value')==1, sauve_mat(handles,true), end

%
% -------------------------------------------------------------------------
%
function pointe_ecart_type(hObject, eventdata, handles)
h = getappdata(handles.fig_bh_tomo_tt, 'h');

% button = 1;
% tcmp = 0;
% acmp = 1;
pointe_actif=1;
xlim = get(handles.trace_simple,'XLim');
ylim = get(handles.trace_simple,'YLim');
n_trace_traite = 0;
tt = h.pick(h.no_trace);
et = h.pick_et(h.no_trace);
while pointe_actif
	update_tout(hObject, eventdata, handles)
	[t1,a,button] = ginput(1);
	pushOnce = true;
	while button ~= 2
		if gca ~= handles.trace_simple || t1<xlim(1) || t1>xlim(2) || a<ylim(1) || a>ylim(2)
			%pointe_actif=0;
			setappdata(handles.fig_bh_tomo_tt, 'h', h)
			return
		end
		if button==1
			tt = t1;
			update_images(hObject, eventdata, handles, tt, et)
		elseif button==3
			et = abs(t1-tt);
			update_images(hObject, eventdata, handles, tt, et)
		end
		[t1,a,button] = ginput(1);
		pushOnce=false;
	end
	% button == 2
	if pushOnce
		tt=-1;
	end
	h.done(h.no_trace) = true;
	h.pick(h.no_trace) = tt;
	h.pick_et(h.no_trace) = et;
	
	n_trace_traite = n_trace_traite+1;
	if rem(n_trace_traite,50)==0 && get(handles.checkbox_interm_save,'Value')==1
		sauve_mat(handles,true)
	end
	
    if get(handles.checkbox_jump,'Value')==1
        next = prochaine_trace(handles, h.no_trace, h.no_trace+1);
    else
        next = h.no_trace+1;
    end
	if  next == h.no_trace
		pointe_actif=0;
	else
		h.no_trace = next;
	end	
% 	h.no_trace = h.no_trace+1;
% 	if h.no_trace<size(h.proc_data,2) && sum(h.done)~=length(h.done)
% 	  while h.done(h.no_trace)
% 		h.no_trace = h.no_trace+1;
% 		if h.no_trace>=length(h.done)
% 		  break
% 		end
% 	  end
% 	end
% 	if h.no_trace > size(h.proc_data,2)
% 	  h.no_trace = h.no_trace-1;
% 	  pointe_actif=0;
% 	end
	set(handles.trace_traite,'String',num2str(h.no_trace))
	setappdata(handles.fig_bh_tomo_tt, 'h', h)
	tt = h.pick(h.no_trace);
	et = h.pick_et(h.no_trace);
end

		
%
% -------------------------------------------------------------------------
%
function pointe_simple(hObject, eventdata, handles)
h = getappdata(handles.fig_bh_tomo_tt, 'h');

%button = 1;
%tcmp = 0;
%acmp = 1;

if get(handles.checkbox_pickTx,'Value') == 1
	pick = h.pickTx;
	done = h.doneTx;
else
	pick = h.pick;
	done = h.done;
end

pointe_actif=1;
xlim = get(handles.trace_simple,'XLim');
ylim = get(handles.trace_simple,'YLim');
n_trace_traite = 0;
tt = pick(h.no_trace);
while pointe_actif
	update_tout(hObject, eventdata, handles)
	[t1,a,button] = ginput(1);
	pushOnce=true;
	while button ~= 2
		if gca ~= handles.trace_simple || t1<xlim(1) || t1>xlim(2) || a<ylim(1) || a>ylim(2)
			%pointe_actif=0;
			setappdata(handles.fig_bh_tomo_tt, 'h', h)
			return
		end
		if button==1
			tt = t1;
			update_images(hObject, eventdata, handles, tt, -1)
		end
		[t1,a,button] = ginput(1);
		pushOnce=false;
	end
	if pushOnce
		tt=-1;
	end
	% button == 2
	done(h.no_trace) = true;
	pick(h.no_trace) = tt;
	
	n_trace_traite = n_trace_traite+1;
	if rem(n_trace_traite,50)==0 && get(handles.checkbox_interm_save,'Value')==1
		sauve_mat(handles,true)
    end
	
    if get(handles.checkbox_jump,'Value')==1
        next = prochaine_trace(handles, h.no_trace, h.no_trace+1);
    else
        next = h.no_trace+1;
    end
	if  next == h.no_trace
		pointe_actif=0;
	else
		h.no_trace = next;
	end	

% 	h.no_trace = h.no_trace+1;
% 	if h.no_trace<size(h.proc_data,2)
% 	  while h.done(h.no_trace)
% 		h.no_trace = h.no_trace+1;
% 		if h.no_trace>=length(h.done)
% 		  break
% 		end
% 	  end
% 	end
% 	if h.no_trace > size(h.proc_data,2)
% 	  h.no_trace = h.no_trace-1;
% 	  pointe_actif=0;
% 	end
	set(handles.trace_traite,'String',num2str(h.no_trace))
	if get(handles.checkbox_pickTx,'Value') == 1
		h.pickTx = pick;
		h.doneTx = done;
	else
		h.pick = pick;
		h.done = done;
	end
	setappdata(handles.fig_bh_tomo_tt, 'h', h)
	tt = pick(h.no_trace);
end


%
% -------------------------------------------------------------------------
%
function pushbutton_active_pointe_Callback(hObject, eventdata, handles)
if get(hObject,'Value') == 1
    axes(handles.trace_simple)
	pointe(hObject, eventdata, handles)
end

%
% -------------------------------------------------------------------------
%
function pushbutton_active_pointe_contigues_Callback(hObject, eventdata, handles)
if get(hObject,'Value') == 1
    axes(handles.traces_contigues)
	pointe_contigues(hObject, eventdata, handles)
end


%
% -------------------------------------------------------------------------
%
function sauve_mat(handles, msg)
hh=0;
if msg, hh=msgbox('Sauvegarde intermediaire'); end
replace_data(handles)
h = getappdata(handles.fig_bh_tomo_tt, 'h');
mog = getappdata(handles.fig_bh_tomo_tt, 'mog');
av = getappdata(handles.fig_bh_tomo_tt, 'av');
ap = getappdata(handles.fig_bh_tomo_tt, 'ap');

if strcmp(mog.data.tunits,'ns')==1
    if mog.user_fac_dt == 0
        if mog.data.synthetique==1
            fac_dt_av = 1;
            fac_dt_ap = 1;
        else
            [~, fac_dt_av, fac_dt_ap] = corr_t0(mog.data.ntrace, av, ap, false);
        end
        
        if ~isempty(av), av.fac_dt = fac_dt_av; end
        if ~isempty(ap), ap.fac_dt = fac_dt_ap; end
        if fac_dt_av~=1 && fac_dt_ap ~= 1
            mog.fac_dt = 0.5*(fac_dt_av+fac_dt_ap);
        elseif fac_dt_av~=1
            mog.fac_dt = fac_dt_av;
        elseif fac_dt_ap~=1
            mog.fac_dt = fac_dt_ap;
        else
            mog.fac_dt = 1;
        end
    end
end

load(h.db_file,'mogs','air')
mogs(h.no_mog) = mog; %#ok<NASGU>
if ~isempty(av), air(mog.av) = av; end
if ~isempty(ap), air(mog.ap) = ap; end %#ok<NASGU>
save(h.db_file,'mogs','air','-append')
if msg && ishandle(hh)==1
	delete(hh)
end
h.saved = true;
setappdata(handles.fig_bh_tomo_tt, 'h', h)

%
% -------------------------------------------------------------------------
%
function update_images(hObject, eventdata, handles, tt, et)

h = getappdata(handles.fig_bh_tomo_tt, 'h');
update_trace_simple(hObject, eventdata, handles)
%axes(handles.trace_simple)
ylim = get(handles.trace_simple,'YLim');
hold(handles.trace_simple, 'on')
plot(handles.trace_simple, [tt tt],ylim,'g');
if et ~= -1
	plot(handles.trace_simple, [tt-et tt-et],ylim,'r');
	plot(handles.trace_simple, [tt+et tt+et],ylim,'r');
end
hold(handles.trace_simple, 'off')

update_traces_contigues(hObject, eventdata, handles)
%axes(handles.traces_contigues)
hold(handles.traces_contigues, 'on')
plot(handles.traces_contigues, h.no_trace,tt,'go');
if et ~= -1
	plot(handles.traces_contigues, h.no_trace,tt-et,'ro');
	plot(handles.traces_contigues, h.no_trace,tt+et,'ro');
end
hold(handles.traces_contigues, 'off')


%
% -------------------------------------------------------------------------
%
function pushbutton_trace_suivante_Callback(hObject, eventdata, handles)
h = getappdata(handles.fig_bh_tomo_tt, 'h');
if (get(handles.radiobutton_data,'Value')==1)
	mog = getappdata(handles.fig_bh_tomo_tt, 'mog');
	in = mog.in;
elseif (get(handles.radiobutton_av,'Value')==1)
	av = getappdata(handles.fig_bh_tomo_tt, 'av');
	in = av.in;
elseif (get(handles.radiobutton_ap,'Value')==1)
	ap = getappdata(handles.fig_bh_tomo_tt, 'ap');
	in = ap.in;
end
nos = 1:size(h.proc_data,2);
nos = nos(in);
ind = find(h.no_trace == nos);
if ind < length(nos)
    h.no_trace = nos(ind+1);
    setappdata(handles.fig_bh_tomo_tt, 'h', h)
    update_tout(hObject, eventdata, handles)
end

%
% -------------------------------------------------------------------------
%
function pushbutton_trace_precedente_Callback(hObject, eventdata, handles)
h = getappdata(handles.fig_bh_tomo_tt, 'h');
if (get(handles.radiobutton_data,'Value')==1)
	mog = getappdata(handles.fig_bh_tomo_tt, 'mog');
	in = mog.in;
elseif (get(handles.radiobutton_av,'Value')==1)
	av = getappdata(handles.fig_bh_tomo_tt, 'av');
	in = av.in;
elseif (get(handles.radiobutton_ap,'Value')==1)
	ap = getappdata(handles.fig_bh_tomo_tt, 'ap');
	in = ap.in;
end
nos = 1:size(h.proc_data,2);
nos = nos(in);
ind = find(h.no_trace == nos);
if ind > 1
    h.no_trace = nos(ind-1);
    setappdata(handles.fig_bh_tomo_tt, 'h', h)
    update_tout(hObject, eventdata, handles)
end


%
% -------------------------------------------------------------------------
%
function pushbutton_prochaine_trace_Callback(hObject, eventdata, handles)
h = getappdata(handles.fig_bh_tomo_tt, 'h');
next = prochaine_trace(handles, h.no_trace,1);
if next ~= h.no_trace
	h.no_trace = next;
	setappdata(handles.fig_bh_tomo_tt, 'h', h)
	update_tout(hObject, eventdata, handles)
end


%
% -------------------------------------------------------------------------
%
function next = prochaine_trace(handles, no_trace, debut)
next = no_trace;
h = getappdata(handles.fig_bh_tomo_tt, 'h');
if (get(handles.radiobutton_data,'Value')==1)
	mog = getappdata(handles.fig_bh_tomo_tt, 'mog');
	in = mog.in;
elseif (get(handles.radiobutton_av,'Value')==1)
	av = getappdata(handles.fig_bh_tomo_tt, 'av');
	in = av.in;
elseif (get(handles.radiobutton_ap,'Value')==1)
	ap = getappdata(handles.fig_bh_tomo_tt, 'ap');
	in = ap.in;
end
for n=debut:size(h.proc_data,2)
	if ~h.done(n) && in(n)
		next = n;
		break
	end
end



%
% -------------------------------------------------------------------------
%
function SaveMenuItem_Callback(hObject, eventdata, handles)
sauve_mat(handles,false)


%
% -------------------------------------------------------------------------
%
function replace_data(handles)

h = getappdata(handles.fig_bh_tomo_tt, 'h');
mog = getappdata(handles.fig_bh_tomo_tt, 'mog');
av = getappdata(handles.fig_bh_tomo_tt, 'av');
ap = getappdata(handles.fig_bh_tomo_tt, 'ap');
if get(handles.radiobutton_data,'Value')==1
	mog.tt = h.pick;
	mog.et = h.pick_et;
	mog.tt_done = h.done;
elseif get(handles.radiobutton_av,'Value')==1
	av.tt = h.pick;
	av.et = h.pick_et;
	av.tt_done = h.done;
elseif get(handles.radiobutton_ap,'Value')==1
	ap.tt = h.pick;
	ap.et = h.pick_et;
	ap.tt_done = h.done;
end
if strcmp(get(handles.checkbox_pickTx,'Enable'),'on')
    mog.ttTx = h.pickTx;
    mog.ttTx_done = h.doneTx;
end
setappdata(handles.fig_bh_tomo_tt, 'h', h)
setappdata(handles.fig_bh_tomo_tt, 'mog', mog)
setappdata(handles.fig_bh_tomo_tt, 'av', av)
setappdata(handles.fig_bh_tomo_tt, 'ap', ap)

%
% -------------------------------------------------------------------------
%
function HelpMenuItem_Callback(hObject, eventdata, handles)
h = getappdata(handles.fig_bh_tomo_tt, 'h');
helpdlg(h.str.s61, h.str.s32)


%
% -------------------------------------------------------------------------
%
function filtre_median_Callback(hObject, eventdata, handles)
update_trace_simple(hObject, eventdata, handles)


%
% -------------------------------------------------------------------------
%
function pushbutton_plot_stats_Callback(hObject, eventdata, handles)
h = getappdata(handles.fig_bh_tomo_tt, 'h');
mog = getappdata(handles.fig_bh_tomo_tt, 'mog');

if strcmp(mog.data.tunits,'ns')==1
	if ~isempty(mog.ttTx)
		t0 = mog.ttTx;
	elseif ~(h.init_av & h.init_ap) %#ok<AND2>
        warndlg(h.str.s59, h.str.s53);
        t0 = zeros(1,mog.data.ntrace);
        %t0av = 0;
        %t0ap = 0;
    else
        %[t0, t0av, t0ap] = get_t0(hObject, eventdata, handles);
        %t0 = get_t0(hObject, eventdata, handles);
        av = getappdata(handles.fig_bh_tomo_tt, 'av');
        ap = getappdata(handles.fig_bh_tomo_tt, 'ap');
        if mog.data.synthetique==1
            t0 = zeros(1,mog.data.ntrace);
            fac_dt_av = 1;
            fac_dt_ap = 1;
        else
            [t0, fac_dt_av, fac_dt_ap] = corr_t0(mog.data.ntrace, av, ap, true);
        end
        
        if ~isempty(av), av.fac_dt = fac_dt_av; end
        if ~isempty(ap), ap.fac_dt = fac_dt_ap; end
        if mog.user_fac_dt == 0
            if fac_dt_av~=1 && fac_dt_ap ~= 1
                mog.fac_dt = 0.5*(fac_dt_av+fac_dt_ap);
            elseif fac_dt_av~=1
                mog.fac_dt = fac_dt_av;
            elseif fac_dt_ap~=1
                mog.fac_dt = fac_dt_ap;
            else
                mog.fac_dt = 1;
            end
        elseif fac_dt_av~=1 || fac_dt_ap ~= 1
            disp('dt correction factor NOT edited at user request (set in bh_tomo_db)')
        end
        setappdata(handles.fig_bh_tomo_tt, 'av', av)
        setappdata(handles.fig_bh_tomo_tt, 'ap', ap)
        setappdata(handles.fig_bh_tomo_tt,'mog', mog);
    end
end

ind2 = (mog.tt ~= -1);
ind = mog.tt_done & ind2 & mog.in;

if get(handles.action,'Value')==1  % pointe avec ecart-type
  tt = mog.tt(ind);
  et = mog.et(ind);
  if sum(et==-1) > 0
	no = 1:length(ind);
	disp([h.str.s60, ': traces'])
    no = no(ind&(mog.et==-1));
	disp(no)
	warndlg(h.str.s60, h.str.s53);
  end
else
  tt = mog.tt(ind);
  et = zeros(size(tt));
end

tt = tt*mog.fac_dt;
et = et*mog.fac_dt;
t0 = t0*mog.fac_dt;

hyp = sqrt((mog.data.Tx_x(ind)-mog.data.Rx_x(ind)).^2+...
		   (mog.data.Tx_y(ind)-mog.data.Rx_y(ind)).^2+...
		   (mog.data.Tx_z(ind)-mog.data.Rx_z(ind)).^2);
dz = mog.data.Rx_z(ind)-mog.data.Tx_z(ind);
theta = 180/pi*asin(dz./hyp);

h_stat = figure;
set(h_stat, 'Position',[256 71 512 620]);

subplot(321)
plot(hyp, tt, 'o')
xlabel(h.str.s35)
ylabel([h.str.s22,' [',mog.data.tunits,']'])

subplot(322)
plot(theta, hyp./tt, 'o')
xlabel(h.str.s36)
ylabel([h.str.s37, '[',mog.data.cunits,'/',mog.data.tunits,']'])
title(h.str.s38)

vapp = hyp./(tt-t0(ind));
no = 1:length(ind);
no=no(ind);
ind2 = vapp<0;
if ~isempty(no(ind2))
  disp(h.str.s40)
  disp(no(ind2))
end

ind2 = vapp>0.3;
if sum(ind2)>0
    disp('Apparent velocity higher that 0.3: traces no')
    no2 = 1:length(ind2);
    no2 = no2(ind2);
	disp([no2' vapp(ind2)'])
end

subplot(323)
plot(theta, hyp./(tt-t0(ind)), 'o')
xlabel(h.str.s36)
ylabel([h.str.s37, '[',mog.data.cunits,'/',mog.data.tunits,']'])
title(h.str.s39)

subplot(324)
plot(t0)
xlabel(h.str.s41)
ylabel([h.str.s22,' [',mog.data.tunits,']'])
title(h.str.s42)

subplot(325)
plot(hyp,et,'o')
xlabel(h.str.s35)
ylabel(h.str.s43)

subplot(326)
plot(theta,et,'o')
xlabel(h.str.s36)
ylabel(h.str.s43)

h_stat = figure;
set(h_stat, 'Position',[256 71 512 620]);
rmin = min(vapp);
if rmin<0
    rmin = 1.001*rmin;
else
    rmin = 0.999*rmin;
end
rmax = 1.001*max(vapp);
c=jet;%cmr;
m = (size(c,1)-1)/(rmax-rmin);
b = 1-rmin*m;
p = m*vapp(1)+b;
couleur = interp1(c,p);

n=1;
plot([mog.data.Tx_x(n) mog.data.Rx_x(n)],...
    [mog.data.Tx_z(n) mog.data.Rx_z(n)], 'Color',couleur)
hold on    
for n=2:length(no)
    p = m*vapp(n)+b;
    couleur = interp1(c,p);
    plot([mog.data.Tx_x(n) mog.data.Rx_x(n)],...
        [mog.data.Tx_z(n) mog.data.Rx_z(n)], 'Color',couleur)
end
hold off
set(gca,'DataAspectRatio',[1 1 1],'YDir','normal')
axis tight
caxis([rmin rmax])
hb=colorbar(jet);
set(get(hb,'Title'),'String','App. Vel. (m/ns)')
xlabel(h.str.s119)
ylabel(h.str.s120)


[hyp, iHyp] = sort(hyp);
tt = tt(iHyp);
no = no(iHyp);

npts = round(length(hyp)/20);
out = [];
for n=1:19
	t=tt((1+(n-1)*npts):n*npts);
	tm = mean(t);
	ts = std(t);
	i = t<tm-3*ts | t>tm+3*ts;
	if any(i)
		nn = no((1+(n-1)*npts):n*npts);
		out = [out nn(i)]; %#ok<AGROW>
	end
end
t=tt((1+19*npts):end);
tm = mean(t);
ts = std(t);
i = t<tm-2*ts | t>tm+2*ts;
if any(i)
	nn = no((1+19*npts):end);
	out = [out nn(i)];
end
if ~isempty(out)
	out = sort(out)';
	disp('outliers - 3 sigma')
	disp(out)
end

%
% -------------------------------------------------------------------------
%
function filtre_ondelette_Callback(hObject, eventdata, handles)
if get(handles.radiobutton_data,'Value') == 0
	return
end
h = getappdata(handles.fig_bh_tomo_tt, 'h');
mog = getappdata(handles.fig_bh_tomo_tt, 'mog');
if ~h.init_data
    warndlg(h.str.s52, h.str.s53)
    set(hObject,'Value',0);
    return
end
if get(hObject,'Value')==0
	if get(handles.radiobutton_data,'Value')==1
		h.proc_data = detrend_rad(mog.data.rdata);
	end
else
  if isempty( mog.fw )
    mog.fw = filtrage_wavelet2(mog.data.rdata);
  end
  h.proc_data = detrend_rad(mog.fw);
end
cmin=-max(max(abs(h.proc_data)));
cmax= -cmin;
h.cminmax = [cmin cmax];

setappdata(handles.fig_bh_tomo_tt,'h',h)
setappdata(handles.fig_bh_tomo_tt,'mog',mog)
update_tout(hObject, eventdata, handles)


%
% -------------------------------------------------------------------------
%
function boutons_trace_contigues_SelectionChangeFcn(hObject, eventdata, handles)


%
% -------------------------------------------------------------------------
%
function radiobutton_av_Callback(hObject, eventdata, handles)
% if get(handles.radiobutton_data,'Value')==0 & ...
% 		get(handles.radiobutton_ap,'Value')==0
% 	return
% end
h = getappdata(handles.fig_bh_tomo_tt, 'h');
mog = getappdata(handles.fig_bh_tomo_tt, 'mog');
av = getappdata(handles.fig_bh_tomo_tt, 'av');
ap = getappdata(handles.fig_bh_tomo_tt, 'ap');
if h.init_av == 0
     return
end

if get(handles.radiobutton_data,'Value')==1
	mog.tt = h.pick;
	mog.et = h.pick_et;
	mog.tt_done = h.done;
elseif get(handles.radiobutton_ap,'Value')==1
	ap.tt = h.pick;
	ap.et = h.pick_et;
	ap.tt_done = h.done;
end
h.pick = av.tt;
h.pick_et = av.et;
h.done = av.tt_done;

h.no_trace = 1;
set(handles.trace_traite,'String',num2str(h.no_trace))

off = [handles.radiobutton_data,handles.radiobutton_ap];
mutual_exclude(off)
h.proc_data = detrend_rad(av.data.rdata);
cmin=-max(max(abs(h.proc_data)));
cmax= -cmin;
h.cminmax = [cmin cmax];
setappdata(handles.fig_bh_tomo_tt, 'h', h)
setappdata(handles.fig_bh_tomo_tt, 'mog', mog)
setappdata(handles.fig_bh_tomo_tt, 'ap', ap)
update_tout(hObject, eventdata, handles)


%
% -------------------------------------------------------------------------
%
function radiobutton_ap_Callback(hObject, eventdata, handles)
% if get(handles.radiobutton_data,'Value')==0 & ...
% 		get(handles.radiobutton_av,'Value')==0
% 	return
% end
h = getappdata(handles.fig_bh_tomo_tt, 'h');
mog = getappdata(handles.fig_bh_tomo_tt, 'mog');
av = getappdata(handles.fig_bh_tomo_tt, 'av');
ap = getappdata(handles.fig_bh_tomo_tt, 'ap');
if h.init_ap == 0
    return
end

if get(handles.radiobutton_data,'Value')==1
	mog.tt = h.pick;
	mog.et = h.pick_et;
	mog.tt_done = h.done;
elseif get(handles.radiobutton_av,'Value')==1
	av.tt = h.pick;
	av.et = h.pick_et;
	av.tt_done = h.done;
end
h.pick = ap.tt;
h.pick_et = ap.et;
h.done = ap.tt_done;

h.no_trace = 1;
set(handles.trace_traite,'String',num2str(h.no_trace))

off = [handles.radiobutton_data,handles.radiobutton_av];
mutual_exclude(off)
h.proc_data = detrend_rad(ap.data.rdata);
cmin=-max(max(abs(h.proc_data)));
cmax= -cmin;
h.cminmax = [cmin cmax];
setappdata(handles.fig_bh_tomo_tt, 'h', h)
setappdata(handles.fig_bh_tomo_tt, 'av', av)
setappdata(handles.fig_bh_tomo_tt, 'mog', mog)
update_tout(hObject, eventdata, handles)


%
% -------------------------------------------------------------------------
%
function radiobutton_data_Callback(hObject, eventdata, handles)
% if get(handles.radiobutton_av,'Value')==0 & ...
% 		get(handles.radiobutton_ap,'Value')==0
% 	return
% end
h = getappdata(handles.fig_bh_tomo_tt, 'h');
mog = getappdata(handles.fig_bh_tomo_tt, 'mog');
av = getappdata(handles.fig_bh_tomo_tt, 'av');
ap = getappdata(handles.fig_bh_tomo_tt, 'ap');
if h.init_data == 0
    return
end

if get(handles.radiobutton_av,'Value')==1
	av.tt = h.pick;
	av.et = h.pick_et;
	av.tt_done = h.done;
elseif get(handles.radiobutton_ap,'Value')==1
	ap.tt = h.pick;
	ap.et = h.pick_et;
	ap.tt_done = h.done;
end
h.pick = mog.tt;
h.pick_et = mog.et;
h.done = mog.tt_done;

nos = 1:mog.data.ntrace;
nos = nos(mog.in);
h.no_trace = nos(1);
set(handles.trace_traite,'String',num2str(h.no_trace))

off = [handles.radiobutton_av,handles.radiobutton_ap];
mutual_exclude(off)

h.proc_data = detrend_rad(mog.data.rdata);
cmin=-max(max(abs(h.proc_data)));
cmax= -cmin;
h.cminmax = [cmin cmax];
setappdata(handles.fig_bh_tomo_tt, 'h', h)
setappdata(handles.fig_bh_tomo_tt, 'av', av)
setappdata(handles.fig_bh_tomo_tt, 'ap', ap)
update_tout(hObject, eventdata, handles)


%
% -------------------------------------------------------------------------
%
function pointe_contigues_simple(hObject, eventdata, handles)
h = getappdata(handles.fig_bh_tomo_tt, 'h');
mog = getappdata(handles.fig_bh_tomo_tt, 'mog');

pointe_actif=1;
%xlim = get(handles.traces_contigues,'XLim');
%ylim = get(handles.traces_contigues,'YLim');
n_trace_traite = 0;
while pointe_actif
	update_tout(hObject, eventdata, handles)
	[nt,tt,button] = ginput(1);
	xlim = get(handles.traces_contigues,'XLim');
	ylim = get(handles.traces_contigues,'YLim');
	while button ~= 2
		if gca ~= handles.traces_contigues || nt<xlim(1) || nt>xlim(2) ...
				  || tt<ylim(1) || tt>ylim(2)
			%pointe_actif=0;
			setappdata(handles.fig_bh_tomo_tt, 'h', h)
			return
		end
		h.no_trace = round(nt);
        h.done(h.no_trace) = true;
		if button==1
			h.pick(h.no_trace) = tt;			
		elseif button==3% && mog.in(h.no_trace)
			h.pick(h.no_trace) = -1;
			h.pick_et(h.no_trace) = -1;
		end
		set(handles.trace_traite,'String',num2str(h.no_trace))
		setappdata(handles.fig_bh_tomo_tt, 'h', h)
		update_tout(hObject, eventdata, handles)
		n_trace_traite = n_trace_traite+1;
		if rem(n_trace_traite,50)==0 && get(handles.checkbox_interm_save,'Value')==1
			sauve_mat(handles,true)
		end
		[nt,tt,button] = ginput(1);
	end
	% button == 2
	if h.data_vrp == 0
		while  h.no_trace<=mog.data.ntrace && ...
				mog.data.Tx_z(round(nt)) == mog.data.Tx_z(h.no_trace)
			h.no_trace = h.no_trace+1;
		end
	else
		while  h.no_trace<=mog.data.ntrace && ...
				mog.data.Tx_x(round(nt)) == mog.data.Tx_x(h.no_trace)
			h.no_trace = h.no_trace+1;
		end
	end
	if h.no_trace > mog.data.ntrace
		h.no_trace = h.no_trace-1;
		%pointe_actif=0;
		setappdata(handles.fig_bh_tomo_tt, 'h', h)
		return
	end
		
	set(handles.trace_traite,'String',num2str(h.no_trace))
	setappdata(handles.fig_bh_tomo_tt, 'h', h)
end


%
% -------------------------------------------------------------------------
%
function pointe_contigues_ecart_type(hObject, eventdata, handles)
h = getappdata(handles.fig_bh_tomo_tt, 'h');
mog = getappdata(handles.fig_bh_tomo_tt, 'mog');

pointe_actif=1;
%xlim = get(handles.traces_contigues,'XLim');
%ylim = get(handles.traces_contigues,'YLim');
n_trace_traite = 0;
while pointe_actif==1
  update_tout(hObject, eventdata, handles)
  [nt,tt,button] = ginput(1);
  xlim = get(handles.traces_contigues,'XLim');
  ylim = get(handles.traces_contigues,'YLim');
  while button ~= 2
	if gca ~= handles.traces_contigues || nt<xlim(1) || nt>xlim(2) ...
			  || tt<ylim(1) || tt>ylim(2)
	  %pointe_actif=0;
	  setappdata(handles.fig_bh_tomo_tt, 'h', h)
	  return
	end
	h.no_trace = round(nt);
	h.done(h.no_trace) = true;
	if button==1
	  if h.pick(h.no_trace) == -1
		warndlg(h.str.s62, h.str.s53)
	  else
		h.pick_et(h.no_trace) = abs(h.pick(h.no_trace)-tt);
	  end
	elseif button==3 && mog.in(h.no_trace)
	  h.pick_et(h.no_trace) = -1;
	end
	set(handles.trace_traite,'String',num2str(h.no_trace))
	setappdata(handles.fig_bh_tomo_tt, 'h', h)
	update_tout(hObject, eventdata, handles)
	n_trace_traite = n_trace_traite+1;
	if rem(n_trace_traite,50)==0 && get(handles.checkbox_interm_save,'Value')==1
	  sauve_mat(handles,true)
	end
	[nt,tt,button] = ginput(1);
  end
  % button == 2
  if h.data_vrp == 0
	while  h.no_trace<=mog.data.ntrace && ...
		  mog.data.Tx_z(round(nt)) == mog.data.Tx_z(h.no_trace)
	  h.no_trace = h.no_trace+1;
	end
  else
	while  h.no_trace<=mog.data.ntrace && ...
		  mog.data.Tx_x(round(nt)) == mog.data.Tx_x(h.no_trace)
	  h.no_trace = h.no_trace+1;
	end
  end
  if h.no_trace > mog.data.ntrace
	h.no_trace = h.no_trace-1;
	%pointe_actif=0;
	setappdata(handles.fig_bh_tomo_tt, 'h', h)
	return
  end
  
  set(handles.trace_traite,'String',num2str(h.no_trace))
  setappdata(handles.fig_bh_tomo_tt, 'h', h)
end


%
% -------------------------------------------------------------------------
%
function selectionne_contigues(hObject, eventdata, handles)
h = getappdata(handles.fig_bh_tomo_tt, 'h');
mog = getappdata(handles.fig_bh_tomo_tt, 'mog');
pointe_actif=1;
%xlim = get(handles.traces_contigues,'XLim');
%ylim = get(handles.traces_contigues,'YLim');
while pointe_actif
	update_tout(hObject, eventdata, handles)
	[nt,tt,button] = ginput(1);
	xlim = get(handles.traces_contigues,'XLim');
	ylim = get(handles.traces_contigues,'YLim');
	while button == 1 %#ok<ALIGN>
		if gca ~= handles.traces_contigues || nt<xlim(1) || nt>xlim(2) ...
				  || tt<ylim(1) || tt>ylim(2)
			%pointe_actif=0;
			setappdata(handles.fig_bh_tomo_tt, 'h', h)
			return
		end
		h.no_trace = round(nt);
		set(handles.trace_traite,'String',num2str(h.no_trace))
		setappdata(handles.fig_bh_tomo_tt, 'h', h)
		update_tout(hObject, eventdata, handles)
		[nt,tt,button] = ginput(1);
    end
    if button == 3
        if h.data_vrp == 0
            while h.no_trace<mog.data.ntrace && ...
                    mog.data.Tx_z(round(nt)) == mog.data.Tx_z(h.no_trace)
                h.no_trace = h.no_trace+1;
            end
        else
            while h.no_trace<=mog.data.ntrace && ...
                    mog.data.Tx_x(round(nt)) == mog.data.Tx_x(h.no_trace)
                h.no_trace = h.no_trace+1;
            end
        end
    else
        if h.data_vrp == 0
            while  h.no_trace>1 && ...
                    mog.data.Tx_z(round(nt)) == mog.data.Tx_z(h.no_trace)
                h.no_trace = h.no_trace-1;
            end
        else
            while  h.no_trace>1 && ...
                    mog.data.Tx_x(round(nt)) == mog.data.Tx_x(h.no_trace)
                h.no_trace = h.no_trace-1;
            end
        end
    end
		
	set(handles.trace_traite,'String',num2str(h.no_trace))
	setappdata(handles.fig_bh_tomo_tt, 'h', h)
end



%
% -------------------------------------------------------------------------
%
function pushbutton_reinit_trace_Callback(hObject, eventdata, handles)
h = getappdata(handles.fig_bh_tomo_tt, 'h');
h.pick(h.no_trace) = -1;
h.pick_et(h.no_trace) = -1;
h.done(h.no_trace) = 0;
setappdata(handles.fig_bh_tomo_tt,'h',h)
update_tout(hObject, eventdata, handles)


%
% -------------------------------------------------------------------------
%
function vapp=calculVapp(handles)
%h = getappdata(handles.fig_bh_tomo_tt, 'h');
mog = getappdata(handles.fig_bh_tomo_tt, 'mog');

ind2 = (mog.tt ~= -1);
ind = mog.tt_done & ind2 & mog.in;
if sum(ind)==0
  vapp = 0;
  return
end
hyp = sqrt((mog.data.Tx_x(ind)-mog.data.Rx_x(ind)).^2+...
		   (mog.data.Tx_y(ind)-mog.data.Rx_y(ind)).^2+...
		   (mog.data.Tx_z(ind)-mog.data.Rx_z(ind)).^2);
if get(handles.action,'Value')==1  % pointe avec ecart-type
  tt = mog.tt(ind);
  et = mog.et(ind);
else
  tt = mog.tt(ind);
  et = zeros(size(tt));
end

vapp = hyp./tt;
if et==0
  vapp=mean(vapp);
else
  w=1./et;
  vapp = sum(vapp.*w)/sum(w);
end


%
% -------------------------------------------------------------------------
%
function montre_vapp_Callback(hObject, eventdata, handles)
update_traces_contigues(hObject, eventdata, handles)



%
% -------------------------------------------------------------------------
%
function set_String_locale(handles, str)

set(handles.radiobutton_data,        'String', str.s01);
set(handles.radiobutton_av,          'String', str.s02);
set(handles.radiobutton_ap,          'String', str.s03);
set(handles.action,                  'String', str.s04);
set(handles.pushbutton_active_pointe,           'String', str.s05);
set(handles.pushbutton_active_pointe_contigues, 'String', str.s06);
set(handles.boutons_trace_contigues, 'Title',  str.s07);
set(handles.radiobutton_tt,          'String', str.s08);
set(handles.radiobutton_et,          'String', str.s09);
set(handles.radiobutton_sl,          'String', str.s10);
set(handles.pushbutton_plot_stats,              'String', str.s11);
set(handles.uipanel1,                'Title',  str.s13);
set(handles.trace_traitee_label,     'String', str.s14);
set(handles.pushbutton_trace_precedente,        'String', str.s15);
set(handles.pushbutton_trace_suivante,          'String', str.s16);
set(handles.pushbutton_prochaine_trace,         'String', str.s17);
set(handles.pushbutton_reinit_trace,            'String', str.s18);
set(handles.filtre_ondelette,        'String', str.s19);
set(handles.montre_vapp,             'String', str.s20);
set(handles.FileMenu,                'Label',  str.s25);
set(handles.OpenMenuItem,            'Label',  str.s26);
set(handles.save,                    'Label',  str.s29);
set(handles.CloseMenuItem,           'Label',  str.s31);
set(handles.ImportMenuItem,          'Label',  [str.s215,' ...']);
set(handles.aide,                    'Label',  str.s32);
set(handles.checkbox_interm_save,    'String', str.s56);
%set(handles.,          'String', str.s);


function ImportMenuItem_Callback(hObject, eventdata, handles)
h = getappdata(handles.fig_bh_tomo_tt, 'h');

[file,rep] = uigetfile('*.tt',h.str.s216);
if isequal(file,0) || isequal(rep,0)
    return
end
new_tt = load([rep,file]);
if size(new_tt,1) > length(h.pick) || max(new_tt(:,1)) > length(h.pick)
	uiwait(errordlg(h.str.s217))
	return
end
if size(new_tt,2) == 2
	ind = new_tt(:,1);
	h.done(ind) = true;
	h.pick(ind) = new_tt(:,2);
elseif size(new_tt,2) == 3
	ind = new_tt(:,1);
	h.done(ind) = true;
	h.pick(ind) = new_tt(:,2);
	h.pick_et(ind) = new_tt(:,3);
else
	uiwait(errordlg(h.str.s217))
	return
end
setappdata(handles.fig_bh_tomo_tt, 'h', h)


function checkbox_interm_save_Callback(hObject, eventdata, handles)


function trace_traitee_label_CreateFcn(hObject, eventdata, handles)
set(hObject,'BackgroundColor',get(0,'defaultUicontrolBackgroundColor'));


% --- Executes on button press in pushbutton_plot_stats.
%function pushbutton9_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton_plot_stats (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
%display('ici')


% --- Executes on key press with focus on fig_bh_tomo_tt or any of its controls.
function fig_bh_tomo_tt_WindowKeyPressFcn(hObject, eventdata, handles)
% hObject    handle to fig_bh_tomo_tt (see GCBO)
% eventdata  structure with the following fields (see FIGURE)
%	Key: name of the key that was pressed, in lower case
%	Character: character interpretation of the key(s) that was pressed
%	Modifier: name(s) of the modifier key(s) (i.e., control, shift) pressed
% handles    structure with handles and user data (see GUIDATA)

if eventdata.Key == 's'
    trace_suivante_Callback(hObject, eventdata, handles)
elseif eventdata.Key == 'a'
    trace_precedente_Callback(hObject, eventdata, handles)
elseif eventdata.Key == 't'
    axes(handles.trace_simple)
	pointe(hObject, eventdata, handles)
end


% --- Executes on button press in dynALim.
function dynALim_Callback(hObject, eventdata, handles)
% hObject    handle to dynALim (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of dynALim


% --- Executes on button press in checkbox_pickTx.
function checkbox_pickTx_Callback(hObject, eventdata, handles)
update_tout(hObject, eventdata, handles)


% --- Executes on button press in checkbox_jump.
function checkbox_jump_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox_jump (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkbox_jump
