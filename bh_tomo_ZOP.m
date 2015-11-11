function varargout = bh_tomo_ZOP(varargin)
% BH_TOMO_ZOP M-file for bh_tomo_ZOP.fig
%      BH_TOMO_ZOP, by itself, creates a new BH_TOMO_ZOP or raises the existing
%      singleton*.
%
%      H = BH_TOMO_ZOP returns the handle to a new BH_TOMO_ZOP or the handle to
%      the existing singleton*.
%
%      BH_TOMO_ZOP('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in BH_TOMO_ZOP.M with the given input arguments.
%
%      BH_TOMO_ZOP('Property','Value',...) creates a new BH_TOMO_ZOP or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before bh_tomo_ZOP_OpeningFunction gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to bh_tomo_ZOP_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Copyright 2002-2003 The MathWorks, Inc.

% Edit the above text to modify the response to help bh_tomo_ZOP

% Last Modified by GUIDE v2.5 11-Aug-2009 16:35:25

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @bh_tomo_ZOP_OpeningFcn, ...
                   'gui_OutputFcn',  @bh_tomo_ZOP_OutputFcn, ...
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


% --- Executes just before bh_tomo_ZOP is made visible.
function bh_tomo_ZOP_OpeningFcn(hObject, eventdata, handles, varargin) %#ok<*INUSL>
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to bh_tomo_ZOP (see VARARGIN)

% Choose default command line output for bh_tomo_ZOP
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

h.mog = '';
h.av = [];
h.ap = [];
h.forageTx = [];
h.forageRx = [];
if nargin >=4
	h.mog = varargin{1};
end
if nargin >=5
	air = varargin{2};
    h.av = air( h.mog.av );
    h.ap = air( h.mog.ap );
end
set(handles.checkbox_scont,'Enable','off');
if ~isempty(h.mog)
    [h.indZop, h.t0] = plotZOP(h.mog, h.av, h.ap, handles.axes1, true,...
		get(handles.checkbox_spreadcorr,'Value'), h.mog.data.rstepsz*0.5);
	clim = caxis(handles.axes1);
	Amax = clim(2);
	set(handles.edit_offsetMax,'String',num2str(h.mog.data.rstepsz*0.5))
	set(handles.edit_sensitivity,'String',num2str(Amax))
	xlim = get(handles.axes1,'XLim');
	set(handles.edit_tmin,'String',num2str(xlim(1)))
	set(handles.edit_tmax,'String',num2str(xlim(2)))
end
if nargin >=6
    forages = varargin{3};
    h.forageTx = forages(h.mog.Tx);
    h.forageRx = forages(h.mog.Rx);
    if ~(isempty(h.forageTx) && isempty(h.forageRx))
        set(handles.checkbox_scont,'Enable','on');
    end
end
setappdata(handles.fig_bh_tomo_ZOP,'h',h)

updateAxes2(handles)


% UIWAIT makes bh_tomo_ZOP wait for user response (see UIRESUME)
% uiwait(handles.fig_bh_tomo_ZOP);


% --- Outputs from this function are returned to the command line.
function varargout = bh_tomo_ZOP_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


% --- Executes on selection change in popupmenu_sensitivity.
function popupmenu_sensitivity_Callback(hObject, eventdata, handles)
setCaxis(handles)

% --- Executes during object creation, after setting all properties.
function popupmenu_sensitivity_CreateFcn(hObject, eventdata, handles)
% hObject    handle to popupmenu_sensitivity (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc
    set(hObject,'BackgroundColor','white');
else
    set(hObject,'BackgroundColor',get(0,'defaultUicontrolBackgroundColor'));
end





function edit_tmin_Callback(hObject, eventdata, handles)
tmin = str2double(get(hObject,'String'));
xlim = get(handles.axes1,'XLim');
xlim(1) = tmin;
set(handles.axes1,'XLim',xlim);


% --- Executes during object creation, after setting all properties.
function edit_tmin_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit_tmax_Callback(hObject, eventdata, handles)
tmax = str2double(get(hObject,'String'));
xlim = get(handles.axes1,'XLim');
xlim(2) = tmax;
set(handles.axes1,'XLim',xlim);

% --- Executes during object creation, after setting all properties.
function edit_tmax_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function plotAppVel(handles)
h = getappdata(handles.fig_bh_tomo_ZOP,'h');
if isempty(h.mog)
    return
end

lrai = sqrt( (h.mog.data.Tx_x(h.indZop)-h.mog.data.Rx_x(h.indZop)).^2 + ...
    (h.mog.data.Tx_y(h.indZop)-h.mog.data.Rx_y(h.indZop)).^2 + ...
    (h.mog.data.Tx_z(h.indZop)-h.mog.data.Rx_z(h.indZop)).^2 );
tt = h.mog.tt(h.indZop)-h.t0(h.indZop);
et = h.mog.et(h.indZop);
vApp = lrai./tt;
vMin = lrai./(tt+et);
vMax = lrai./(tt-et);
Tx_z = h.mog.data.Tx_z(h.indZop);
ind = h.mog.tt(h.indZop)~=-1;
plot(handles.axes2, vApp(ind), Tx_z(ind), 'o')
z = Tx_z(ind);
vMin = vMin(ind);
vMax = vMax(ind);
hold(handles.axes2, 'on')
for n=1:length(vMin)
	plot(handles.axes2, [vMin(n) vMax(n)], [z(n) z(n)],'Color',[0.75 0.75 0.75])
end
hold(handles.axes2, 'off')
if get(handles.checkbox_scont,'Value')==1
    plot_scont(handles, handles.axes2)
end
grid(handles.axes2, 'on')
xlabel(handles.axes2, 'Velocity [m/ns]')
ylabel(handles.axes2, 'Elevation [m]')
title(handles.axes2, 'Apparent velocity')
set(handles.axes2,'YLim', get(handles.axes1,'YLim'))


function plotAmp(handles)
h = getappdata(handles.fig_bh_tomo_ZOP,'h');
if isempty(h.mog)
    return
end

lrai = sqrt( (h.mog.data.Tx_x(h.indZop)-h.mog.data.Rx_x(h.indZop)).^2 + ...
    (h.mog.data.Tx_y(h.indZop)-h.mog.data.Rx_y(h.indZop)).^2 + ...
    (h.mog.data.Tx_z(h.indZop)-h.mog.data.Rx_z(h.indZop)).^2 );
tmin = h.mog.amp_tmin(h.indZop);
tmax = h.mog.amp_tmax(h.indZop);
imin = zeros(size(tmin));
imax = zeros(size(tmax));
App = nan(size(tmin));
for n=1:length(tmin), imin(n) = findnear(tmin(n), h.mog.data.timestp); end
for n=1:length(tmax), imax(n) = findnear(tmax(n), h.mog.data.timestp); end
for n=1:length(tmin)
	if imax(n)>imin(n)
		App(n) = lrai(n)*(max(h.mog.data.rdata(imin:imax,h.indZop(n))) - ...
			min(h.mog.data.rdata(imin:imax,h.indZop(n))));
	end
end
Tx_z = h.mog.data.Tx_z(h.indZop);
axes(handles.axes2)
semilogx(App, Tx_z,'LineStyle','none','Marker','o','MarkerEdgeColor',[0 0.5 0])
xlabel('Amplitude')
ylabel('Elevation [m]')
title('Peak-to-peak amplitude')
if get(handles.checkbox_scont,'Value')==1
    plot_scont(handles, handles.axes2)
end
set(handles.axes2,'YLim', get(handles.axes1,'YLim'))




function plotVelAmp(handles)
h = getappdata(handles.fig_bh_tomo_ZOP,'h');
if isempty(h.mog)
    return
end

lrai = sqrt( (h.mog.data.Tx_x(h.indZop)-h.mog.data.Rx_x(h.indZop)).^2 + ...
    (h.mog.data.Tx_y(h.indZop)-h.mog.data.Rx_y(h.indZop)).^2 + ...
    (h.mog.data.Tx_z(h.indZop)-h.mog.data.Rx_z(h.indZop)).^2 );
tt = h.mog.tt(h.indZop)-h.t0(h.indZop);
et = h.mog.et(h.indZop);
vApp = lrai./tt;
vMin = lrai./(tt+et);
vMax = lrai./(tt-et);
Tx_z = h.mog.data.Tx_z(h.indZop);
ind = h.mog.tt(h.indZop)~=-1;
tmin = h.mog.amp_tmin(h.indZop);
tmax = h.mog.amp_tmax(h.indZop);
imin = zeros(size(tmin));
imax = zeros(size(tmax));
App = nan(size(tmin));
for n=1:length(tmin), imin(n) = findnear(tmin(n), h.mog.data.timestp); end
for n=1:length(tmax), imax(n) = findnear(tmax(n), h.mog.data.timestp); end
for n=1:length(tmin)
	if imax(n)>imin(n)
		App(n) = lrai(n)*(max(h.mog.data.rdata(imin(n):imax(n),h.indZop(n))) - ...
			min(h.mog.data.rdata(imin(n):imax(n),h.indZop(n))));
	end
end

axes(handles.axes2)
[ax,h1,h2] = plotxx(vApp(ind), Tx_z(ind), App, Tx_z, 'plot', 'semilogx');
set(get(ax(1),'Xlabel'),'String','Apparent velocity [m/ns]')
set(get(ax(2),'Xlabel'),'String','Peak-to-peak amplitude')
set(h1,'LineStyle','none','Marker','o','MarkerEdgeColor','b')
set(h2,'LineStyle','none','Marker','o','MarkerEdgeColor',[0 0.5 0])
set(ax(2),'YAxisLocation','right')
ylabel('Elevation [m]')

z = Tx_z(ind);
vMin = vMin(ind);
vMax = vMax(ind);
hold(ax(1),'on')
for n=1:length(vMin)
	plot([vMin(n) vMax(n)], [z(n) z(n)],'Color',[0.75 0.75 0.75])
end
hold(ax(1),'off')
if get(handles.checkbox_scont,'Value')==1
    plot_scont(handles, ax(1))
end
set(ax(1),'YLim', get(handles.axes1,'YLim'))
set(ax(2),'YLim', get(handles.axes1,'YLim'))


function setCaxis(handles)
axes(handles.axes1)
if get(handles.popupmenu_sensitivity,'Value')==1
	caxis([-28000 28000])
	set(handles.edit_sensitivity,'String','28000')
elseif get(handles.popupmenu_sensitivity,'Value')==2
	caxis([-7000 7000])
	set(handles.edit_sensitivity,'String','7000')
elseif get(handles.popupmenu_sensitivity,'Value')==3
	caxis([-700 700])
	set(handles.edit_sensitivity,'String','700')
end


function edit_sensitivity_Callback(hObject, eventdata, handles)
axes(handles.axes1)
clim = str2double(get(hObject,'String'));
caxis([-clim clim])

function edit_sensitivity_CreateFcn(hObject, eventdata, handles)

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end




function pushbutton_print_Callback(hObject, eventdata, handles)
% printdlg(handles.fig_bh_tomo_ZOP, [.13 .11 .62 .815],...
% 	get(handles.fig_bh_tomo_ZOP,'defaultaxesposition'),...
% 	[handles.uipanel_control handles.text_tmin handles.text_tmax ...
% 	handles.edit_tmin handles.edit_tmax handles.edit_sensitivity ...
% 	handles.popupmenu_sensitivity handles.pushbutton_print])
h = getappdata(handles.fig_bh_tomo_ZOP,'h');
nomzop = ['ZOP_',h.mog.name];
print('-depsc2','-noui',nomzop) 


function checkbox_spreadcorr_Callback(hObject, eventdata, handles)
h = getappdata(handles.fig_bh_tomo_ZOP,'h');
[h.indZop, h.t0] = plotZOP(h.mog, h.av, h.ap, handles.axes1, true,...
	get(hObject,'Value'));
clim = caxis(handles.axes1);
Amax = clim(2);
set(handles.edit_sensitivity,'String',num2str(Amax))
xlim(1) = str2double(get(handles.edit_tmin,'String'));
xlim(2) = str2double(get(handles.edit_tmax,'String'));
set(handles.axes1,'XLim',xlim)


function edit_offsetMax_Callback(hObject, eventdata, handles)
offMax = str2double(get(hObject,'String'));
h = getappdata(handles.fig_bh_tomo_ZOP,'h');
[h.indZop, h.t0] = plotZOP(h.mog, h.av, h.ap, handles.axes1, true,...
	get(handles.checkbox_spreadcorr,'Value'), offMax);
clim = caxis(handles.axes1);
Amax = clim(2);
set(handles.edit_sensitivity,'String',num2str(Amax))
setappdata(handles.fig_bh_tomo_ZOP,'h',h)
updateAxes2(handles)
xlim(1) = str2double(get(handles.edit_tmin,'String'));
xlim(2) = str2double(get(handles.edit_tmax,'String'));
set(handles.axes1,'XLim',xlim)


function edit_offsetMax_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end




function pushbutton_rais_Callback(hObject, eventdata, handles)
h = getappdata(handles.fig_bh_tomo_ZOP,'h');
ind = h.indZop;
xy = sqrt( (h.mog.data.Tx_x-h.mog.data.Rx_x).^2 + ...
    (h.mog.data.Tx_y-h.mog.data.Rx_y).^2 );

figure
plot([zeros(size(xy(ind)))' xy(ind)']', ...
    [h.mog.data.Tx_z(ind)' h.mog.data.Rx_z(ind)']','g')
set(gca, 'DataAspectRatio',[1 1 1])
axis tight
xlabel('Tx-Rx Distance [m]')
ylabel('Elevation [m]')




function edit_zmin_Callback(hObject, eventdata, handles)
zmin = str2double(get(hObject,'String'));
if isempty(zmin), return, end
zlim = get(handles.axes1,'YLim');
zlim(1) = zmin;
set(handles.axes1,'YLim',zlim);
set(handles.axes2,'YLim',zlim);


function edit_zmin_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


function edit_zmax_Callback(hObject, eventdata, handles)
zmax = str2double(get(hObject,'String'));
if isempty(zmax), return, end
zlim = get(handles.axes1,'YLim');
zlim(2) = zmax;
set(handles.axes1,'YLim',zlim);
set(handles.axes2,'YLim',zlim);


function edit_zmax_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


function updateAxes2(handles)
if get(handles.checkbox_velocity,'Value')==1 && ...
		get(handles.checkbox_attenuation,'Value')==1
	plotVelAmp(handles)
elseif get(handles.checkbox_velocity,'Value')==1
	plotAppVel(handles)
elseif get(handles.checkbox_attenuation,'Value')==1
	plotAmp(handles)
else
	cla(handles.axes2)
end



function checkbox_attenuation_Callback(hObject, eventdata, handles)
updateAxes2(handles)

function checkbox_velocity_Callback(hObject, eventdata, handles)
updateAxes2(handles)



% --- Executes on button press in checkbox_scont.
function checkbox_scont_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox_scont (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkbox_scont
updateAxes2(handles)


function plot_scont(handles, axes)
h = getappdata(handles.fig_bh_tomo_ZOP,'h');

hold(axes,'on')
if ~isempty( h.forageTx.scont )
    plot(axes, 1./h.forageTx.scont.valeur, h.forageTx.scont.z,'b')
end
if ~isempty( h.forageRx.scont )
    plot(axes, 1./h.forageRx.scont.valeur, h.forageRx.scont.z,'g')
end    
hold(axes,'off')
