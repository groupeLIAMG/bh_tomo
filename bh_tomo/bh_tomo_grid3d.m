function varargout = bh_tomo_grid3d(varargin)
% BH_TOMO_GRID3D MATLAB code for bh_tomo_grid3d.fig
%      BH_TOMO_GRID3D, by itself, creates a new BH_TOMO_GRID3D or raises the existing
%      singleton*.
%
%      H = BH_TOMO_GRID3D returns the handle to a new BH_TOMO_GRID3D or the handle to
%      the existing singleton*.
%
%      BH_TOMO_GRID3D('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in BH_TOMO_GRID3D.M with the given input arguments.
%
%      BH_TOMO_GRID3D('Property','Value',...) creates a new BH_TOMO_GRID3D or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before bh_tomo_grid3d_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to bh_tomo_grid3d_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES
% Edit the above text to modify the response to help bh_tomo_grid3d

% Last Modified by GUIDE v2.5 07-Jan-2013 14:52:10

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @bh_tomo_grid3d_OpeningFcn, ...
                   'gui_OutputFcn',  @bh_tomo_grid3d_OutputFcn, ...
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


% --- Executes just before bh_tomo_grid3d is made visible.
function bh_tomo_grid3d_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to bh_tomo_grid3d (see VARARGIN)

% Choose default command line output for bh_tomo_grid3d
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

data = [];
g.grx = [];
g.gry = [];
g.grz = [];
g.cont.slowness.data = [];
g.cont.slowness.data_xi = [];
g.cont.attenuation.data = [];
g.Tx = [];
g.Rx = [];
g.TxCosDir = [];
g.RxCosDir = [];
g.x0 = [];
g.bord = [1 1 1 1];  % [nxm nxp nzm nzp]
g.Tx_Z_water = [];
g.Rx_Z_water = [];
g.in = [];
if nargin >= 4
    data = varargin{1};
end
if nargin >= 5
    g = varargin{2};
end
setappdata(handles.fig_bh_grid3d,'data',data)
setappdata(handles.fig_bh_grid3d,'g',g)
update_plot3D(handles)

set(handles.edit_nxm,'String',num2str(g.bord(1)))
set(handles.edit_nxp,'String',num2str(g.bord(2)))
set(handles.edit_nym,'String',num2str(g.bord(3)))
set(handles.edit_nyp,'String',num2str(g.bord(4)))
if isempty(g.x0)
    set(handles.edit_nx,'String','10')
    set(handles.edit_ny,'String','10')
    set(handles.edit_nz,'String','10')
else
    set(handles.edit_dx,'String',num2str(g.grx(2)-g.grx(1)))
    set(handles.edit_dy,'String',num2str(g.gry(2)-g.gry(1)))
    set(handles.edit_dz,'String',num2str(g.grz(2)-g.grz(1)))
    set(handles.edit_nx,'String',num2str(length(g.grx)))
    set(handles.edit_ny,'String',num2str(length(g.gry)))
    set(handles.edit_nz,'String',num2str(length(g.grz)))
end

str = get_str_locale();
setappdata(handles.fig_bh_grid3d,'str',str)
set_String_locale(handles, str)

if ~isempty(data)
    if isempty(g.x0)
        xmin = min([data.Tx(data.in,1);data.Rx(data.in,1)]);
        ymin = min([data.Tx(data.in,2);data.Rx(data.in,2)]);
        zmin = min([data.Tx(data.in,3);data.Rx(data.in,3)]);
        
        set(handles.edit_x0,'String',num2str(xmin))
        set(handles.edit_y0,'String',num2str(ymin))
        set(handles.edit_z0,'String',num2str(zmin))
		g.x0 = [xmin ymin zmin];
        setappdata(handles.fig_bh_grid3d,'g',g)
    else
        set(handles.edit_x0,'String',num2str(g.x0(1)))
        set(handles.edit_y0,'String',num2str(g.x0(2)))
        set(handles.edit_z0,'String',num2str(g.x0(3)))
    end
    if isempty(g.grx)
        init_grid(handles)
    end
    update_info(handles)
end

% UIWAIT makes bh_tomo_grid3d wait for user response (see UIRESUME)
uiwait(handles.fig_bh_grid3d);

% --- Outputs from this function are returned to the command line.
function varargout = bh_tomo_grid3d_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = hObject;
if isfield(handles, 'fig_bh_grid3d')
    g = getappdata(handles.fig_bh_grid3d,'g');
    varargout{2} = g;
else
    varargout{2} = [];
end


function init_grid(handles)
g = getappdata(handles.fig_bh_grid3d,'g');
data = getappdata(handles.fig_bh_grid3d,'data');
g.Tx = data.Tx;
g.Rx = data.Rx;
g.in = data.in;

xmin = str2double(get(handles.edit_x0,'String'));
xmax = max([g.Tx(g.in,1);g.Rx(g.in,1)]);
ymin = str2double(get(handles.edit_y0,'String'));
ymax = max([g.Tx(g.in,2);g.Rx(g.in,2)]);
zmin = str2double(get(handles.edit_z0,'String'));
zmax = max([g.Tx(g.in,3);g.Rx(g.in,3)]);

nx = str2double(get(handles.edit_nx,'String'));
ny = str2double(get(handles.edit_ny,'String'));
nz = str2double(get(handles.edit_nz,'String'));

dx = (xmax-xmin)/nx;
dy = (ymax-ymin)/ny;
dz = (zmax-zmin)/nz;

nxm = str2double(get(handles.edit_nxm,'String'));
nxp = str2double(get(handles.edit_nxp,'String'));
g.grx = xmin+dx*((0-nxm):(nx+nxp));
g.grx = g.grx';

nym = str2double(get(handles.edit_nym,'String'));
nyp = str2double(get(handles.edit_nyp,'String'));
g.gry = ymin+dy*((0-nym):(ny+nyp));
g.gry = g.gry';

nzm = str2double(get(handles.edit_nzm,'String'));
nzp = str2double(get(handles.edit_nzp,'String'));
g.grz = zmin+dz*((0-nzm):(nz+nzp));
g.grz = g.grz';

setappdata(handles.fig_bh_grid3d,'g',g)



function update_plot3D(handles)
data = getappdata(handles.fig_bh_grid3d,'data');
if isempty(data)
    return
end
axes(handles.axes_plot3D)
set(handles.axes_plot3D,'NextPlot','add')

for n=1:length(data.boreholes)
    plot3([data.boreholes(n).X data.boreholes(n).Xmax], ...
        [data.boreholes(n).Y data.boreholes(n).Ymax], ...
        [data.boreholes(n).Z data.boreholes(n).Zmax], 'r')
end
plot3(data.Tx(data.in,1), data.Tx(data.in,2), data.Tx(data.in,3),'b.')
plot3(data.Rx(data.in,1), data.Rx(data.in,2), data.Rx(data.in,3),'g.')
view(3)
axis equal
grid on
xlabel('X')
ylabel('Y')
zlabel('Z')
set(handles.axes_plot3D,'NextPlot','replace')


function togglebutton_rotate3d_Callback(hObject, eventdata, handles)
button_state = get(hObject,'Value');
if button_state == get(hObject,'Max')
    % toggle button is pressed
	rotate3d(handles.axes_plot3D,'on')
elseif button_state == get(hObject,'Min')
    % toggle button is not pressed
	rotate3d(handles.axes_plot3D,'off')
end

function togglebutton_zoom_Callback(hObject, eventdata, handles)
button_state = get(hObject,'Value');
if button_state == get(hObject,'Max')
    % toggle button is pressed
	zoom(handles.axes_plot3D,'on')
elseif button_state == get(hObject,'Min')
    % toggle button is not pressed
	zoom(handles.axes_plot3D,'off')
end


function set_String_locale(handles, str)
set(handles.uipanel_info,             'Title',  str.s133)
set(handles.uipanel_grid,           'Title',  str.s92)
set(handles.text_pas,                 'String', lower(str.s162))
set(handles.pushbutton_contraintes,   'String', str.s190)
set(handles.pushbutton_quitter,       'String', str.s193)
set(handles.text_cell_p,              'String', [str.s191,' +'])
set(handles.text_cell_m,              'String', [str.s191,' -'])
set(handles.text_origine,             'String', str.s194)

load rotate
set(handles.togglebutton_rotate3d,'CData',cdata)
load zoom
set(handles.togglebutton_zoom,'CData',zoomCData)
clear cdata zoomCData



function edit_nxm_Callback(hObject, eventdata, handles)
nxm = str2double(get(hObject,'string'));
if isnan(nxm)
    str = getappdata(handles.fig_bh_grid3d,'str');
    errordlg(str.s54,str.s45,'modal')
    return
end
g = getappdata(handles.fig_bh_grid3d,'g');
g.bord(1) = nxm;

nxp = str2double(get(handles.edit_nxp,'String'));
xmin = str2double(get(handles.edit_x0,'String'));
nx = str2double(get(handles.edit_nx,'String'));
dx = str2double(get(handles.edit_dx,'String'));
g.grx = xmin+dx*((0-nxm):(nx+nxp));
g.grx = g.grx';

setappdata(handles.fig_bh_grid3d,'g',g)
update_info(handles)

function edit_nxm_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


function edit_nxp_Callback(hObject, eventdata, handles)
nxp = str2double(get(hObject,'string'));
if isnan(nxp)
    str = getappdata(handles.fig_bh_grid3d,'str');
    errordlg(str.s54,str.s45,'modal')
    return
end
g = getappdata(handles.fig_bh_grid3d,'g');
g.bord(2) = nxp;

nxm = str2double(get(handles.edit_nxm,'String'));
xmin = str2double(get(handles.edit_x0,'String'));
nx = str2double(get(handles.edit_nx,'String'));
dx = str2double(get(handles.edit_dx,'String'));
g.grx = xmin+dx*((0-nxm):(nx+nxp));
g.grx = g.grx';

setappdata(handles.fig_bh_grid3d,'g',g)
update_info(handles)

function edit_nxp_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


function edit_dx_Callback(hObject, eventdata, handles)
dx = str2double(get(hObject,'string'));
if isnan(dx)
    str = getappdata(handles.fig_bh_grid3d,'str');
    errordlg(str.s54,str.s45,'modal')
    return
end
g = getappdata(handles.fig_bh_grid3d,'g');
xmin = str2double(get(handles.edit_x0,'String'));
xmax = max([g.Tx(g.in,1);g.Rx(g.in,1)]) + 0.5*dx;
nx = ceil((xmax-xmin)/dx);
set(handles.edit_nx,'String',num2str(nx))

nxm = str2double(get(handles.edit_nxm,'String'));
nxp = str2double(get(handles.edit_nxp,'String'));
g.grx = xmin+dx*((0-nxm):(nx+nxp));
g.grx = g.grx';
setappdata(handles.fig_bh_grid3d,'g',g)
update_info(handles)

function edit_dx_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


function edit_nym_Callback(hObject, eventdata, handles)
nym = str2double(get(hObject,'string'));
if isnan(nym)
    str = getappdata(handles.fig_bh_grid3d,'str');
    errordlg(str.s54,str.s45,'modal')
    return
end
g = getappdata(handles.fig_bh_grid3d,'g');
g.bord(3) = nym;

nyp = str2double(get(handles.edit_nyp,'String'));
ymin = str2double(get(handles.edit_y0,'String'));
ny = str2double(get(handles.edit_ny,'String'));
dy = str2double(get(handles.edit_dy,'String'));
g.gry = ymin+dy*((0-nym):(ny+nyp));
g.gry = g.gry';

setappdata(handles.fig_bh_grid3d,'g',g)
update_info(handles)

function edit_nym_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


function edit_nyp_Callback(hObject, eventdata, handles)
nyp = str2double(get(hObject,'string'));
if isnan(nyp)
    str = getappdata(handles.fig_bh_grid3d,'str');
    errordlg(str.s54,str.s45,'modal')
    return
end
g = getappdata(handles.fig_bh_grid3d,'g');
g.bord(4) = nyp;

nym = str2double(get(handles.edit_nym,'String'));
ymin = str2double(get(handles.edit_y0,'String'));
ny = str2double(get(handles.edit_ny,'String'));
dy = str2double(get(handles.edit_dy,'String'));
g.gry = ymin+dy*((0-nym):(ny+nyp));
g.gry = g.gry';

setappdata(handles.fig_bh_grid3d,'g',g)
update_info(handles)

function edit_nyp_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


function edit_dy_Callback(hObject, eventdata, handles)
dy = str2double(get(hObject,'string'));
if isnan(dy)
    str = getappdata(handles.fig_bh_grid3d,'str');
    errordlg(str.s54,str.s45,'modal')
    return
end
g = getappdata(handles.fig_bh_grid3d,'g');
ymin = str2double(get(handles.edit_y0,'String'));
ymax = max([g.Tx(g.in,2);g.Rx(g.in,2)]) + 0.5*dy;
ny = ceil((ymax-ymin)/dy);
set(handles.edit_ny,'String',num2str(ny))

nym = str2double(get(handles.edit_nym,'String'));
nyp = str2double(get(handles.edit_nyp,'String'));
g.gry = ymin+dy*((0-nym):(ny+nyp));
g.gry = g.gry';
setappdata(handles.fig_bh_grid3d,'g',g)
update_info(handles)

function edit_dy_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function edit_nzm_Callback(hObject, eventdata, handles)
nzm = str2double(get(hObject,'string'));
if isnan(nzm)
    str = getappdata(handles.fig_bh_grid3d,'str');
    errordlg(str.s54,str.s45,'modal')
    return
end
g = getappdata(handles.fig_bh_grid3d,'g');
g.bord(3) = nzm;

nzp = str2double(get(handles.edit_nzp,'String'));
zmin = str2double(get(handles.edit_z0,'String'));
nz = str2double(get(handles.edit_nz,'String'));
dz = str2double(get(handles.edit_dz,'String'));
g.grz = zmin+dz*((0-nzm):(nz+nzp));
g.grz = g.grz';

setappdata(handles.fig_bh_grid3d,'g',g)
update_info(handles)


function edit_nzm_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit_nzp_Callback(hObject, eventdata, handles)
nzp = str2double(get(hObject,'string'));
if isnan(nzp)
    str = getappdata(handles.fig_bh_grid3d,'str');
    errordlg(str.s54,str.s45,'modal')
    return
end
g = getappdata(handles.fig_bh_grid3d,'g');
g.bord(4) = nzp;

nzm = str2double(get(handles.edit_nzm,'String'));
zmin = str2double(get(handles.edit_z0,'String'));
nz = str2double(get(handles.edit_nz,'String'));
dz = str2double(get(handles.edit_dz,'String'));
g.grz = zmin+dz*((0-nzm):(nz+nzp));
g.grz = g.grz';

setappdata(handles.fig_bh_grid3d,'g',g)
update_info(handles)


function edit_nzp_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


function edit_dz_Callback(hObject, eventdata, handles)
dz = str2double(get(hObject,'string'));
if isnan(dz)
    str = getappdata(handles.fig_bh_grid3d,'str');
    errordlg(str.s54,str.s45,'modal')
    return
end
g = getappdata(handles.fig_bh_grid3d,'g');
zmin = str2double(get(handles.edit_z0,'String'));
zmax = max([g.Tx(g.in,3);g.Rx(g.in,3)]) + 0.5*dz;
nz = ceil((zmax-zmin)/dz);
set(handles.edit_nz,'String',num2str(nz))

nzm = str2double(get(handles.edit_nzm,'String'));
nzp = str2double(get(handles.edit_nzp,'String'));
g.grz = zmin+dz*((0-nzm):(nz+nzp));
g.grz=g.grz';
setappdata(handles.fig_bh_grid3d,'g',g)
update_info(handles)

function edit_dz_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end




function edit_x0_Callback(hObject, eventdata, handles)
val = str2double(get(hObject,'string'));
if isnan(val)
    errordlg(str.s54,str.s45,'modal')
    return
end
xmin = min([g.Tx(g.in,1);g.Rx(g.in,1)]);
if val>xmin
    errordlg('Value larger than borehole x min',str.s45,'modal')
    set(hobject,'String',num2str(xmin))
end

xmax = max([g.Tx(g.in,1);g.Rx(g.in,1)]) + 0.5*dx;
nx = ceil((xmax-xmin)/dx);
set(handles.edit_nx,'String',num2str(nx))

nxm = str2double(get(handles.edit_nxm,'String'));
nxp = str2double(get(handles.edit_nxp,'String'));
g.grx = xmin+dx*((0-nxm):(nx+nxp));
g.grx = g.grx';
setappdata(handles.fig_bh_grid3d,'g',g)
update_info(handles)


function edit_x0_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


function edit_y0_Callback(hObject, eventdata, handles)
val = str2double(get(hObject,'string'));
if isnan(val)
    errordlg(str.s54,str.s45,'modal')
    return
end
ymin = min([g.Tx(g.in,2);g.Rx(g.in,2)]);
if val>ymin
    errordlg('Value larger than borehole y min',str.s45,'modal')
    set(hobject,'String',num2str(ymin))
end

ymax = max([g.Tx(g.in,2);g.Rx(g.in,2)]) + 0.5*dy;
ny = ceil((ymax-ymin)/dy);
set(handles.edit_ny,'String',num2str(ny))

nym = str2double(get(handles.edit_nym,'String'));
nyp = str2double(get(handles.edit_nyp,'String'));
g.gry = ymin+dy*((0-nym):(ny+nyp));
g.gry = g.gry';
setappdata(handles.fig_bh_grid3d,'g',g)
update_info(handles)


function edit_y0_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


function edit_z0_Callback(hObject, eventdata, handles)
val = str2double(get(hObject,'string'));
if isnan(val)
    errordlg(str.s54,str.s45,'modal')
    return
end
zmin = min([g.Tx(g.in,3);g.Rx(g.in,3)]);
if val>zmin
    errordlg('Value larger than borehole z min',str.s45,'modal')
    set(hobject,'String',num2str(zmin))
end

zmax = max([g.Tx(g.in,3);g.Rx(g.in,3)]) + 0.5*dz;
nz = ceil((zmax-zmin)/dz);
set(handles.edit_nz,'String',num2str(nz))

nzm = str2double(get(handles.edit_nzm,'String'));
nzp = str2double(get(handles.edit_nzp,'String'));
g.grz = zmin+dz*((0-nzm):(nz+nzp));
g.grz=g.grz';
setappdata(handles.fig_bh_grid3d,'g',g)
update_info(handles)


function edit_z0_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end




function update_info(handles)
str = getappdata(handles.fig_bh_grid3d,'str');
g = getappdata(handles.fig_bh_grid3d,'g');
data = getappdata(handles.fig_bh_grid3d,'data');

nx = length(g.grx);
ny = length(g.gry);
nz = length(g.grz);
texte = {'';
    [str.s34,': ',num2str(nx*ny*nz)];
    [str.s130,': ',num2str(size(data.Tx,1))]};
set(handles.text_info,'String',texte)
set(handles.edit_dx,'String',num2str(g.grx(2)-g.grx(1)))
set(handles.edit_dy,'String',num2str(g.gry(2)-g.gry(1)))
set(handles.edit_dz,'String',num2str(g.grz(2)-g.grz(1)))
set(handles.edit_nx,'String',num2str(nx))
set(handles.edit_ny,'String',num2str(ny))
set(handles.edit_nz,'String',num2str(nz))


function pushbutton_contraintes_Callback(hObject, eventdata, handles)
data = getappdata(handles.fig_bh_grid3d,'data');
g = getappdata(handles.fig_bh_grid3d,'g');
h = getappdata(handles.fig_bh_grid3d,'h');
grx = [g.grx(1) g.grx(2)-g.grx(1) g.grx(length(g.grx))];
gry = [g.gry(1) g.gry(2)-g.gry(1) g.gry(length(g.gry))];
grz = [g.grz(1) g.grz(2)-g.grz(1) g.grz(length(g.grz))];

for n=numel(data.boreholes)
	if isfield( data.boreholes(n), 'scont' )
		g.cont.slowness = [g.cont.slowness; ...
			data.boreholes(n).scont];
	end
	if isfield( data.boreholes(n), 'acont' )
		g.cont.attenuation = [g.cont.attenuation; ...
			data.boreholes(n).acont];
	end
end

[hh, g.cont] = bh_tomo_contraintes3d(grx, gry, grz, g.cont);
if ishandle(hh), close(hh); end
setappdata(handles.fig_bh_grid3d,'g',g)


function pushbutton_quitter_Callback(hObject, eventdata, handles)
guidata(hObject, handles);
uiresume;


function fig_bh_grid3d_CloseRequestFcn(hObject, eventdata, handles)
guidata(hObject, handles);
uiresume(handles.fig_bh_grid3d);
%delete(hObject);




function edit_nx_Callback(hObject, eventdata, handles)
nx = str2double(get(hObject,'String'));
if isnan(nx)
    str = getappdata(handles.fig_bh_grid3d,'str');
    errordlg(str.s54,str.s45,'modal')
    return
end
g = getappdata(handles.fig_bh_grid3d,'g');
xmin = str2double(get(handles.edit_x0,'String'));
xmax = max([g.Tx(g.in,1);g.Rx(g.in,1)]);
dx = (xmax-xmin)/nx;
xmax = xmax+0.5*dx;
dx = (xmax-xmin)/nx;
set(handles.edit_dx,'String',num2str(dx))

nxm = str2double(get(handles.edit_nxm,'String'));
nxp = str2double(get(handles.edit_nxp,'String'));
g.grx = xmin+dx*((0-nxm):(nx+nxp));
g.grx=g.grx';
setappdata(handles.fig_bh_grid3d,'g',g)
update_info(handles)

function edit_nx_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function edit_ny_Callback(hObject, eventdata, handles)
ny = str2double(get(hObject,'String'));
if isnan(ny)
    str = getappdata(handles.fig_bh_grid3d,'str');
    errordlg(str.s54,str.s45,'modal')
    return
end
g = getappdata(handles.fig_bh_grid3d,'g');
ymin = str2double(get(handles.edit_y0,'String'));
ymax = max([g.Tx(g.in,2);g.Rx(g.in,2)]);
dy = (ymax-ymin)/ny;
ymax = ymax+0.5*dy;
dy = (ymax-ymin)/ny;
set(handles.edit_dy,'String',num2str(dy))

nym = str2double(get(handles.edit_nym,'String'));
nyp = str2double(get(handles.edit_nyp,'String'));
g.gry = ymin+dy*((0-nym):(ny+nyp));
g.gry=g.gry';
setappdata(handles.fig_bh_grid3d,'g',g)
update_info(handles)


function edit_ny_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function edit_nz_Callback(hObject, eventdata, handles)
nz = str2double(get(hObject,'String'));
if isnan(nz)
    str = getappdata(handles.fig_bh_grid3d,'str');
    errordlg(str.s54,str.s45,'modal')
    return
end
g = getappdata(handles.fig_bh_grid3d,'g');
zmin = str2double(get(handles.edit_z0,'String'));
zmax = max([g.Tx(g.in,3);g.Rx(g.in,3)]);
dz = (zmax-zmin)/nz;
zmax = zmax + 0.5*dz;
dz = (zmax-zmin)/nz;
set(handles.edit_dz,'String',num2str(dz))

nzm = str2double(get(handles.edit_nzm,'String'));
nzp = str2double(get(handles.edit_nzp,'String'));
g.grz = zmin+dz*((0-nzm):(nz+nzp));
g.grz=g.grz';
setappdata(handles.fig_bh_grid3d,'g',g)
%plot_proj(handles)
update_info(handles)

function edit_nz_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
