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
g.borehole_x0 = [];
g.bord = [1 1 1 1];  % [nxm nxp nzm nzp]
g.flip = 0;
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
    tmp = cell(1,length(data.boreholes));
    for n=1:length(data.boreholes)
        tmp{n} = char( data.boreholes(n).name );
    end
    set(handles.popupmenu_origine,'String',tmp)
    if ~isempty(g.borehole_x0)
        set(handles.popupmenu_origine,'Value',g.borehole_x0)
    end
    if isempty(g.x0)
        set(handles.edit_x0,'String',num2str(data.boreholes(1).X))
        set(handles.edit_y0,'String',num2str(data.boreholes(1).Y))
        set(handles.edit_z0,'String',num2str(data.boreholes(1).Z))
		g.borehole_x0 = 1;
		g.x0 = [data.boreholes(1).X data.boreholes(1).Y data.boreholes(1).Z];
        setappdata(handles.fig_bh_grid3d,'g',g)
    else
        set(handles.edit_x0,'String',num2str(g.x0(1)))
        set(handles.edit_y0,'String',num2str(g.x0(2)))
        set(handles.edit_z0,'String',num2str(g.x0(3)))
    end
    trouve_plan(handles)
    if isempty(g.grx)
        update_proj(handles)
%     else
%         plot_proj(handles)
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
set(handles.text_f_origine,           'String', str.s192)
set(handles.text_origine,             'String', str.s194)

load rotate
set(handles.togglebutton_rotate3d,'CData',cdata)
load zoom
set(handles.togglebutton_zoom,'CData',zoomCData)
clear cdata zoomCData


function popupmenu_origine_Callback(hObject, eventdata, handles)
data = getappdata(handles.fig_bh_grid3d,'data');
g = getappdata(handles.fig_bh_grid3d,'g');
no = get(hObject,'Value');
g.x0(1) = data.boreholes(no).X;
g.x0(2) = data.boreholes(no).Y;
g.x0(3) = data.boreholes(no).Z;
g.borehole_x0 = no;
set(handles.edit_x0,'String',num2str(data.boreholes(no).X))
set(handles.edit_y0,'String',num2str(data.boreholes(no).Y))
set(handles.edit_z0,'String',num2str(data.boreholes(no).Z))
setappdata(handles.fig_bh_grid3d,'g',g)
update_proj(handles)
update_info(handles)

function popupmenu_origine_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


function edit_nxm_Callback(hObject, eventdata, handles)
val = str2double(get(hObject,'string'));
if isnan(val)
    str = getappdata(handles.fig_bh_grid3d,'str');
    errordlg(str.s54,str.s45,'modal')
    return
end
g = getappdata(handles.fig_bh_grid3d,'g');
g.bord(1) = val;
setappdata(handles.fig_bh_grid3d,'g',g)
update_proj(handles)
update_info(handles)

function edit_nxm_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


function edit_nxp_Callback(hObject, eventdata, handles)
val = str2double(get(hObject,'string'));
if isnan(val)
    str = getappdata(handles.fig_bh_grid3d,'str');
    errordlg(str.s54,str.s45,'modal')
    return
end
g = getappdata(handles.fig_bh_grid3d,'g');
g.bord(2) = val;
setappdata(handles.fig_bh_grid3d,'g',g)
update_proj(handles)
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
xmin = min([g.Tx(g.in,1);g.Rx(g.in,1)]) - 0.5*dx;
xmax = max([g.Tx(g.in,1);g.Rx(g.in,1)]) + 0.5*dx;
nx = ceil((xmax-xmin)/dx);
set(handles.edit_nx,'String',num2str(nx))

nxm = str2double(get(handles.edit_nxm,'String'));
nxp = str2double(get(handles.edit_nxp,'String'));
g.grx = xmin+dx*((0-nxm):(nx+nxp));
g.grx = g.grx';
setappdata(handles.fig_bh_grid3d,'g',g)
%plot_proj(handles)
update_info(handles)

function edit_dx_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


function edit_nym_Callback(hObject, eventdata, handles)
val = str2double(get(hObject,'string'));
if isnan(val)
    str = getappdata(handles.fig_bh_grid3d,'str');
    errordlg(str.s54,str.s45,'modal')
    return
end
g = getappdata(handles.fig_bh_grid3d,'g');
g.bord(3) = val;
setappdata(handles.fig_bh_grid3d,'g',g)
update_proj(handles)
update_info(handles)

function edit_nym_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


function edit_nyp_Callback(hObject, eventdata, handles)
val = str2double(get(hObject,'string'));
if isnan(val)
    str = getappdata(handles.fig_bh_grid3d,'str');
    errordlg(str.s54,str.s45,'modal')
    return
end
g = getappdata(handles.fig_bh_grid3d,'g');
g.bord(4) = val;
setappdata(handles.fig_bh_grid3d,'g',g)
update_proj(handles)
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
ymin = min([g.Tx(g.in,1);g.Rx(g.in,1)]) - 0.5*dy;
ymax = max([g.Tx(g.in,1);g.Rx(g.in,1)]) + 0.5*dy;
ny = ceil((ymax-ymin)/dy);
set(handles.edit_ny,'String',num2str(ny))

nym = str2double(get(handles.edit_nym,'String'));
nyp = str2double(get(handles.edit_nyp,'String'));
g.gry = ymin+dy*((0-nym):(ny+nyp));
g.gry = g.gry';
setappdata(handles.fig_bh_grid3d,'g',g)
%plot_proj(handles)
update_info(handles)

function edit_dy_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function edit_nzm_Callback(hObject, eventdata, handles)
val = str2double(get(hObject,'string'));
if isnan(val)
    str = getappdata(handles.fig_bh_grid3d,'str');
    errordlg(str.s54,str.s45,'modal')
    return
end
g = getappdata(handles.fig_bh_grid3d,'g');
g.bord(3) = val;
setappdata(handles.fig_bh_grid3d,'g',g)
update_proj(handles)
update_info(handles)


function edit_nzm_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit_nzp_Callback(hObject, eventdata, handles)
val = str2double(get(hObject,'string'));
if isnan(val)
    str = getappdata(handles.fig_bh_grid3d,'str');
    errordlg(str.s54,str.s45,'modal')
    return
end
g = getappdata(handles.fig_bh_grid3d,'g');
g.bord(4) = val;
setappdata(handles.fig_bh_grid3d,'g',g)
update_proj(handles)
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
zmin = min([g.Tx(g.in,3);g.Rx(g.in,3)]) - 0.5*dz;
zmax = max([g.Tx(g.in,3);g.Rx(g.in,3)]) + 0.5*dz;
nz = ceil((zmax-zmin)/dz);
set(handles.edit_nz,'String',num2str(nz))

nzm = str2double(get(handles.edit_nzm,'String'));
nzp = str2double(get(handles.edit_nzp,'String'));
g.grz = zmin+dz*((0-nzm):(nz+nzp));
g.grz=g.grz';
setappdata(handles.fig_bh_grid3d,'g',g)
%plot_proj(handles)
update_info(handles)

function edit_dz_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


function trouve_plan(handles)
data = getappdata(handles.fig_bh_grid3d,'data');
if isempty(data)
    return
end
% X = [data.boreholes(1).X data.boreholes(1).Y data.boreholes(1).Z;
%     data.boreholes(1).Xmax data.boreholes(1).Ymax data.boreholes(1).Zmax;
%     data.boreholes(2).X data.boreholes(2).Y data.boreholes(2).Z;
%     data.boreholes(2).Xmax data.boreholes(2).Ymax data.boreholes(2).Zmax];

%X = [data.boreholes(1).fdata; data.boreholes(2).fdata];

% trouve le plan qui passe par les points
uTx = unique(data.Tx(data.in,:), 'rows');
uRx = unique(data.Rx(data.in,:), 'rows');
%[h.x0, h.a]=lsplane([data.boreholes(1).fdata; data.boreholes(2).fdata]);

[h.x0,h.a]=lsplane([uTx; uRx]);
h.x0 = h.x0';  % Centroid of the data = point on the best-fit plane
h.a = h.a';    % Direction cosines of the normal to the best-fit plane
if h.a(3)<0, h.a=-h.a; end

% on projète les positions des Tx et Rx
h.Tx = proj_plan(data.Tx, h.x0, h.a);
h.Rx = proj_plan(data.Rx, h.x0, h.a);

% figure
% plot3(data.Tx(:,1),data.Tx(:,2),data.Tx(:,3), 'g+')
% hold on
% plot3(data.Rx(:,1),data.Rx(:,2),data.Rx(:,3), '+')
% plot3(h.Tx(:,1),h.Tx(:,2),h.Tx(:,3), 'go')
% plot3(h.Rx(:,1),h.Rx(:,2),h.Rx(:,3), 'o')
% hold off

setappdata(handles.fig_bh_grid3d,'h',h)

function update_proj(handles)
h = getappdata(handles.fig_bh_grid3d,'h');
g = getappdata(handles.fig_bh_grid3d,'g');
data = getappdata(handles.fig_bh_grid3d,'data');
% origine(1) = str2double(get(handles.edit_x0,'String'));
% origine(2) = str2double(get(handles.edit_y0,'String'));
% origine(3) = str2double(get(handles.edit_z0,'String'));
% [az,dip] = calcul_azimuth_dip(handles);

g.Tx = h.Tx;%transl_rotat(h.Tx, origine, az, dip);
g.Rx = h.Rx;%transl_rotat(h.Rx, origine, az, dip);
% g.TxCosDir = transl_rotat(data.TxCosDir, [0 0 0], az, dip);
% g.RxCosDir = transl_rotat(data.RxCosDir, [0 0 0], az, dip);
% if ~isnan(data.Tx_Z_water(1,1))
%     g.Tx_Z_water = transl_rotat(data.Tx_Z_water, origine, az, dip);
% end
% if ~isnan(data.Rx_Z_water(1,1))
%     g.Rx_Z_water = transl_rotat(data.Rx_Z_water, origine, az, dip);
% end
g.in = data.in;
% figure
% plot3(g.Tx(:,1), g.Tx(:,2), g.Tx(:,3), 'o', g.Rx(:,1), g.Rx(:,2), g.Rx(:,3), '*')
% set(gca,'DataAspectRatio',[1 1 1])

dx = str2double(get(handles.edit_dx,'String'));
dy = str2double(get(handles.edit_dy,'String'));
dz = str2double(get(handles.edit_dz,'String'));

if isnan(dx)
	xmin = min([g.Tx(g.in,1);g.Rx(g.in,1)]);
	xmax = max([g.Tx(g.in,1);g.Rx(g.in,1)]);
	nx = str2double(get(handles.edit_nx,'String'));
	dx = (xmax-xmin)/nx;
	xmin = xmin - 0.5*dx;
	xmax = xmax + 0.5*dx;
	dx = (xmax-xmin)/nx;
	set(handles.edit_dx,'String',num2str(dx))
else
	xmin = min([g.Tx(g.in,1);g.Rx(g.in,1)]) - 0.5*dx;
	xmax = max([g.Tx(g.in,1);g.Rx(g.in,1)]) + 0.5*dx;
	nx = ceil((xmax-xmin)/dx);
end
if isnan(dy)
	ymin = min([g.Tx(g.in,1);g.Rx(g.in,1)]);
	ymax = max([g.Tx(g.in,1);g.Rx(g.in,1)]);
	ny = str2double(get(handles.edit_nx,'String'));
	dy = (ymax-ymin)/ny;
	ymin = ymin - 0.5*dy;
	ymax = ymax + 0.5*dy;
	dy = (xmax-xmin)/nx;
	set(handles.edit_dy,'String',num2str(dy))
else
	ymin = min([g.Tx(g.in,1);g.Rx(g.in,1)]) - 0.5*dy;
	ymax = max([g.Tx(g.in,1);g.Rx(g.in,1)]) + 0.5*dy;
	ny = ceil((ymax-ymin)/dy);
end
if isnan(dz)
	zmin = min([g.Tx(g.in,3);g.Rx(g.in,3)]);
	zmax = max([g.Tx(g.in,3);g.Rx(g.in,3)]);
	nz = str2double(get(handles.edit_nz,'String'));
	dz = (zmax-zmin)/nz;
	zmin = zmin - 0.5*dz;
	zmax = zmax + 0.5*dz;
	dz = (zmax-zmin)/nz;
	set(handles.edit_dz,'String',num2str(dz))
else
	zmin = min([g.Tx(g.in,3);g.Rx(g.in,3)]) - 0.5*dz;
	zmax = max([g.Tx(g.in,3);g.Rx(g.in,3)]) + 0.5*dz;
	nz = ceil((zmax-zmin)/dz);
end
nxm = str2double(get(handles.edit_nxm,'String'));
nxp = str2double(get(handles.edit_nxp,'String'));
nym = str2double(get(handles.edit_nym,'String'));
nyp = str2double(get(handles.edit_nyp,'String'));
nzm = str2double(get(handles.edit_nzm,'String'));
nzp = str2double(get(handles.edit_nzp,'String'));
g.grx = xmin+dx*((0-nxm):(nx+nxp));
g.gry = ymin+dy*((0-nym):(ny+nyp));
g.grz = zmin+dz*((0-nzm):(nz+nzp));
g.grx = g.grx';
g.gry = g.gry';
g.grz = g.grz';
setappdata(handles.fig_bh_grid3d,'g',g)
%plot_proj(handles)

function plot_proj(handles)
g = getappdata(handles.fig_bh_grid3d,'g');

xlim = [min(g.grx) max(g.grx)];
zlim = [min(g.grz) max(g.grz)];
xx1 = repmat(g.grx(:).',3,1) ;
zz1 = repmat([zlim(:) ; nan],1,numel(g.grx)) ;
xx2 = repmat([xlim(:) ; nan],1,numel(g.grz)) ;
zz2 = repmat(g.grz(:).',3,1) ;
xx1 = [xx1 xx2] ;
zz1 = [zz1 zz2] ;

%aa=[min(g.grx)*ones(length(g.grz),1) max(g.grx)*ones(length(g.grz),1);g.grx g.grx];
%bb=[g.grz g.grz; min(g.grz)*ones(length(g.grx),1) max(g.grz)*ones(length(g.grx),1)];

axes(handles.axes_proj);
%plot(aa',bb','Color',[0.5 0.5 0.5])
%plot3(Tx(:,1),Tx(:,2),Tx(:,3),'b.')
plot(xx1, zz1,'Color',[0.5 0.5 0.5])
hold on
plot(g.Tx(g.in,1),g.Tx(g.in,3),'b.')
%plot3(Rx(:,1),Rx(:,2),Rx(:,3),'g.')
plot(g.Rx(g.in,1),g.Rx(g.in,3),'g.')

hold off
axis equal
axis tight
xlabel('X')
zlabel('Z')
%ylabel('Y')
%set(handles.axes_proj,'ZDir','reverse')

function edit_x0_Callback(hObject, eventdata, handles)
val = str2double(get(hObject,'string'));
if isnan(val)
    errordlg(str.s54,str.s45,'modal')
    return
end
update_proj(handles)
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
update_proj(handles)
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
update_proj(handles)
update_info(handles)


function edit_z0_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


function [az,dip] = calcul_azimuth_dip(handles)
h = getappdata(handles.fig_bh_grid3d,'h');
% az p/r a l'axe des x
% dip est l'angle p/r a l'horizontal


%
% forme normale de l'équation du plan
%
%   x*a(1) + y*a(2) + z*a(3) = d
%
%   avec d la distance de l'origine au plan
%

d = sum( h.x0 .* h.a );
x = d/h.a(1);
y = d/h.a(2);
az = atan2(y,x);
%a=h.a;
% rotation p/r azimuth
%if abs(az) > (pi/720)  % si plus grand que 1/4 degre
%    rot = [cos(az) -sin(az); sin(az) cos(az)];
%    a(1:2) = a(1:2)*rot';
%end
dip = asin(h.a(3));
%if a(2)<0; dip=-dip; end
g = getappdata(handles.fig_bh_grid3d,'g');
az = az + g.flip*pi;


function update_info(handles)
str = getappdata(handles.fig_bh_grid3d,'str');
g = getappdata(handles.fig_bh_grid3d,'g');
data = getappdata(handles.fig_bh_grid3d,'data');

nx = length(g.grx);
ny = length(g.gry);
nz = length(g.grz);
texte = {'';
    [str.s34,': ',num2str(nx*nz)];
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
no = get(handles.popupmenu_origine,'Value');
tmp = 1:2;
no2 = tmp(no~=tmp);

f1(1) = data.boreholes(no).X;
f1(2) = data.boreholes(no).Y;
f1(3) = data.boreholes(no).Z_surf;
f2(1) = data.boreholes(no2).X;
f2(2) = data.boreholes(no2).Y;
f2(3) = data.boreholes(no2).Z_surf;

f1 = proj_plan(f1, h.x0, h.a);
f2 = proj_plan(f2, h.x0, h.a);
origine(1) = str2double(get(handles.edit_x0,'String'));
origine(2) = str2double(get(handles.edit_y0,'String'));
origine(3) = str2double(get(handles.edit_z0,'String'));
[az,dip] = calcul_azimuth_dip(handles);

f1 = transl_rotat(f1, origine, az, dip);
f2 = transl_rotat(f2, origine, az, dip);

tmp = [f1(1) f1(3) f2(1) f2(3)]; %  [xTx zTx xRx zRx]
for n=[no no2]
	if isfield( data.boreholes(n), 'scont' )
		g.cont.slowness = ajouteContraintesBH(g.cont.slowness, ...
			data.boreholes(n).scont, h.x0, h.a, origine, az, dip);
	end
	if isfield( data.boreholes(n), 'acont' )
		g.cont.attenuation = ajouteContraintesBH(g.cont.attenuation, ...
			data.boreholes(n).acont, h.x0, h.a, origine, az, dip);
	end
end

plan.x0      = h.x0;
plan.a       = h.a;
plan.origine = origine;
plan.az      = az;
plan.dip     = dip;
	
[hh, g.cont] = bh_tomo_contraintes3d(grx, gry, grz, tmp, g.cont, plan);
if ishandle(hh), close(hh); end
setappdata(handles.fig_bh_grid3d,'g',g)


function pushbutton_quitter_Callback(hObject, eventdata, handles)
guidata(hObject, handles);
uiresume;


function fig_bh_grid3d_CloseRequestFcn(hObject, eventdata, handles)
guidata(hObject, handles);
uiresume(handles.fig_bh_grid3d);
%delete(hObject);



function pushbutton_qualite_Callback(hObject, eventdata, handles)
h = getappdata(handles.fig_bh_grid3d,'h');
g = getappdata(handles.fig_bh_grid3d,'g');
data = getappdata(handles.fig_bh_grid3d,'data');

figure

% distance entre les pts originaux et les pts projetÃ©s dans le plan

dTx = sqrt( sum( (data.Tx-h.Tx).^2, 2) );

subplot(321)
plot(dTx,'o')
title('Distance between original and projected Tx')

dRx = sqrt( sum( (data.Rx-h.Rx).^2, 2) );

subplot(322)
plot(dRx,'o')
title('Distance between original and projected Rx')

l_orig = sqrt( sum( (data.Tx-data.Rx).^2,2 ) );
l_new = sqrt( sum( (h.Tx-h.Rx).^2,2 ) );

% Cosinus directeurs des antennes après rotation

subplot(323)
plot(g.TxCosDir(:,1),'b')
hold on
plot(g.TxCosDir(:,2),'r')
plot(g.TxCosDir(:,3),'g')
hold off
title('Tx direction cosines after rotation')

subplot(324)
plot(g.RxCosDir(:,1),'b')
hold on
plot(g.RxCosDir(:,2),'r')
plot(g.RxCosDir(:,3),'g')
hold off
title('Rx direction cosines after rotation')


% Erreur relative sur la longueur des rais aprÃ¨s projection
diff = 100*(l_orig - l_new)./l_orig;
subplot(325)
plot(diff)
str = getappdata(handles.fig_bh_grid3d,'str');
title(str.s196)


function edit_nx_Callback(hObject, eventdata, handles)
nx = str2double(get(hObject,'String'));
if isnan(nx)
    str = getappdata(handles.fig_bh_grid3d,'str');
    errordlg(str.s54,str.s45,'modal')
    return
end
g = getappdata(handles.fig_bh_grid3d,'g');
xmin = min([g.Tx(g.in,1);g.Rx(g.in,1)]);
xmax = max([g.Tx(g.in,1);g.Rx(g.in,1)]);
dx = (xmax-xmin)/nx;
xmin = xmin-0.5*dx;
xmax = xmax+0.5*dx;
dx = (xmax-xmin)/nx;
set(handles.edit_dx,'String',num2str(dx))

nxm = str2double(get(handles.edit_nxm,'String'));
nxp = str2double(get(handles.edit_nxp,'String'));
g.grx = xmin+dx*((0-nxm):(nx+nxp));
g.grx=g.grx';
setappdata(handles.fig_bh_grid3d,'g',g)
%plot_proj(handles)
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
ymin = min([g.Tx(g.in,1);g.Rx(g.in,1)]);
ymax = max([g.Tx(g.in,1);g.Rx(g.in,1)]);
dy = (ymax-ymin)/ny;
ymin = ymin-0.5*dy;
ymax = ymax+0.5*dy;
dy = (ymax-ymin)/ny;
set(handles.edit_dy,'String',num2str(dy))

nym = str2double(get(handles.edit_nym,'String'));
nyp = str2double(get(handles.edit_nyp,'String'));
g.gry = ymin+dy*((0-nym):(ny+nyp));
g.gry=g.gry';
setappdata(handles.fig_bh_grid3d,'g',g)
%plot_proj(handles)
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
zmin = min([g.Tx(g.in,3);g.Rx(g.in,3)]);
zmax = max([g.Tx(g.in,3);g.Rx(g.in,3)]);
dz = (zmax-zmin)/nz;
zmin = zmin - 0.5*dz;
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
