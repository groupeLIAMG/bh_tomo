function varargout = bh_tomo_Sgrille(varargin)
% BH_TOMO_SGRILLE M-file for bh_tomo_Sgrille.fig
%      BH_TOMO_SGRILLE, by itself, creates a new BH_TOMO_SGRILLE or raises the existing
%      singleton*.
%
%      H = BH_TOMO_SGRILLE returns the handle to a new BH_TOMO_SGRILLE or the handle to
%      the existing singleton*.
%
%      BH_TOMO_SGRILLE('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in BH_TOMO_SGRILLE.M with the given input arguments.
%
%      BH_TOMO_SGRILLE('Property','Value',...) creates a new BH_TOMO_SGRILLE or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before bh_tomo_Sgrille_OpeningFunction gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to bh_tomo_Sgrille_OpeningFcn via varargin.
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

% Edit the above text to modify the response to help bh_tomo_Sgrille

% Last Modified by GUIDE v2.5 21-Dec-2012 15:28:22

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @bh_tomo_Sgrille_OpeningFcn, ...
                   'gui_OutputFcn',  @bh_tomo_Sgrille_OutputFcn, ...
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


% --- Executes just before bh_tomo_Sgrille is made visible.
function bh_tomo_Sgrille_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to bh_tomo_Sgrille (see VARARGIN)

% Choose default command line output for bh_tomo_Sgrille
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

data = [];
g.grx = [];
g.grz = [];
g.cont.slowness.data = [];
g.cont.attenuation.data = [];
g.Tx = [];
g.Rx = [];
g.x0 = [];
g.borehole_x0 = [];
g.bord = [1 1 1 1];  % [nxm nxp nzm nzp]
if nargin >= 4
    data = varargin{1};
end
if nargin >= 5
    g = varargin{2};
end
setappdata(handles.fig_bh_grille,'data',data)
setappdata(handles.fig_bh_grille,'g',g)
update_plot3D(handles)

set(handles.edit_nxm,'String',num2str(g.bord(1)))
set(handles.edit_nxp,'String',num2str(g.bord(2)))
set(handles.edit_nzm,'String',num2str(g.bord(3)))
set(handles.edit_nzp,'String',num2str(g.bord(4)))
if isempty(g.x0)
    set(handles.edit_dx,'String','0.2')
    set(handles.edit_dz,'String','0.2')
else
    set(handles.edit_dx,'String',num2str(g.grx(2)-g.grx(1)))
    set(handles.edit_dz,'String',num2str(g.grz(2)-g.grz(1)))
end

str = get_str_locale();
setappdata(handles.fig_bh_grille,'str',str)
set_String_locale(handles, str)

if ~isempty(data)
    h.order_boreholes = boreholes_order( data.boreholes );
    h.n_plans = length(data.boreholes)-1;
    h.ordre_plans = 1:h.n_plans;
    f1 = h.order_boreholes(1);
    f2 = h.order_boreholes( length(h.order_boreholes) );
    tmp{1} = data.boreholes(f1).name;
    tmp{2} = data.boreholes(f2).name;
    set(handles.popupmenu_origine,'String',tmp)
    if ~isempty(g.borehole_x0)
        set(handles.popupmenu_origine,'Value',g.borehole_x0)
        if g.borehole_x0~=1
            h.order_boreholes = fliplr(h.order_boreholes);
            h.ordre_plans = fliplr(h.ordre_plans);
        end
    end
    setappdata(handles.fig_bh_grille,'h',h)
    if isempty(g.x0)
        set(handles.edit_x0,'String',num2str(data.boreholes(f1).X))
        set(handles.edit_y0,'String',num2str(data.boreholes(f1).Y))
        set(handles.edit_z0,'String',num2str(data.boreholes(f1).Z))
    else
        set(handles.edit_x0,'String',num2str(g.x0(1)))
        set(handles.edit_y0,'String',num2str(g.x0(2)))
        set(handles.edit_z0,'String',num2str(g.x0(3)))
    end
    construit_grd(handles)
    if isempty(g.grx)
        update_proj(handles)
    else
        plot_proj(handles)
    end
    update_info(handles)
end

% UIWAIT makes bh_tomo_Sgrille wait for user response (see UIRESUME)
uiwait(handles.fig_bh_grille);

% --- Outputs from this function are returned to the command line.
function varargout = bh_tomo_Sgrille_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = hObject;
if isfield(handles, 'fig_bh_grille')
    g = getappdata(handles.fig_bh_grille,'g');
    varargout{2} = g;
else
    varargout{2} = [];
end

function construit_grd(handles)
data = getappdata(handles.fig_bh_grille,'data');
if isempty(data)
    return
end
h = getappdata(handles.fig_bh_grille,'h');

for n=1:h.n_plans
    n1 = h.order_boreholes(n);
    n2 = h.order_boreholes(n+1);
    np = h.ordre_plans(n);
    X = [data.boreholes(n1).X   data.boreholes(n1).Y    data.boreholes(n1).Z;
        data.boreholes(n1).Xmax data.boreholes(n1).Ymax data.boreholes(n1).Zmax;
        data.boreholes(n2).X    data.boreholes(n2).Y    data.boreholes(n2).Z;
        data.boreholes(n2).Xmax data.boreholes(n2).Ymax data.boreholes(n2).Zmax];
    % trouve le plan qui passe par les 4 points
    [h.plans(np).x0, h.plans(np).a]=lsplane(X);
    h.plans(np).x0 = h.plans(np).x0';  % Centroid of the data = point on the best-fit plane
    h.plans(np).a = h.plans(np).a';    % Direction cosines of the normal to the best-fit plane
    if h.plans(np).a(3)<0, h.plans(np).a=-h.plans(np).a; end
    % distance horizontale entre les trous
    h.plans(np).l = sqrt( (data.boreholes(n2).X-data.boreholes(n1).X)^2 +...
        (data.boreholes(n2).Y-data.boreholes(n1).Y)^2 );
end
% for n=2:h.n_plans
%    theta = acos( dot( h.plans(n).a, h.plans(1).a) / ...
%        (norm(h.plans(n).a)*norm(h.plans(1).a)) );
%    if theta>pi/2
%        h.plans(n).a = -1*h.plans(n).a;
%    end
% end


[h.Tx, h.Tx_no_plan] = proj_plans(data.Tx, h.plans);
[h.Rx, h.Rx_no_plan] = proj_plans(data.Rx, h.plans);
h.Tx = data.Tx;
h.Rx = data.Rx;

setappdata(handles.fig_bh_grille,'h',h)


function update_plot3D(handles)
data = getappdata(handles.fig_bh_grille,'data');
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
plot3(data.Tx(:,1), data.Tx(:,2), data.Tx(:,3),'b*')
plot3(data.Rx(:,1), data.Rx(:,2), data.Rx(:,3),'gx')
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
set(handles.uipanel_grille,           'Title',  str.s92)
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
data = getappdata(handles.fig_bh_grille,'data');
g = getappdata(handles.fig_bh_grille,'g');
h = getappdata(handles.fig_bh_grille,'h');
no = get(hObject,'Value');

if g.borehole_x0 == no
    return
else
    g.borehole_x0 = no;
    h.order_boreholes = fliplr(h.order_boreholes);
    h.ordre_plans = fliplr(h.ordre_plans);
    h.plans = fliplr( h.plans );
end
g.x0(1) = data.boreholes(h.order_boreholes(1)).X;
g.x0(2) = data.boreholes(h.order_boreholes(1)).Y;
g.x0(3) = data.boreholes(h.order_boreholes(1)).Z;
set(handles.edit_x0,'String',num2str(g.x0(1)))
set(handles.edit_y0,'String',num2str(g.x0(2)))
set(handles.edit_z0,'String',num2str(g.x0(3)))
setappdata(handles.fig_bh_grille,'g',g)
setappdata(handles.fig_bh_grille,'h',h)
update_proj(handles)
update_info(handles)

function popupmenu_origine_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


function edit_nxm_Callback(hObject, eventdata, handles)
val = str2double(get(hObject,'string'));
if isnan(val)
    str = getappdata(handles.fig_bh_grille,'str');
    errordlg(str.s54,str.s45,'modal')
    return
end
g = getappdata(handles.fig_bh_grille,'g');
g.bord(1) = val;
setappdata(handles.fig_bh_grille,'g',g)
update_proj(handles)
update_info(handles)

function edit_nxm_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


function edit_nxp_Callback(hObject, eventdata, handles)
val = str2double(get(hObject,'string'));
if isnan(val)
    str = getappdata(handles.fig_bh_grille,'str');
    errordlg(str.s54,str.s45,'modal')
    return
end
g = getappdata(handles.fig_bh_grille,'g');
g.bord(2) = val;
setappdata(handles.fig_bh_grille,'g',g)
update_proj(handles)
update_info(handles)

function edit_nxp_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


function edit_dx_Callback(hObject, eventdata, handles)
val = str2double(get(hObject,'string'));
if isnan(val)
    str = getappdata(handles.fig_bh_grille,'str');
    errordlg(str.s54,str.s45,'modal')
    return
end
update_proj(handles)
update_info(handles)

function edit_dx_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


function edit_nzm_Callback(hObject, eventdata, handles)
val = str2double(get(hObject,'string'));
if isnan(val)
    str = getappdata(handles.fig_bh_grille,'str');
    errordlg(str.s54,str.s45,'modal')
    return
end
g = getappdata(handles.fig_bh_grille,'g');
g.bord(3) = val;
setappdata(handles.fig_bh_grille,'g',g)
update_proj(handles)
update_info(handles)

function edit_nzm_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


function edit_nzp_Callback(hObject, eventdata, handles)
val = str2double(get(hObject,'string'));
if isnan(val)
    str = getappdata(handles.fig_bh_grille,'str');
    errordlg(str.s54,str.s45,'modal')
    return
end
g = getappdata(handles.fig_bh_grille,'g');
g.bord(4) = val;
setappdata(handles.fig_bh_grille,'g',g)
update_proj(handles)
update_info(handles)

function edit_nzp_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


function edit_dz_Callback(hObject, eventdata, handles)
val = str2double(get(hObject,'string'));
if isnan(val)
    str = getappdata(handles.fig_bh_grille,'str');
    errordlg(str.s54,str.s45,'modal')
    return
end
update_proj(handles)
update_info(handles)

function edit_dz_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


function update_proj(handles)
h = getappdata(handles.fig_bh_grille,'h');
g = getappdata(handles.fig_bh_grille,'g');
data = getappdata(handles.fig_bh_grille,'data');

origine(3) = str2double(get(handles.edit_z0,'String'));

[az,dip] = calcul_azimuth_dip(handles);
for n=1:h.n_plans
    no_f = h.order_boreholes(n);
    no_plan = h.ordre_plans(n);
    ind = h.Tx_no_plan==no_plan;
    oo = [data.boreholes(no_f).X data.boreholes(no_f).Y 0];
    g.Tx(ind,:) = transl_rotat(h.Tx(ind,:), oo, az(no_plan), dip(no_plan));
    if n>1
        g.Tx(ind,1) = g.Tx(ind,1)+h.plans(n-1).l;
    end
    
%     figure
%     subplot(121)
%     plot(h.Tx(ind,1),h.Tx(ind,3),'ro')
%     hold on
%     
%     subplot(122)
%     plot(g.Tx(ind,1),g.Tx(ind,3),'ro')
%     hold on

    ind = h.Rx_no_plan==no_plan;
    g.Rx(ind,:) = transl_rotat(h.Rx(ind,:), oo, az(no_plan), dip(no_plan));
    if n>1
        g.Rx(ind,1) = g.Rx(ind,1)+h.plans(n-1).l;
    end
    
%     subplot(121)
%     plot(h.Rx(ind,1),h.Rx(ind,3),'x')
%     hold off
%     
%     subplot(122)
%     plot(g.Rx(ind,1),g.Rx(ind,3),'x')
%     hold off
end
g.Tx(:,3) = g.Tx(:,3);
g.Rx(:,3) = g.Rx(:,3);

% y = 0 en 2D
g.Tx(:,2) = 0;
g.Rx(:,2) = 0;

xmin = min([g.Tx(:,1);g.Rx(:,1)]);
g.Tx(:,1) = g.Tx(:,1)-xmin;
g.Rx(:,1) = g.Rx(:,1)-xmin;
xmin = 0;
xmax = max([g.Tx(:,1);g.Rx(:,1)]);


if false
    l_orig = sqrt( sum( (data.Tx-data.Rx).^2,2 ) );
    l_new = sqrt( sum( (h.Tx-h.Rx).^2,2 ) );
	diff = 100*(l_orig - l_new)./l_orig;
	figure
	plot(diff)
	str = getappdata(handles.fig_bh_grille,'str');
	title(str.s196)
end

% zf = zeros(length(data.boreholes),1);
% for n=1:length(data.boreholes)
%     zf(n) = data.boreholes(n).Z;
% end

dx = str2double(get(handles.edit_dx,'String'));
dz = str2double(get(handles.edit_dz,'String'));

zmin = min([g.Tx(:,3);g.Rx(:,3)]);
zmax = max([g.Tx(:,3);g.Rx(:,3)]);

zmin = origine(3) - ceil( (origine(3)-zmin)/dz )*dz;
zmax = origine(3) + ceil( (zmax-origine(3))/dz )*dz;

nx = ceil((xmax-xmin)/dx);
nz = ceil((zmax-zmin)/dz);
nxm = str2double(get(handles.edit_nxm,'String'));
nxp = str2double(get(handles.edit_nxp,'String'));
nzm = str2double(get(handles.edit_nzm,'String'));
nzp = str2double(get(handles.edit_nzp,'String'));
g.grx = xmin+dx*((0-nxm):(nx+nxp));
g.grz = zmin+dz*((0-nzm):(nz+nzp));
g.grx=g.grx';
g.grz=g.grz';
setappdata(handles.fig_bh_grille,'g',g)
plot_proj(handles)

function plot_proj(handles)
g = getappdata(handles.fig_bh_grille,'g');
aa=[min(g.grx)*ones(length(g.grz),1) max(g.grx)*ones(length(g.grz),1);g.grx g.grx];
bb=[g.grz g.grz; min(g.grz)*ones(length(g.grx),1) max(g.grz)*ones(length(g.grx),1)];

axes(handles.axes_proj);
%plot3(Tx(:,1),Tx(:,2),Tx(:,3),'b.')
plot(g.Tx(:,1),g.Tx(:,3),'b*')
hold on
%plot3(Rx(:,1),Rx(:,2),Rx(:,3),'g.')
plot(g.Rx(:,1),g.Rx(:,3),'gx')
plot(aa',bb','Color',[0.5 0.5 0.5])

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

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
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

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


function edit_z0_Callback(hObject, eventdata, handles)
val = str2double(get(hObject,'string'));
if isnan(val)
    errordlg(str.s54,str.s45,'modal')
    return
end
g = getappdata(handles.fig_bh_grille,'g');
g.x0(3) = val;
setappdata(handles.fig_bh_grille,'g',g)
update_proj(handles)
update_info(handles)


function edit_z0_CreateFcn(hObject, eventdata, handles)

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


function [az,dip] = calcul_azimuth_dip(handles)
% az  est l'angle p/r a l'axe des x
% dip est l'angle p/r a l'hrizontal
data = getappdata(handles.fig_bh_grille,'data');
h = getappdata(handles.fig_bh_grille,'h');

az = zeros(1,h.n_plans);
dip = zeros(1,h.n_plans);
for nn=1:h.n_plans
    no1 = h.order_boreholes(nn);
    no2 = h.order_boreholes(nn+1);
    
    % determination de l'azimuth du plan
    origine(1) = data.boreholes(no1).X;
    origine(2) = data.boreholes(no1).Y;
    origine(3) = data.boreholes(no1).Z;
    
    % eq de la droite du borehole no2
    d=sqrt( sum([data.boreholes(no2).Xmax-data.boreholes(no2).X ...
        data.boreholes(no2).Ymax-data.boreholes(no2).Y ...
        data.boreholes(no2).Zmax-data.boreholes(no2).Z].^2) );
    % cosinus directeurs
    l=(data.boreholes(no2).Xmax-data.boreholes(no2).X)/d;
    m=(data.boreholes(no2).Ymax-data.boreholes(no2).Y)/d;
    n=(data.boreholes(no2).Zmax-data.boreholes(no2).Z)/d;
    % intersection  Z = origine(3)
    x2 = data.boreholes(no2).X + l*(origine(3)-data.boreholes(no2).Z)/n;
    y2 = data.boreholes(no2).Y + m*(origine(3)-data.boreholes(no2).Z)/n;
%     az(nn) = atan2(y2-origine(2), x2-origine(1));

    no_plan = h.ordre_plans(nn);
    d = sum( h.plans(no_plan).x0 .* h.plans(no_plan).a );
    x = d/h.plans(no_plan).a(1);
    y = d/h.plans(no_plan).a(2);
    az(no_plan) = atan2(y,x);

    
    rot = [cos(az(no_plan)) -sin(az(no_plan)); sin(az(no_plan)) cos(az(no_plan))];
    xr = [x2-origine(1) y2-origine(2)]*rot';
    if ( xr(1) < 0 ) 
        az(no_plan) = az(no_plan) + pi;
    end
    
    % boreholes verticaux
    dip(nn) = 0;
end


function update_info(handles)
str = getappdata(handles.fig_bh_grille,'str');
g = getappdata(handles.fig_bh_grille,'g');
data = getappdata(handles.fig_bh_grille,'data');

nx = length(g.grx);
nz = length(g.grz);
texte = {'';
    [str.s34,': ',num2str(nx*nz)];
    [str.s130,': ',num2str(size(data.Tx,1))];
	'';
	['Xmin: ',num2str(g.grx(1)),', Xmax: ',num2str(g.grx(end))];
	['Zmin: ',num2str(g.grz(1)),', Zmax: ',num2str(g.grz(end))]
	};
set(handles.text_info,'String',texte)




function pushbutton_contraintes_Callback(hObject, eventdata, handles)
data = getappdata(handles.fig_bh_grille,'data');
g = getappdata(handles.fig_bh_grille,'g');
h = getappdata(handles.fig_bh_grille,'h');
grx = [g.grx(1) g.grx(2)-g.grx(1) g.grx(end)];
grz = [g.grz(1) g.grz(2)-g.grz(1) g.grz(end)];

tmp = [0 data.boreholes( h.order_boreholes(1) ).Z_surf];
l = 0;
for n=1:h.n_plans
    nf = h.order_boreholes(n+1);
    l = l+h.plans(n).l;
    tmp = [tmp l data.boreholes(nf).Z_surf];
end

origine(1) = str2double(get(handles.edit_x0,'String'));
origine(2) = str2double(get(handles.edit_y0,'String'));
origine(3) = str2double(get(handles.edit_z0,'String'));
[az,dip] = calcul_azimuth_dip(handles);

for n=0:h.n_plans
    nf = h.order_boreholes(n+1);
	if n==0
		no_plan = h.ordre_plans(1);
	else
		no_plan = h.ordre_plans(n);
	end
    oo = [data.boreholes(nf).X data.boreholes(nf).Y 0];
	
	if isfield( data.boreholes(nf), 'scont' )
		g.cont.slowness = ajouteContraintesBH(g.cont.slowness, ...
			data.boreholes(nf).scont, [], [], oo, az(no_plan), dip(no_plan));
	end
	if isfield( data.boreholes(nf), 'acont' )
		g.cont.attenuation = ajouteContraintesBH(g.cont.attenuation, ...
			data.boreholes(nf).acont, [], [], oo, az(no_plan), dip(no_plan));
	end
end

plan.x0      = h.x0;
plan.a       = h.a;
plan.origine = origine;
plan.az      = az;
plan.dip     = dip;
	
[hh, g.cont] = bh_tomo_contraintes(grx, grz, tmp, g.cont, plan);
if ishandle(hh), close(hh); end
setappdata(handles.fig_bh_grille,'g',g)


function pushbutton_quitter_Callback(hObject, eventdata, handles)
guidata(hObject, handles);
uiresume;
