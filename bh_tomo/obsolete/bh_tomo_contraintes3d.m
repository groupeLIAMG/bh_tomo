function varargout = bh_tomo_contraintes3d(varargin)
% BH_TOMO_CONTRAINTES3D M-file for bh_tomo_contraintes3d.fig
%      BH_TOMO_CONTRAINTES3D, by itself, creates a new BH_TOMO_CONTRAINTES3D or raises the existing
%      singleton*.
%
%      H = BH_TOMO_CONTRAINTES3D returns the handle to a new BH_TOMO_CONTRAINTES3D or the handle to
%      the existing singleton*.
%
%      BH_TOMO_CONTRAINTES3D('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in BH_TOMO_CONTRAINTES3D.M with the given input arguments.
%
%      BH_TOMO_CONTRAINTES3D('Property','Value',...) creates a new BH_TOMO_CONTRAINTES3D or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before bh_tomo_contraintes_OpeningFunction gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to bh_tomo_contraintes3d_OpeningFcn via varargin.
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

% Edit the above text to modify the response to help bh_tomo_contraintes3d

% Last Modified by GUIDE v2.5 07-Feb-2013 14:03:56

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
    'gui_Singleton',  gui_Singleton, ...
    'gui_OpeningFcn', @bh_tomo_contraintes3d_OpeningFcn, ...
    'gui_OutputFcn',  @bh_tomo_contraintes3d_OutputFcn, ...
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


% --- Executes just before bh_tomo_contraintes3d is made visible.
function bh_tomo_contraintes3d_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to bh_tomo_contraintes3d (see VARARGIN)

% Choose default command line output for bh_tomo_contraintes3d
handles.output = hObject;
handles.second_output = [];

% Update handles structure
guidata(hObject, handles);

h.grx = [-0.2 0.2 1];
h.gry = [-0.2 0.2 1];
h.grz = [-0.2 0.2 2];
h.n_vitesse = 0.1;
h.n_att = 1;
h.v_air = 0.2998;    % vitesse ds air
h.a_air = 0;      % attenuation ds air
h.f = [];
h.v_couche_inf = 0.12;
h.a_couche_inf = 0.5;
h.change_plot = true;
h.cont_orig = [];

if nargin>=7
    h.cont_orig = varargin{4};
end
if nargin>=6
    h.grz = varargin{3};
end
if nargin>=5
    h.gry = varargin{2};
end
if nargin>=4
    h.grx = varargin{1};
end
str = get_str_locale();

x=(h.grx(1):h.grx(2):h.grx(3))';
y=(h.gry(1):h.gry(2):h.gry(3))';
z=(h.grz(1):h.grz(2):h.grz(3))';

ind1=2:length(x);
ind2=1:length(x)-1;
h.gridx=(x(ind2)+x(ind1))/2;
ind1=2:length(y);
ind2=1:length(y)-1;
h.gridy=(y(ind2)+y(ind1))/2;
ind1=2:length(z);
ind2=1:length(z)-1;
h.gridz=(z(ind2)+z(ind1))/2;

set(handles.popupmenu_y,'String',h.gridy);

h.vitesse = nan(length(h.gridx), length(h.gridy),length(h.gridz));

DX = h.grx(2);
DY = h.gry(2);
DZ = h.grz(2);
h.xx = (h.grx(1)+DX/2):DX:(h.grx(3)-DX/2);
h.yy = (h.gry(1)+DY/2):DY:(h.gry(3)-DY/2);
h.zz = (h.grz(1)+DZ/2):DZ:(h.grz(3)-DZ/2);
h.SW = nan(length(h.xx),length(h.yy),length(h.zz));
h.XI = nan(length(h.xx),length(h.yy),length(h.zz));
h.ATT = h.SW;
h.VAR_S = h.SW;
h.VAR_xi = h.SW;
h.VAR_A = h.SW;
h.variance_s = 1;
h.variance_a = 0.25;

[h.XX,h.YY,h.ZZ] = meshgrid(h.xx,h.yy,h.zz);

h.aa = [min(x)*ones(length(z),1) max(x)*ones(length(z),1);x x]';
h.bb = [z z; min(z)*ones(length(x),1) max(z)*ones(length(x),1)]';

if ~isempty( h.cont_orig.slowness )
    s = h.cont_orig.slowness.data;
    for n=1:size(s,1)
        ix = findnear(s(n,1), h.xx);
        iy = findnear(s(n,2), h.yy);
        iz = findnear(s(n,3), h.zz);
        h.SW(ix,iy,iz) = s(n,4);
        if size(s,2)==5
            h.VAR_S(ix,iy,iz) = s(n,5);
        else
            h.VAR_S(ix,iy,iz) = 0;
        end
    end
    if isfield( h.cont_orig.slowness, 'data_xi' )
        s = h.cont_orig.slowness.data_xi;
        for n=1:size(s,1)
            ix = findnear(s(n,1), h.xx);
            iy = findnear(s(n,2), h.yy);
            iz = findnear(s(n,3), h.zz);
            h.XI(ix,iy,iz) = s(n,4);
            if size(s,2)==5
                h.VAR_xi(ix,iy,iz) = s(n,5);
            else
                h.VAR_xi(ix,iy,iz) = 0;
            end
        end
    end
end

if ~isempty( h.cont_orig.attenuation )
    s = h.cont_orig.attenuation.data;
    for n=1:size(s,1)
        ix = findnear(s(n,1), h.xx);
        iy = findnear(s(n,2), h.yy);
        iz = findnear(s(n,3), h.zz);
        h.ATT(ix,iy,iz) = s(n,4);
        if size(s,2)==5
            h.VAR_A(ix,iy,iz) = s(n,5);
        else
            h.VAR_A(ix,iy,iz) = 0;
        end
    end
    if isfield( h.cont_orig.slowness, 'variance' )
        h.variance_s = h.cont_orig.slowness.variance;
    end
    if isfield( h.cont_orig.attenuation, 'variance' )
        h.variance_a = h.cont_orig.attenuation.variance;
    end
end

for i = 1:length(h.SW(1,1,:))
    h.SW(:,:,i) = reshape(h.SW(:,:,i)',length(h.SW(:,1,1)),length(h.SW(1,:,1)));
    h.VAR_S(:,:,i) = reshape(h.VAR_S(:,:,i)',length(h.SW(:,1,1)),length(h.SW(1,:,1)));
    h.ATT(:,:,i) = reshape(h.ATT(:,:,i)',length(h.SW(:,1,1)),length(h.SW(1,:,1)));
    h.VAR_A(:,:,i) = reshape(h.VAR_A(:,:,i)',length(h.SW(:,1,1)),length(h.SW(1,:,1)));
end
set(handles.edit_variance_cont,'String',num2str(h.variance_s));

setappdata(handles.fig_bh_cont3d, 'h', h);
setappdata(handles.fig_bh_cont3d, 'str', str);
set_str_locale(handles);
update_fig(handles);

% UIWAIT makes bh_tomo_contraintes3d wait for user response (see UIRESUME)
uiwait(handles.fig_bh_cont3d);

% --- Outputs from this function are returned to the command line.
function varargout = bh_tomo_contraintes3d_OutputFcn(hObject, eventdata, handles)
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure

%varargout{1} = handles.output;
varargout{1} = hObject;
varargout{2} = handles.second_output;
delete(hObject);



function edit_valeur_Callback(hObject, eventdata, handles)
h = getappdata(handles.fig_bh_cont3d, 'h');
h.n_vitesse = str2double(get(hObject,'String'));
setappdata(handles.fig_bh_cont3d, 'h', h)


function edit_valeur_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


function pushbutton_edit_Callback(hObject, eventdata, handles)
if get(handles.radiobutton_vel,'Value')==1
	edit_SW(handles)
else
    edit_ATT(handles)
end


function pushbutton_save_Callback(hObject, eventdata, handles)
cont = prepare_cont(handles);
[fichier,rep]=uiputfile('*.dat','Fichier contraintes');
fid=fopen([rep,fichier],'wt');
if get(handles.radiobutton_vel,'Value')==1
    fprintf(fid,'%f     %f     %f     %d\n', cont.slowness.data);
else
    fprintf(fid,'%f     %f     %f     %d\n', cont.attenuation.data);
end
fclose(fid);


function pushbutton_quit_Callback(hObject, eventdata, handles)
handles.second_output = prepare_cont(handles);
guidata(hObject, handles);
uiresume;%(handles.fig_bh_cont3d);

function fig_bh_cont3d_CloseRequestFcn(hObject, eventdata, handles)
handles.second_output = prepare_cont(handles);
guidata(hObject, handles);
uiresume;
%delete(hObject);



function togglebutton_zoom_Callback(hObject, eventdata, handles)
button_state = get(hObject,'Value');
if button_state == get(hObject,'Max')
    zoom(handles.fig_bh_cont3d,'on')
elseif button_state == get(hObject,'Min')
    zoom(handles.fig_bh_cont3d,'off')
end

function pushbutton_zReset_Callback(hObject, eventdata, handles)
zoom(handles.fig_bh_cont3d,'out')


function update_fig(handles)
%h = getappdata(handles.fig_bh_cont3d, 'h');
if get(handles.radiobutton_vel,'Value')==1
    update_fig_vitesse(handles)
else
    update_fig_att(handles)
end


function update_fig_vitesse(handles)
h = getappdata(handles.fig_bh_cont3d, 'h');
iy = get(handles.popupmenu_y,'Value');
SW = reshape(h.SW(:,iy,:),length(h.SW(:,1,1)),length(h.SW(1,1,:)))';
if h.change_plot
    str = getappdata(handles.fig_bh_cont3d, 'str');
    axes(handles.axes1)
    edit_cmax_Callback(handles);
    hold(handles.axes1,'on')
    plot(handles.axes1,h.aa,h.bb,'Color',[0.5 0.5 0.5])
    xlabel(str.s119)
    ylabel(str.s120)
    set(handles.axes1,'YDir','normal')
    hold(handles.axes1,'off')
    h.change_plot = false;
else
    set(h.h1,'CData',1./SW)
end
setappdata(handles.fig_bh_cont3d, 'h', h)


function update_fig_att(handles)
h = getappdata(handles.fig_bh_cont3d, 'h');
iy = get(handles.popupmenu_y,'Value');
ATT = reshape(h.ATT(:,iy,:),length(h.ATT(:,1,1)),length(h.ATT(1,1,:)))';
if h.change_plot
    str = getappdata(handles.fig_bh_cont3d, 'str');
    axes(handles.axes1)
    h.h1 = imagesc(h.xx,h.zz,ATT);
    %caxis(handles.axes1,[0 1]);
    edit_cmax_Callback(handles)
    %colorbar('peer',handles.axes1);
    hold(handles.axes1,'on')
    plot(handles.axes1,h.aa,h.bb,'Color',[0.5 0.5 0.5])
    xlabel(str.s119)
    ylabel(str.s120)
    set(handles.axes1,'YDir','normal');
    hold(handles.axes1,'off');
    h.change_plot = false;
else
    set(h.h1,'CData',ATT);
end
setappdata(handles.fig_bh_cont3d, 'h', h)

% function edit_vitesse(handles)
% h = getappdata(handles.fig_bh_cont3d, 'h');
% axes(handles.axes1);
% [x,z,b] = ginput(1);
% get(handles.pushbutton_edit,'Value')
% iy = get(handles.popupmenu_y,'Value');
% while b==1
%     ix=findnear(x,h.gridx);
%     iz=findnear(z,h.gridz);
%     h.vitesse(ix,iy,iz) = h.n_vitesse;
%     setappdata(handles.fig_bh_cont3d, 'h', h);
%     update_fig(handles);
%     [x,z,b] = ginput(1);
% end

function edit_SW(handles)
h = getappdata(handles.fig_bh_cont3d, 'h');
axes(handles.axes1);
[x,z,b] = ginput(1);
%get(handles.pushbutton_edit,'Value')
iy = get(handles.popupmenu_y,'Value');
while b==1
    ix=findnear(x,h.xx);
    iz=findnear(z,h.zz);
    h.SW(ix,iy,iz) = 1/h.n_vitesse;
    h.VAR_S(ix,iy,iz) = h.variance_s;
    setappdata(handles.fig_bh_cont3d, 'h', h);
    update_fig(handles);
    [x,z,b] = ginput(1);
end

function edit_ATT(handles)
h = getappdata(handles.fig_bh_cont3d, 'h');
axes(handles.axes1);
[x,z,b] = ginput(1);
get(handles.pushbutton_edit,'Value')
iy = get(handles.popupmenu_y,'Value');
while b==1
    ix=findnear(x,h.xx);
    iz=findnear(z,h.zz);
    h.ATT(ix,iy,iz) = h.n_vitesse;
    h.VAR_A(ix,iy,iz) = h.variance_a;
    setappdata(handles.fig_bh_cont3d, 'h', h);
    update_fig(handles);
    [x,z,b] = ginput(1);
end

function cont = prepare_cont(handles)
h = getappdata(handles.fig_bh_cont3d, 'h');

ind = find( ~isnan(h.SW) );
slown = h.SW(ind);
ZZ = h.ZZ(ind);
YY = h.YY(ind);
XX = h.XX(ind);
var = h.VAR_S(ind);
cont.slowness.data = [XX YY ZZ slown var];
cont.slowness.variance = h.variance_s;

ind = find( ~isnan(h.ATT) );
att = h.ATT(ind);
ZZ = h.ZZ(ind);
YY = h.YY(ind);
XX = h.XX(ind);
var = h.VAR_A(ind);
cont.attenuation.data = [XX YY ZZ att var];
cont.attenuation.variance = h.variance_a;


function pushbutton_annuler_Callback(hObject, eventdata, handles)
h = getappdata(handles.fig_bh_cont3d, 'h');
handles.second_output = h.cont_orig;
guidata(hObject, handles);
uiresume(handles.fig_bh_cont3d);


function radiobutton_vel_Callback(hObject, eventdata, handles)
if get(hObject,'Value')==1
    h = getappdata(handles.fig_bh_cont3d, 'h');
    str = getappdata(handles.fig_bh_cont3d, 'str');
    set(handles.text_valeur,'String',[str.s121,' [m/ns]'])
    set(handles.edit_valeur,'String',num2str(h.n_vitesse))
    set(handles.edit_variance_cont,'String',num2str(h.variance_s))
    h.change_plot = true;
    setappdata(handles.fig_bh_cont3d, 'h', h)
end
update_fig(handles)

function radiobutton_att_Callback(hObject, eventdata, handles)
if get(hObject,'Value')==1
    h = getappdata(handles.fig_bh_cont3d, 'h');
    str = getappdata(handles.fig_bh_cont3d, 'str');
    set(handles.text_valeur,'String',[str.s177,' [Np/m]'])
    set(handles.edit_valeur,'String',num2str(h.n_att))
    set(handles.edit_variance_cont,'String',num2str(h.variance_a))
    h.change_plot = true;
    setappdata(handles.fig_bh_cont3d, 'h', h)
end
update_fig(handles)


function pushbutton_importer_Callback(hObject, eventdata, handles)
str = getappdata(handles.fig_bh_cont3d, 'str');
[file, rep, index] = uigetfile('*.con',str.s252);
if file==0 | index==0
    return
end
h = getappdata(handles.fig_bh_cont3d, 'h');
cont = load([rep,file]);
if get(handles.radiobutton_vel,'Value')==1
    for n=1:size(cont,1)
        ix = findnear(cont(n,1), h.xx);
        iy = findnear(cont(n,2), h.yy);
        iz = findnear(cont(n,3), h.zz);
        h.SW(ix,iy,iz) = 1/cont(n,4);
        if size(cont,2)==5
            h.VAR_S(ix,iy,iz) = cont(n,5);
        else
            h.VAR_S(ix,iy,iz) = 0;
        end
    end
else
    for n=1:size(cont,1)
        ix = findnear(cont(n,1), h.xx);
        iy = findnear(cont(n,2), h.yy);
        iz = findnear(cont(n,3), h.zz);
        h.ATT(ix,iy,iz) = cont(n,4);
        if size(cont,2)==5
            h.VAR_A(ix,iy,iz) = cont(n,5);
        else
            h.VAR_A(ix,iy,iz) = 0;
        end
    end
end
h.change_plot = true;
setappdata(handles.fig_bh_cont3d, 'h', h)

update_fig(handles)

function pushbutton_var_cont_Callback(hObject, eventdata, handles)
h = getappdata(handles.fig_bh_cont3d, 'h');
str = getappdata(handles.fig_bh_cont3d, 'str');
iy = get(handles.popupmenu_y,'Value');
if get(handles.radiobutton_vel,'Value')==1
    variance = h.variance_s;
    %VAR = h.VAR_S;
    VAR = reshape(h.VAR_S(:,iy,:),length(h.VAR_S(:,1,1)),length(h.VAR_S(1,1,:)))';
else
    variance = h.variance_a;
    %VAR = h.VAR_A;
    VAR = reshape(h.VAR_A(:,iy,:),length(h.VAR_A(:,1,1)),length(h.VAR_A(1,1,:)))';
end
axes(handles.axes1)
h.h1 = imagesc(h.xx,h.zz,VAR);
edit_cmax_Callback(handles);
hold(handles.axes1,'on')
plot(handles.axes1,h.aa,h.bb,'Color',[0.5 0.5 0.5])
xlabel(str.s119)
ylabel(str.s120)
set(handles.axes1,'YDir','normal')
hold(handles.axes1,'off')

[x,z,b] = ginput(1);
while b==1
    ix=findnear(x,h.xx);
    iz=findnear(z,h.zz);
    if VAR(iz,ix) == 0
        VAR(iz,ix) = variance;
    else
        VAR(iz,ix) = 0;
    end
    set(h.h1,'CData',VAR);
    [x,z,b] = ginput(1);
end
if get(handles.radiobutton_vel,'Value')==1
    h.VAR_S(:,iy,:) = VAR';
else
    h.VAR_A(:,iy,:) = VAR';
end

h.change_plot = true;
setappdata(handles.fig_bh_cont3d, 'h', h)

update_fig(handles)

function edit_variance_cont_Callback(hObject, eventdata, handles)
h = getappdata(handles.fig_bh_cont3d, 'h');
if get(handles.radiobutton_vel,'Value')==1
    h.variance_s = str2double(get(hObject,'String'));
else
    h.variance_a = str2double(get(hObject,'String'));
end
setappdata(handles.fig_bh_cont3d, 'h', h)


function edit_variance_cont_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit_cmin_Callback(hObject, eventdata, handles)
cmin = str2double(get(hObject,'String'));
cmax = str2double(get(handles.edit_cmax,'String'));
if cmax<cmin
    cmax = cmin+0.01;
    set(handles.edit_cmax, 'String', num2str(cmax))
end
caxis(handles.axes1,[cmin cmax]);
colorbar('peer',handles.axes1);


function edit_cmin_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit_cmax_Callback(handles)
cmax = str2double(get(handles.edit_cmax,'String'));
cmin = str2double(get(handles.edit_cmin,'String'));
if cmin>cmax
    cmin = cmax-0.01;
    set(handles.edit_cmin, 'String', num2str(cmin))
end
caxis(handles.axes1,[cmin cmax]);
colorbar('peer',handles.axes1);


function edit_cmax_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function pushbutton_reinit_Callback(hObject, eventdata, handles)
h = getappdata(handles.fig_bh_cont3d, 'h');
if get(handles.radiobutton_vel,'Value')==1
	h.SW(:) = nan;
else
    h.ATT(:) = nan;
end
h.change_plot = true;
setappdata(handles.fig_bh_cont3d, 'h', h)

update_fig(handles)

function popupmenu_y_Callback(hObject, eventdata, handles)
h = getappdata(handles.fig_bh_cont3d, 'h');
h.change_plot = true;
setappdata(handles.fig_bh_cont3d, 'h', h)
update_fig(handles);

function popupmenu_y_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function set_str_locale(handles)
str = getappdata(handles.fig_bh_cont3d, 'str');
h = getappdata(handles.fig_bh_cont3d, 'h');

set(handles.text_valeur,'String',[str.s121,' [m/ns]'])
set(handles.edit_valeur,'String',num2str(h.n_vitesse))
set(handles.pushbutton_edit,'String',str.s122)
set(handles.pushbutton_save,'String',str.s29)
set(handles.pushbutton_quit,'String',str.s193)
set(handles.pushbutton_annuler,'String',str.s91)
set(handles.pushbutton_importer,'String',[str.s215,' ...'])
set(handles.uipanel_prop,'Title',str.s222)
set(handles.radiobutton_vel,'String',str.s121)
set(handles.radiobutton_att,'String',str.s177)
set(handles.uipanel_prop,'Title',str.s298)
set(handles.radiobutton_vel,'String',str.s299)
set(handles.pushbutton_var_cont,'String',str.s302)
set(handles.text_variance_cont,'String',str.s303)
