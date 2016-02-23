function varargout = bh_tomo_decime(varargin)
% BH_TOMO_DECIME M-file for bh_tomo_decime.fig
%      BH_TOMO_DECIME, by itself, creates a new BH_TOMO_DECIME or raises the existing
%      singleton*.
%
%      H = BH_TOMO_DECIME returns the handle to a new BH_TOMO_DECIME or the handle to
%      the existing singleton*.
%
%      BH_TOMO_DECIME('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in BH_TOMO_DECIME.M with the given input arguments.
%
%      BH_TOMO_DECIME('Property','Value',...) creates a new BH_TOMO_DECIME or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before bh_tomo_decime_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to bh_tomo_decime_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help bh_tomo_decime

% Last Modified by GUIDE v2.5 19-May-2013 22:51:51

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @bh_tomo_decime_OpeningFcn, ...
                   'gui_OutputFcn',  @bh_tomo_decime_OutputFcn, ...
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


% --- Executes just before bh_tomo_decime is made visible.
function bh_tomo_decime_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to bh_tomo_decime (see VARARGIN)

% Choose default command line output for bh_tomo_decime
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

h.mog = varargin{1};
h.in = h.mog.in;
h.params = varargin{2};
h.SB = [];
setappdata(handles.fig_bh_tomo_decime,'h',h)

str = get_str_locale();
set_String_locale(handles, str)

set(handles.edit_sautTx,'String', num2str(h.params.sautTx))
set(handles.edit_sautRx,'String', num2str(h.params.sautRx))
set(handles.edit_arrondi,'String', num2str(h.params.arrondi))
set(handles.edit_seuil_SB,'String', num2str(h.params.seuil_SB))
if isfield( h.params, 'zmin' )
    set(handles.edit_zmin,'String',num2str(h.params.zmin))
    set(handles.edit_zmax,'String',num2str(h.params.zmax))
end
set(handles.checkbox_seuil_SB,'Value', h.params.use_SB)


update(handles)

% UIWAIT makes bh_tomo_decime wait for user response (see UIRESUME)
uiwait(handles.fig_bh_tomo_decime);


% --- Outputs from this function are returned to the command line.
function varargout = bh_tomo_decime_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;
if isfield(handles, 'fig_bh_tomo_decime')
    h = getappdata(handles.fig_bh_tomo_decime,'h');
    varargout{2} = h.in;
	p.sautTx = str2double(get(handles.edit_sautTx,'String'));
	p.sautRx = str2double(get(handles.edit_sautRx,'String'));
	p.arrondi = str2double(get(handles.edit_arrondi,'String'));
	p.use_SB = get(handles.checkbox_seuil_SB,'Value');
	p.seuil_SB = str2double(get(handles.edit_seuil_SB,'String'));
    p.zmin = str2double(get(handles.edit_zmin,'String'));
    p.zmax = str2double(get(handles.edit_zmax,'String'));
	varargout{3} = p;
else
    varargout{2} = [];
	varargout{3} = [];
end


function pushbuttonOK_Callback(hObject, eventdata, handles)
guidata(hObject, handles);
uiresume;


function edit_sautTx_Callback(hObject, eventdata, handles)
update(handles)


function edit_sautTx_CreateFcn(hObject, eventdata, handles)

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


function edit_sautRx_Callback(hObject, eventdata, handles)
update(handles)

function edit_sautRx_CreateFcn(hObject, eventdata, handles)

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function update(handles)
h = getappdata(handles.fig_bh_tomo_decime,'h');
sautTx = str2double(get(handles.edit_sautTx,'String'));
sautRx = str2double(get(handles.edit_sautRx,'String'));
zmin = str2double(get(handles.edit_zmin,'String'));
zmax = str2double(get(handles.edit_zmax,'String'));

thetamin = str2double(get(handles.edit_angle_min,'String'))*pi/180;
thetamax = str2double(get(handles.edit_angle_max,'String'))*pi/180;

inTx = h.mog.data.Tx_z >= zmin & h.mog.data.Tx_z <= zmax;
inRx = h.mog.data.Rx_z >= zmin & h.mog.data.Rx_z <= zmax;

if ~any(inTx) || ~any(inRx)
    zmin = min([h.mog.data.Tx_z h.mog.data.Rx_z]);
    zmax = max([h.mog.data.Tx_z h.mog.data.Rx_z]);
    set(handles.edit_zmin,'String',num2str(zmin));
    set(handles.edit_zmax,'String',num2str(zmax));
    inTx = h.mog.data.Tx_z >= zmin & h.mog.data.Tx_z <= zmax;
    inRx = h.mog.data.Rx_z >= zmin & h.mog.data.Rx_z <= zmax;
end

Tx = [h.mog.data.Tx_x; h.mog.data.Tx_y; h.mog.data.Tx_z]';
Rx = [h.mog.data.Rx_x; h.mog.data.Rx_y; h.mog.data.Rx_z]';

dr = sqrt( (h.mog.data.Tx_x-h.mog.data.Rx_x).^2 +...
    (h.mog.data.Tx_y-h.mog.data.Rx_y).^2 );
theta = atan2(h.mog.data.Rx_z-h.mog.data.Tx_z, dr);

inTheta = theta>thetamin & theta<thetamax;

arrondi = str2double(get(handles.edit_arrondi,'String'));
if arrondi > 0
	Tx = arrondi*(round(Tx/arrondi));
	Rx = arrondi*(round(Rx/arrondi));
end
uTx = unique(Tx(inTx,:), 'rows');
uRx = unique(Rx(inRx,:), 'rows');

[m,n] = max(var(uTx));
[m,n] = sort(uTx(:,n));
uTx = uTx(n,:);
[m,n] = max(var(uRx));
[m,n] = sort(uRx(:,n));
uRx = uRx(n,:);

uTx = uTx(1:(1+sautTx):end,:);
uRx = uRx(1:(1+sautRx):end,:);



inTx = false(1,h.mog.data.ntrace);
nTx = size(uTx,1);
for n1=1:length(inTx)
	inTx(n1) = ~isempty(find(all(repmat(Tx(n1,:),nTx,1)==uTx,2),1));
end

inRx = false(1,h.mog.data.ntrace);
nRx = size(uRx,1);
for n1=1:length(inRx)
	inRx(n1) = ~isempty(find(all(repmat(Rx(n1,:),nRx,1)==uRx,2),1));
%	for n2=1:size(uRx,1)
%		inRx(n1) = inRx(n1) | all(Rx(n1,:)==uRx(n2,:));
%	end
end

if isempty(h.SB)
	h.SB = calcul_SB(h.mog);
end
seuil_SB = str2double(get(handles.edit_seuil_SB,'String'));

h.in = inTx & inRx & h.SB>seuil_SB & inTheta;
setappdata(handles.fig_bh_tomo_decime,'h',h)

[x0, a]=lsplane([uTx;uRx]);
el = (pi-a(3))*180/pi;
az = atan( cos(a(2))/cos(a(1)) )*180/pi;

if isempty(h.SB)
	h.SB = calcul_SB(h.mog);
end
seuil_SB = str2double(get(handles.edit_seuil_SB,'String'));


texte{9} = '';
texte{1} = 'Infos';
texte{2} = '';
texte{3} = [num2str(size(uTx,1)), ' Tx'];
texte{4} = [num2str(size(uRx,1)), ' Rx'];
texte{5} = [num2str(round(100*sum(~(inTx & inRx))/length(h.in)),'%d'), '% removed - Tx & Rx'];
texte{6} = [num2str(round(100*sum(h.SB<=seuil_SB)/length(h.in)),'%d'), '% removed - S/M ratio'];
texte{7} = [num2str(round(100*sum(~inTheta)/length(h.in)),'%d'), '% removed - ray angle'];
texte{8} = [num2str(round(100*sum(h.in)/length(h.in)),'%d'), '% of traces kept'];

set(handles.text_info,'String',texte)

axes(handles.axes1)
plot3(uTx(:,1), uTx(:,2), uTx(:,3),'b*')
hold on
plot3(uRx(:,1), uRx(:,2), uRx(:,3),'gx')
axis equal
grid on
xlabel('X')
ylabel('Y')
zlabel('Z')
view(az,el)
legend('Tx','Rx')
hold off


function edit_arrondi_Callback(hObject, eventdata, handles)
update(handles)

function edit_arrondi_CreateFcn(hObject, eventdata, handles)

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


function togglebutton_rotate3d_Callback(hObject, eventdata, handles)
button_state = get(hObject,'Value');
if button_state == get(hObject,'Max')
    % toggle button is pressed
	rotate3d(handles.axes1,'on')
elseif button_state == get(hObject,'Min')
    % toggle button is not pressed
	rotate3d(handles.axes1,'off')
end

function togglebutton_zoom_Callback(hObject, eventdata, handles)
button_state = get(hObject,'Value');
if button_state == get(hObject,'Max')
    % toggle button is pressed
	zoom(handles.axes1,'on')
elseif button_state == get(hObject,'Min')
    % toggle button is not pressed
	zoom(handles.axes1,'off')
end


function set_String_locale(handles, str)

load rotate
set(handles.togglebutton_rotate3d,'CData',cdata)
load zoom
set(handles.togglebutton_zoom,'CData',zoomCData)
clear cdata zoomCData


function edit_seuil_SB_Callback(hObject, eventdata, handles)
if get(handles.checkbox_seuil_SB,'Value')==1
	update(handles)
end


function edit_seuil_SB_CreateFcn(hObject, eventdata, handles)

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


function checkbox_seuil_SB_Callback(hObject, eventdata, handles)
update(handles)


function SB = calcul_SB(mog)
SB = ones(size(1,mog.data.ntrace));
[m,i] = max(abs(detrend_rad(mog.data.rdata)));

largeur = 60;

i1 = i-largeur/2;
i2 = i+largeur/2;
i1(i1<1) = 1;
i2(i2>mog.data.nptsptrc) = mog.data.nptsptrc;
for n=1:mog.data.ntrace
	SB(n) = std(mog.data.rdata(i1(n):i2(n), n))/std(mog.data.rdata(1:largeur, n));
end



function edit_zmin_Callback(hObject, eventdata, handles)
update(handles)


function edit_zmin_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit_zmax_Callback(hObject, eventdata, handles)
update(handles)


function edit_zmax_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


function edit_angle_min_Callback(hObject, eventdata, handles)
update(handles)

function edit_angle_min_CreateFcn(hObject, eventdata, handles)

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


function edit_angle_max_Callback(hObject, eventdata, handles)
update(handles)

function edit_angle_max_CreateFcn(hObject, eventdata, handles)

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
