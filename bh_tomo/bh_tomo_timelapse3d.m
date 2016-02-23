function varargout = bh_tomo_timelapse3d(varargin)
% BH_TOMO_TIMELAPSE3D MATLAB code for bh_tomo_timelapse3d.fig
%      BH_TOMO_TIMELAPSE3D, by itself, creates a new BH_TOMO_TIMELAPSE3D or raises the existing
%      singleton*.
%
%      H = BH_TOMO_TIMELAPSE3D returns the handle to a new BH_TOMO_TIMELAPSE3D or the handle to
%      the existing singleton*.
%
%      BH_TOMO_TIMELAPSE3D('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in BH_TOMO_TIMELAPSE3D.M with the given input arguments.
%
%      BH_TOMO_TIMELAPSE3D('Property','Value',...) creates a new BH_TOMO_TIMELAPSE3D or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before bh_tomo_timelapse3d_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to bh_tomo_timelapse3d_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help bh_tomo_timelapse3d

% Last Modified by GUIDE v2.5 26-Feb-2013 10:10:06

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @bh_tomo_timelapse3d_OpeningFcn, ...
                   'gui_OutputFcn',  @bh_tomo_timelapse3d_OutputFcn, ...
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


% --- Executes just before bh_tomo_timelapse3d is made visible.
function bh_tomo_timelapse3d_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to bh_tomo_timelapse3d (see VARARGIN)

% Choose default command line output for bh_tomo_timelapse3d
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes bh_tomo_timelapse3d wait for user response (see UIRESUME)
% uiwait(handles.fig_bh_timelapse);

h.db_file = get(handles.fig_bh_timelapse,'UserData');
h.tt_b_coul_min = 0.06;
h.tt_b_coul_max = 0.12;

h.tt_b_coul_min_diff = -0.02;
h.tt_b_coul_max_diff = 0.02;

h.amp_b_coul_min = 0;
h.amp_b_coul_max = 3;
h.modele_initial = [];
str = get_str_locale();

h.axes_s1 = subplot(4,3,[4 7 10]);
h.axes_s2 = subplot(4,3,[5 8 11]);
h.axes_ds = subplot(4,3,[6 9 12]);

setappdata(handles.fig_bh_timelapse, 'h', h);
setappdata(handles.fig_bh_timelapse, 'str', str);
set_string_locale(handles);

% --- Outputs from this function are returned to the command line.
function varargout = bh_tomo_timelapse3d_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


function popupmenu_model3d_Callback(hObject, eventdata, handles)
h = getappdata(handles.fig_bh_timelapse,'h');
model3d = h.model3d;
str_model3d = get(handles.popupmenu_model3d,'String');
no_model3d = get(handles.popupmenu_model3d,'Value');
if isempty(model3d(no_model3d).inv_res)
    msgbox(['there is no inversion result for ' str_model3d(no_model3d) '.'],'warning','modal');
    set(handles.popupmenu_model3d,'Value',1);
    return;
end

names_tomo = {};
for n=1:length(model3d(no_model3d).inv_res)
	names_tomo{n} = [char( model3d(no_model3d).inv_res(n).name ), ', ',char( model3d(no_model3d).inv_res(n).tomo.date)] ;
end

set(handles.popupmenu_tomo1,'String',names_tomo,'Value',1);
set(handles.popupmenu_tomo2,'String',names_tomo,'Value',1);
popupmenu_Callback(hObject, eventdata, handles);
draw_cb(handles);

function popupmenu_model3d_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


function popupmenu_tomo1_Callback(hObject, eventdata, handles)
popupmenu_Callback(hObject, eventdata, handles);
draw_cb(handles);


function popupmenu_tomo1_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function popupmenu_tomo2_Callback(hObject, eventdata, handles)
popupmenu_Callback(hObject, eventdata, handles);
draw_cb(handles);

function popupmenu_tomo2_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function checkbox_colorbar_Callback(hObject, eventdata, handles)
draw_cb(handles);


function edit_min_Callback(hObject, eventdata, handles)
% val = str2double(get(hObject,'String'));
% h = getappdata(handles.fig_bh_timelapse, 'h');
% 
% h.tt_b_coul_min = val;
% if get(handles.checkbox_colorbar,'Value')==1
%     caxis(handles.axes_s1,[h.tt_b_coul_min h.tt_b_coul_max])
%     caxis(handles.axes_s2,[h.tt_b_coul_min h.tt_b_coul_max])
%     caxis(handles.axes_ds,[h.tt_b_coul_min h.tt_b_coul_max])
% end
% 
% setappdata(handles.fig_bh_timelapse, 'h',h);
draw_cb(handles)

function edit_min_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit_max_Callback(hObject, eventdata, handles)
% val = str2double(get(hObject,'String'));
% h = getappdata(handles.fig_bh_timelapse, 'h');
% 
% h.tt_b_coul_max = val;
% if get(handles.checkbox_colorbar,'Value')==1
% 	caxis(handles.axes_s1,[h.tt_b_coul_min h.tt_b_coul_max]);
%     caxis(handles.axes_s2,[h.tt_b_coul_min h.tt_b_coul_max]);
%     caxis(handles.axes_ds,[h.tt_b_coul_min h.tt_b_coul_max]);
% end
% 
% setappdata(handles.fig_bh_timelapse, 'h',h);
draw_cb(handles)

function edit_max_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function popupmenu_color_Callback(hObject, eventdata, handles)
maps = get(hObject,'String');
eval(['colormap(',maps{get(hObject,'Value')},')'])

function popupmenu_color_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function menu_file_Callback(hObject, eventdata, handles)


function menu_file_choose_db_Callback(hObject, eventdata, handles)
%[file, rep] = uigetfile('*.mat','Open DB');
%load([rep,file],'model3d');
h = getappdata(handles.fig_bh_timelapse, 'h');
[h.no_model3d,h.db_file,hh] = choisirPanneau('UserData', h.db_file,'3d');
if ishandle(hh), delete(hh), end
if isequal(h.no_model3d,0)
	return
end
load(h.db_file,'model3d')
set(handles.popupmenu_model3d,'Value',1);
set(handles.popupmenu_tomo1,'Value',1);
set(handles.popupmenu_tomo2,'Value',1);
if isempty(model3d)
    msgbox('Please define 3d panels and carry out inversion','warning','modal');
    set(handles.popupmenu_model3d,'String','model3d','Value',1);
    set(handles.popupmenu_tomo1,'String','Tomo1','Value',1);
    set(handles.popupmenu_tomo2,'String','Tomo2','Value',1);
    return;
else
    h.model3d = model3d;
    setappdata(handles.fig_bh_timelapse,'h',h);
    names_model3d = {};
    for n = 1:length(model3d)
        names_model3d{n} = char( model3d(n).name );
    end
    set(handles.popupmenu_model3d,'String',names_model3d);

    popupmenu_model3d_Callback(hObject, eventdata, handles);
    draw_cb(handles);
end
% --------------------------------------------------------------------
function menu_file_save_Callback(hObject, eventdata, handles)


function set_string_locale(handles)
h = getappdata(handles.fig_bh_timelapse, 'h');
set(handles.edit_min,'String',h.tt_b_coul_min);
set(handles.edit_max,'String',h.tt_b_coul_max);
set(handles.edit_min_diff,'String',h.tt_b_coul_min_diff);
set(handles.edit_max_diff,'String',h.tt_b_coul_max_diff);
set(handles.popupmenu_color,'Value',8);


function popupmenu_Callback(hObject, eventdata, handles)
h = getappdata(handles.fig_bh_timelapse,'h');
str = getappdata(handles.fig_bh_timelapse, 'str');
no_tomo = [get(handles.popupmenu_tomo1,'Value') get(handles.popupmenu_tomo2,'Value')];
panneau = h.model3d(get(handles.popupmenu_model3d,'Value'));
perc = get(handles.checkbox_absolut,'Value');
% handles.axes_s1 = subplot(4,3,[4 7 10]);
% display_s(handles,no_tomo(1),panneau,str,handles.axes_s1);
% handles.axes_s2 = subplot(4,3,[5 8 11]);
% display_s(handles,no_tomo(2),panneau,str,handles.axes_s2);
% handles.axes_ds = subplot(4,3,[6 9 12]);
% display_ds(handles,no_tomo,panneau,str,handles.axes_ds,perc);
% guidata(hObject,handles);
display_s(handles,no_tomo(1),panneau,str,h.axes_s1);
display_s(handles,no_tomo(2),panneau,str,h.axes_s2);
display_ds(handles,no_tomo,panneau,str,h.axes_ds,perc);
draw_cb(handles);

function display_s(handles,no_tomo,panneau,str,hAxes)
tomo = panneau.inv_res(no_tomo).tomo;
param = panneau.inv_res(no_tomo).param;
ssplot = reshape((tomo.s),length(tomo.x),length(tomo.y),length(tomo.z));
[nx,ny,nz] = size(ssplot);
ix = round(nx/2);
iy = round(ny/2);
iz = round(nz/2);
axes(hAxes);
if param.tomo_amp==0
	slice(1./ssplot,iy,ix,iz,'Parent',hAxes);
	titre = [panneau.inv_res(no_tomo).name,', ', str.s121,' [m/ns]'];
else
	slice(ssplot,iy,ix,iz,'Parent',hAxes);
	titre = [panneau.inv_res(no_tomo).name,', ', str.s177,' [Np/m]'];
end

set(gca,'DataAspectRatio',[1 1 1],'YDir','normal')
axis tight
title(titre,'FontSize',10)
xlabel(str.s119,'FontSize',8)
ylabel(str.s119,'FontSize',8)
zlabel(str.s120,'FontSize',8)

function display_ds(handles,no_tomo,panneau,str,hAxes,perc)
tomo1 = panneau.inv_res(no_tomo(1)).tomo;
tomo2 = panneau.inv_res(no_tomo(2)).tomo;
param = panneau.inv_res(no_tomo).param;
splot = reshape((tomo1.s),length(tomo1.x),length(tomo1.y),length(tomo1.z));
[nx,ny,nz] = size(splot);
ix = round(nx/2);
iy = round(ny/2);
iz = round(nz/2);
axes(hAxes);
if param.tomo_amp==0
    dsplot = reshape((1./tomo1.s-1./tomo2.s),length(tomo1.x),length(tomo1.y),length(tomo1.z)); 
%     if sum(sum(sum(tomo1.s))) ~= 0 || sum(sum(sum(tomo2.s))) ~= 0 || isnan(dsplot)
%         dsplot = zeros(length(tomo1.x),length(tomo1.y),length(tomo1.z));
%     end
	slice(dsplot,iy,ix,iz,'Parent',hAxes);
	titre = ['difference, ', str.s121,' [m/ns]'];
else
    dsplot = reshape((tomo1.s-tomo2.s),length(tomo1.x),length(tomo1.y),length(tomo1.z));
	slice(dsplot,iy,ix,iz,'Parent',hAxes);
	titre = ['difference, ', str.s177,' [Np/m]'];
end
if perc
    if sum(sum(sum(dsplot))) ~= 0
    dsplot = 100*dsplot.\reshape(1./tomo1.s,length(tomo1.x),length(tomo1.y), length(tomo1.z));         
    % else
       % msgbox('ds = 0, can''t display in percentage.','Message','modal');
        %return;
    end
    slice(dsplot,iy,ix,iz,'Parent',hAxes);
    titre = 'difference (%)';
end

set(gca,'DataAspectRatio',[1 1 1],'YDir','normal')
axis tight
title(titre,'FontSize',10)
xlabel(str.s119,'FontSize',8)
ylabel(str.s119,'FontSize',8)
zlabel(str.s120,'FontSize',8)

function edit_max_diff_Callback(hObject, eventdata, handles)
draw_cb(handles)

function edit_max_diff_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit_min_diff_Callback(hObject, eventdata, handles)
draw_cb(handles)

function edit_min_diff_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function draw_cb(handles)
h = getappdata(handles.fig_bh_timelapse,'h');
cbar = get(handles.checkbox_colorbar,'Value');
max1 = str2double(get(handles.edit_max,'String'));
max2 = str2double(get(handles.edit_max_diff,'String'));
min1 = str2double(get(handles.edit_min,'String'));
min2 = str2double(get(handles.edit_min_diff,'String'));

if cbar
    caxis(h.axes_s1,[min1 max1])
    caxis(h.axes_s2,[min1 max1])
    caxis(h.axes_ds,[min2 max2])
    axes(h.axes_s1); 
    colorbar
    axes(h.axes_s2); 
    colorbar
    axes(h.axes_ds); 
    colorbar
else
    axes(h.axes_s1); 
    colorbar off
    axes(h.axes_s2); 
    colorbar off
    axes(h.axes_ds); 
    colorbar off
end


% --- Executes on button press in checkbox_absolut.
function checkbox_absolut_Callback(hObject, eventdata, handles)
h = getappdata(handles.fig_bh_timelapse,'h');
str = getappdata(handles.fig_bh_timelapse, 'str');
no_tomo = [get(handles.popupmenu_tomo1,'Value') get(handles.popupmenu_tomo2,'Value')];
panneau = h.model3d(get(handles.popupmenu_model3d,'Value'));
hAxes = h.axes_ds;
perc = get(handles.checkbox_absolut,'Value');
display_ds(handles,no_tomo,panneau,str,hAxes,perc);
draw_cb(handles);

% --- Executes on button press in pushbutton_movie.
function pushbutton_movie_Callback(hObject, eventdata, handles)
h = getappdata(handles.fig_bh_timelapse,'h');
panneau = h.model3d(get(handles.popupmenu_model3d,'Value')); 
nbTomo = length(panneau.inv_res);

% str = getappdata(handles.fig_bh_timelapse, 'str');
% perc = get(handles.checkbox_absolut,'Value');

%hFig = figure('Visible', 'on', 'Units','characters');

% for frame=1:nbTomo
%      display_s(handles,frame,panneau,str,handles.axes_s1);
%      F(1,frame) = getframe;
%      display_s(handles,frame,panneau,str,handles.axes_s2);
%      F(2,frame) = getframe;
%      display_ds(handles,[1 frame],panneau,str,handles.axes_ds,perc);
%     F(3,frame) = getframe;
% end

% movie(handles.axes_s1,F(1,:),1,4);
% movie(handles.axes_s2,F(2,:),1,4);
% movie(handles.axes_ds,F(3,:),1,4);

%if isfield(h,'movieName') 
%movieName = h.movieName;
[FileName,PathName,FilterIndex] = uigetfile('*.jpeg','Display movie');
if FilterIndex
    [movieName rest]= strtok(FileName,'_');
    for frame = 1:nbTomo
        fName = [PathName movieName '_tomo_' num2str(frame) '.jpeg'];
        F.cdata = imread(fName);
        F.colormap = [];
        FT(1,frame) = F;
        fName = [PathName movieName '_diff_' num2str(frame) '.jpeg'];
        F.cdata = imread(fName);
        F.colormap = [];
        FD(1,frame) = F;
    end
    % axes_s1 = subplot(4,3,[4 7 10]);
    % axes_ds = subplot(4,3,[6 9 12]);
    movie(h.axes_s1,FT(1,:),1,4);
    movie(h.axes_ds,FD(1,:),1,4);
% else
%    msgbox('Please first save graphics.','Message','modal');
%    return;
end


% --- Executes on button press in pushbutton_save_fig.
function pushbutton_save_fig_Callback(hObject, eventdata, handles)
h = getappdata(handles.fig_bh_timelapse,'h');
if ~isfield(h,'model3d')
    msgbox('Please import data.','Message','modal');
    return;
elseif isempty(h.model3d)
    msgbox('Please add 3d panels in database module.','Message','modal');
    return;
end
nameDefault = 'Movies3d';
prompt = {'Enter folder name:                                                               '};
title = 'Save graphics';
nblines = 1;
answer = inputdlg(prompt,title,nblines,{nameDefault});

if ~isempty(answer)
    if length(answer{1})==0
        folderName='Movies3d';
    else
        folderName=answer{1};    
    end
else
    return;
end
DefaultFileName = 'movie3d';
[fileName, pathName, filterindex] = uiputfile('.jpeg','Save graphics',DefaultFileName);
if filterindex   
    A = exist([pathName folderName],'dir');
    if A == 0
        [status,message,messageid] = mkdir(pathName,folderName);
    else
%         pathName = [pathName folderName '\'];
        status = 1;
    end
    movieName = fileName(1:end-5);
    if status
        str = getappdata(handles.fig_bh_timelapse, 'str');
        model3d = h.model3d(get(handles.popupmenu_model3d,'Value')); 
        perc = get(handles.checkbox_absolut,'Value');
        nbTomo = length(model3d.inv_res);
        hFig = figure('Visible', 'off', 'Units','characters');
        % axes_s1 = subplot(1,2,1);
        % axes_ds = subplot(1,2,2);
        for frame=1:nbTomo
             display_s(handles,frame,model3d,str,h.axes_s1);
             display_s(handles,frame,model3d,str,h.axes_s2);
             F = getframe;
             fName = [pathName folderName '\' movieName '_tomo_' num2str(frame) '.jpeg'];
             imwrite(F.cdata,fName);
             %imwrite(F,'movie_tomo','fmt','jpeg','Compression','none','WriteMode','append');
             display_ds(handles,[1 frame],model3d,str,h.axes_ds,perc);
             F = getframe;
             fName = [pathName folderName '\' movieName '_diff_' num2str(frame) '.jpeg'];
             imwrite(F.cdata,fName);
        end
        delete(hFig);
        h.movieName = movieName;
        setappdata(handles.fig_bh_timelapse,'h',h);
    end
end
