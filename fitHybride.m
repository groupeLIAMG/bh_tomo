function varargout = fitHybride(varargin)
% FITHYBRIDE M-file for fitHybride.fig
%      FITHYBRIDE, by itself, creates a new FITHYBRIDE or raises the existing
%      singleton*.
%
%      H = FITHYBRIDE returns the handle to a new FITHYBRIDE or the handle to
%      the existing singleton*.
%
%      FITHYBRIDE('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in FITHYBRIDE.M with the given input arguments.
%
%      FITHYBRIDE('Property','Value',...) creates a new FITHYBRIDE or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before fitHybride_OpeningFunction gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to fitHybride_OpeningFcn via varargin.
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

% Edit the above text to modify the response to help fitHybride

% Last Modified by GUIDE v2.5 17-Dec-2012 14:37:16

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @fitHybride_OpeningFcn, ...
                   'gui_OutputFcn',  @fitHybride_OutputFcn, ...
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



% --- Executes just before fitHybride is made visible.
function fitHybride_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to fitHybride (see VARARGIN)

% Choose default command line output for fitHybride
handles.output = 0;
handles.second_output = hObject;
handles.third_output = [];

% Update handles structure
guidata(hObject, handles);

h = varargin{1};
h.uTx = unique(h.Tx,'rows');
nTx = size(h.uTx,1);
mA = [h.lrai ones(size(h.lrai))]\h.Acorr';
h.A0 = mA(2)*ones(1,nTx);
h.data = zeros(size(h.Acorr));
h.sigma_s2 = zeros(1,nTx);
h.tau = zeros(length(h.data),2);
for n=1:nTx
	ind = ( h.Tx(:,3) == h.uTx(n,3) & ...
		h.Tx(:,2) == h.uTx(n,2) & h.Tx(:,1) == h.uTx(n,1));
	h.sigma_s2(n) = median(h.scentroid(ind));
	h.data(ind) = h.f * h.fcentroid(ind) / h.sigma_s2(n);
end
h.tau(:,2) = (h.A0(1) - h.Acorr)';

mF = [h.lrai ones(size(h.lrai))]\h.data';
h.f0 = mF(2)*ones(1,nTx);
h.f = h.f*ones(1,nTx);
h.att_app = 0.5*mA(1)+0.5*mF(1);
for n=1:nTx
	ind = ( h.Tx(:,3) == h.uTx(n,3) & ...
		h.Tx(:,2) == h.uTx(n,2) & h.Tx(:,1) == h.uTx(n,1));
	h.tau(ind,1) = h.f0(1) - h.data(ind);
end
h.poids = [0.5;   % fc
	       0.5];  % App
h.tauDC = [];

h.noTx = 1;

h.str = get_str_locale();
set_String_locale(handles, h.str)
if nargin > 4
	params = varargin{2};
	handles.third_output = params;
	if ~isempty(params)
		if isfield(params,'att_app'), h.att_app = params.att_app; end
		if isfield(params,'A0')
			if length(params.A0)==1
				h.A0(:) = params.A0;
			else
				h.A0 = params.A0;
			end
		end
		if isfield(params,'f0')
			if length(params.f0)==1
				h.f0(:) = params.f0;
			else
				h.f0 = params.f0;
			end
		end
		if isfield(params,'f')
			if length(params.f)==1
				h.f(:) = params.f;
			else
				h.f = params.f;
			end
		end
		if isfield(params,'poids_fc'), h.poids(1) = params.poids_fc; end
		if isfield(params,'poids_App'), h.poids(2) = params.poids_App; end
	end
	setappdata(handles.fig_hybrid, 'h', h)
	calculeTau(handles);
	update_plots(handles)
else
	setappdata(handles.fig_hybrid, 'h', h)
	update_plots(handles)
end

set(handles.listbox_no_Tx,'String',num2str((1:nTx)'))
set(handles.edit_poids_App,'string',num2str(h.poids(2)))
set(handles.edit_poids_fc,'string',num2str(h.poids(1)))
set(handles.text_xTx,'String',num2str(h.uTx(1,1)))
set(handles.text_yTx,'String',num2str(h.uTx(1,2)))
set(handles.text_zTx,'String',num2str(h.uTx(1,3)))

update_infos(handles)

% UIWAIT makes fitHybride wait for user response (see UIRESUME)
uiwait(handles.fig_hybrid);



% --- Outputs from this function are returned to the command line.
function varargout = fitHybride_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;
varargout{2} = handles.second_output;
varargout{3} = handles.third_output;

function update_infos(handles)
h = getappdata(handles.fig_hybrid, 'h');
set(handles.edit_A0,'string',num2str(h.A0(h.noTx)))
set(handles.edit_f0,'string',num2str(h.f0(h.noTx)*...
	h.sigma_s2(h.noTx)/(h.ffac*h.f(h.noTx))))
set(handles.edit_f,'String',num2str(h.f(h.noTx)/h.ffac))
set(handles.edit_pente,'string',num2str(-h.att_app))


function update_plots(handles)
h = getappdata(handles.fig_hybrid, 'h');
sigma_s2 = h.sigma_s2;
if get(handles.checkbox_Tx,'Value')==1
	update_plots_ind(handles)
else
	update_plots_coll(handles)
	sigma_s2(:) = median(h.scentroid);
end
nTx = 1:size(h.uTx,1);
axes(handles.axes7)
plot(nTx, h.A0)
hold on
plot(nTx(h.noTx), h.A0(h.noTx),'ro')
hold off
ylabel('A_0')

axes(handles.axes8)
plot(nTx, h.f0.*sigma_s2./(h.ffac*h.f))
hold on
plot(nTx(h.noTx), h.f0(h.noTx).*sigma_s2(h.noTx)./(h.ffac*h.f(h.noTx)),'ro')
hold off
ylabel('f_S')
axes(handles.axes9)
plot(nTx, h.f/h.ffac)
hold on
plot(nTx(h.noTx), h.f(h.noTx)/h.ffac,'ro')
hold off
xlabel('Tx no')
ylabel('f')


function update_plots_coll(handles)
h = getappdata(handles.fig_hybrid, 'h');

axes(handles.axes1)
plot(h.lrai, h.Acorr, 'gd')
xlim=get(gca,'XLim');
hold on
plot(xlim,[xlim' [1 1]']*[h.att_app; h.A0(h.noTx)],'k-')
hold off
xlabel(h.str.s86)
ylabel(h.str.s85)
title(['A_0 = ',num2str(h.A0(h.noTx))])

axes(handles.axes2)
plot(h.lrai, h.data,'co')
xlim=get(gca,'XLim');
hold on
plot(xlim,[xlim' [1 1]']*[h.att_app; h.f0(h.noTx)],'k-')
hold off
xlabel(h.str.s86)
ylabel([str.s82,' [',h.fUnits,'^{-1}]'])
title(['f_S = ',num2str(h.f0(h.noTx)*h.sigma_s2(h.noTx)/(h.ffac*h.f(h.noTx)))])

axes(handles.axes3)
plot(h.lrai, h.tau(:,1),'co')
hold on
plot(h.lrai, h.tau(:,2),'gd')
hold off
xlabel(h.str.s86)
ylabel('\tau')
title(h.str.s200)

axes(handles.axes4)
if ~isempty(h.tauDC)
	plot(h.lrai, h.tauDC./h.lrai,'m*')
	hold on
	plot(h.lrai, (h.tau*h.poids)./h.lrai,'cs')
	hold off
else
	plot(h.lrai, (h.tau*h.poids)./h.lrai,'cs')
end
ylabel('\alpha_a')
xlabel(h.str.s86)
title(h.str.s201)

axes(handles.axes10)
if ~isempty(h.tauDC)
	plot(h.lrai, h.tauDC,'r*')
	hold on
	plot(h.lrai, (h.tau*h.poids),'k+')
	hold off
else
	plot(h.lrai, (h.tau*h.poids),'k+')
end
ylabel('\tau')
xlabel(h.str.s86)


function update_plots_ind(handles)
h = getappdata(handles.fig_hybrid, 'h');
[ind,indTx] = getInd(handles);

axes(handles.axes1)
plot(h.lrai(~ind), h.Acorr(~ind), 'gd')
hold on
plot(h.lrai(ind), h.Acorr(ind), 'kd')
hold off
xlim=get(gca,'XLim');
hold on
plot(xlim,[xlim' [1 1]']*[h.att_app; h.A0(h.noTx)],'k-')
hold off
xlabel(h.str.s86)
ylabel(h.str.s85)
title(['A_0 = ',num2str(h.A0(h.noTx))])

axes(handles.axes2)
plot(h.lrai(~ind), h.data(~ind),'co')
hold on
plot(h.lrai(ind), h.data(ind),'ko')
hold off
xlim=get(gca,'XLim');
hold on
plot(xlim,[xlim' [1 1]']*[h.att_app; h.f0(h.noTx)],'k-')
hold off
xlabel(h.str.s86)
ylabel([str.s82,' [',h.fUnits,'^{-1}]'])
title(['f_0 = ',num2str(h.f0(h.noTx)*h.sigma_s2(h.noTx)/(h.ffac*h.f(h.noTx)))])

axes(handles.axes3)
plot(h.lrai(~ind), h.tau(~ind,1),'co')
hold on
plot(h.lrai(~ind), h.tau(~ind,2),'gd')
plot(h.lrai(ind), h.tau(ind,1),'ko')
plot(h.lrai(ind), h.tau(ind,2),'kd')
hold off
xlabel(h.str.s86)
ylabel('\tau')
title(h.str.s200)

axes(handles.axes4)
if ~isempty(h.tauDC)
	plot(h.lrai(~ind), h.tauDC(~ind)./h.lrai(~ind),'m*')
	hold on
	plot(h.lrai(~ind), (h.tau(~ind,:)*h.poids)./h.lrai(~ind),'cs')
	plot(h.lrai(ind), h.tauDC(ind)./h.lrai(ind),'k*')
	plot(h.lrai(ind), (h.tau(ind,:)*h.poids)./h.lrai(ind),'ks')
	hold off
else
	plot(h.lrai(~ind), (h.tau(~ind,:)*h.poids)./h.lrai(~ind),'cs')
	hold on
	plot(h.lrai(ind), (h.tau(ind,:)*h.poids)./h.lrai(ind),'ks')
	hold off
end
ylabel('\alpha_a')
xlabel(h.str.s86)
title(h.str.s201)

axes(handles.axes10)
if ~isempty(h.tauDC)
	plot(h.lrai(~ind), h.tauDC(~ind),'r*')
	hold on
	plot(h.lrai(~ind), (h.tau(~ind,:)*h.poids),'k+')
	plot(h.lrai(ind), h.tauDC(ind),'y*')
	plot(h.lrai(ind), (h.tau(ind,:)*h.poids),'y+')
	hold off
else
	plot(h.lrai(~ind), (h.tau(~ind,:)*h.poids),'k+')
	hold on
	plot(h.lrai(ind), (h.tau(ind,:)*h.poids),'y+')
	hold off
end
ylabel('\tau')
xlabel(h.str.s86)


function pushbutton_OK_Callback(hObject, eventdata, handles)
h = getappdata(handles.fig_hybrid, 'h');
handles.output = h.tau*h.poids;
params.att_app = h.att_app;
params.A0 = h.A0;
params.f0 = h.f0;
params.f = h.f;
params.poids_fc = h.poids(1);
params.poids_App = h.poids(2);
handles.third_output = params;

% Update handles structure
guidata(hObject, handles);
uiresume;


function pushbutton_cancel_Callback(hObject, eventdata, handles)
handles.output = 0;

% Update handles structure
guidata(hObject, handles);
uiresume;


function set_String_locale(handles, str)
%set(handles.uipanel_A,              'Title',  str.s89);
%set(handles.text_penteApp,         'String', str.s90);
set(handles.pushbutton_cancel,     'String', str.s91);



function calculeTau(handles,varargin)
h = getappdata(handles.fig_hybrid, 'h');
if nargin >=2
	no = varargin{1};
else
	no = h.noTx;
end
[ind,indTx] = getInd(handles,no);
if get(handles.checkbox_Tx,'Value')==1
	h.data(indTx) = h.f(no) * h.fcentroid(indTx) / h.sigma_s2(no);
else
%	for n=1:size(h.uTx,1)
%		ind = ( h.Tx(:,3) == h.uTx(n,3) & ...
%			h.Tx(:,2) == h.uTx(n,2) & h.Tx(:,1) == h.uTx(n,1));
%		h.data(ind) = h.f(n) * h.fcentroid(ind) / h.sigma_s2(n);
%	end
	sigma_s2 = median(h.scentroid);
	h.data = h.f(1) * h.fcentroid / sigma_s2;
end
h.tau(indTx,1) = (h.f0(no) - h.data(indTx))';
h.tau(indTx,2) = (h.A0(no) - h.Acorr(indTx))';
setappdata(handles.fig_hybrid, 'h', h)


function edit_A0_Callback(hObject, eventdata, handles)
h = getappdata(handles.fig_hybrid, 'h');
if get(handles.checkbox_Tx,'Value')==1
	h.A0(h.noTx) = str2num(get(hObject,'string'));
else
	h.A0(:) = str2num(get(hObject,'string'));
end
setappdata(handles.fig_hybrid, 'h', h)
calculeTau(handles)
update_plots(handles)


function edit_A0_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


function edit_pente_Callback(hObject, eventdata, handles)
h = getappdata(handles.fig_hybrid, 'h');
h.att_app = -str2num(get(hObject,'string'));
setappdata(handles.fig_hybrid, 'h', h)
calculeTau(handles)
update_plots(handles)


function edit_pente_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


function edit_poids_App_Callback(hObject, eventdata, handles)
h = getappdata(handles.fig_hybrid, 'h');
h.poids(2) = str2num(get(hObject,'string'));
h.poids(1) = 1-h.poids(2);
set(handles.edit_poids_fc,'string',num2str(h.poids(1)))
setappdata(handles.fig_hybrid, 'h', h)
calculeTau(handles)
update_plots(handles)

function edit_poids_App_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


function edit_f0_Callback(hObject, eventdata, handles)
h = getappdata(handles.fig_hybrid, 'h');
h.f0(h.noTx) = h.ffac*h.f(h.noTx)*str2num(get(hObject,'string'))/h.sigma_s2(h.noTx);
setappdata(handles.fig_hybrid, 'h', h)
calculeTau(handles)
update_plots(handles)


function edit_f0_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


function edit_poids_fc_Callback(hObject, eventdata, handles)
h = getappdata(handles.fig_hybrid, 'h');
h.poids(1) = str2num(get(hObject,'string'));
h.poids(2) = 1-h.poids(1);
set(handles.edit_poids_App,'string',num2str(h.poids(2)))
setappdata(handles.fig_hybrid, 'h', h)
calculeTau(handles)
update_plots(handles)


function edit_poids_fc_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


function edit_f_Callback(hObject, eventdata, handles)
h = getappdata(handles.fig_hybrid, 'h');
old_f = h.f;
%if get(handles.checkbox_Tx,'Value')==1
%	h.f(h.noTx) = h.ffac*str2num(get(hObject,'string'));
%else
	h.f(:) = h.ffac*str2num(get(hObject,'string'));
%end
h.f0 = h.f0.*h.f./old_f;
set(handles.edit_f0,'string',num2str(h.f0(h.noTx)*h.sigma_s2(h.noTx)/(h.ffac*h.f(h.noTx))))
setappdata(handles.fig_hybrid, 'h', h)
calculeTau(handles)
update_plots(handles)


function edit_f_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


function pushbutton_alphaDC_Callback(hObject, eventdata, handles)
h = getappdata(handles.fig_hybrid, 'h');
load(h.db_file,'panels');
if isempty(panels), return, end
tomo = [];
for n=1:length(panels)
	if strcmp(panels(n).name, h.nom_Ldc{1})  %nom-->name
		break;
	end
end
for nn=1:length(panels(n).inv_res)
	if strcmp(panels(n).inv_res(nn).nom,h.nom_Ldc{2})
		tomo = panels(n).inv_res(nn).tomo;
		break;
	end
end
if isempty(tomo), return, end
[filename, rep] = uigetfile({'*.xyz';'*.*'},'Modele de conductivite DC');
if isequal(filename,0) || isequal(rep,0)
	return
end
sigmaDC = load([rep,filename]);
if ~exist('sigmaDC')
	return
end
[xi,zi] = meshgrid(tomo.x,tomo.z);

sigmaDC = griddata(sigmaDC(:,1), sigmaDC(:,2), sigmaDC(:,3), xi, zi);
sigmaDC = reshape(sigmaDC,numel(sigmaDC),1);

mu0 = 4e-7*pi*1e9;  % 1e9 -> m/ns  m/s
alphaDC = 0.5*mu0*sigmaDC./tomo.s;

figure
subplot(121)
imagesc(tomo.x,tomo.z,reshape(sigmaDC,length(tomo.z), length(tomo.x)))
set(gca,'DataAspectRatio',[1 1 1],'YDir','normal')
colorbar
xlabel(h.str.s119)
ylabel(h.str.s120)
title('\sigma_{DC}','FontSize',14)

subplot(122)
imagesc(tomo.x,tomo.z,reshape(alphaDC,length(tomo.z), length(tomo.x)))
set(gca,'DataAspectRatio',[1 1 1],'YDir','normal')
colorbar
xlabel(h.str.s119)
ylabel(h.str.s120)
title('\alpha_{DC}','FontSize',14)

ind = [];
for n=1:length(h.no_traces)
	ii = find( tomo.no_trace==h.no_traces(n) );
	if isempty(ii)
		uiwait(warndlg([h.str.s186{1}, num2str(h.no_traces(n)), h.str.s186{2}]))
		return
	else
		ind = [ind ii];
	end
end
alphaDC(isnan(alphaDC)) = 0;

h.tauDC = tomo.L(ind,:)*alphaDC;
setappdata(handles.fig_hybrid, 'h', h)
update_plots(handles)


function listbox_no_Tx_Callback(hObject, eventdata, handles)
h = getappdata(handles.fig_hybrid, 'h');
h.noTx = get(hObject,'Value');
set(handles.text_xTx,'String',num2str(h.uTx(h.noTx,1)))
set(handles.text_yTx,'String',num2str(h.uTx(h.noTx,2)))
set(handles.text_zTx,'String',num2str(h.uTx(h.noTx,3)))
setappdata(handles.fig_hybrid, 'h', h)
calculeTau(handles)
update_plots(handles)
update_infos(handles)

function listbox_no_Tx_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


function checkbox_Tx_Callback(hObject, eventdata, handles)
h = getappdata(handles.fig_hybrid, 'h');
if get(hObject,'Value')==1
	set(handles.pushbutton_fitA0_1,'Enable','on')
	set(handles.pushbutton_fitF0_1,'Enable','on')
	for n=1:size(h.uTx,1)
		calculeTau(handles,n)
	end
else
	set(handles.pushbutton_fitA0_1,'Enable','off')
	set(handles.pushbutton_fitF0_1,'Enable','off')
	calculeTau(handles)
end
update_plots(handles)

function popupmenu_bin_Callback(hObject, eventdata, handles)
calculeTau(handles)
update_plots(handles)


function popupmenu_bin_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function [ind,indTx] = getInd(handles,varargin)
h = getappdata(handles.fig_hybrid, 'h');
if get(handles.checkbox_Tx,'Value')==0
	ind = true(size(h.Tx(:,3)));
	indTx = true(size(h.Tx(:,3)));
	return
end
if nargin>=2
	nos = varargin{1};
else
	nos = h.noTx;
end
no = nos;
s=get(handles.popupmenu_bin,'String');
n=str2num(s{get(handles.popupmenu_bin,'Value')});
n = floor(n/2);
nos = (nos-n : n+nos);
nos = nos(nos>=1 & nos<=numel(h.uTx));

ind = false(size(h.Tx(:,3)));
for n=nos
	ind = ind | ( h.Tx(:,3) == h.uTx(n,3) & ...
		h.Tx(:,2) == h.uTx(n,2) & h.Tx(:,1) == h.uTx(n,1));
end
indTx = h.Tx(:,3) == h.uTx(no,3) & ...
		h.Tx(:,2) == h.uTx(no,2) & h.Tx(:,1) == h.uTx(no,1);

function pushbutton_fitF0_1_Callback(hObject, eventdata, handles)
h = getappdata(handles.fig_hybrid, 'h');
[ind,indTx] = getInd(handles);
h.f(h.noTx) = h.ffac*str2num(get(handles.edit_f,'String'));
data = h.f(h.noTx) * h.fcentroid(ind) / h.sigma_s2(h.noTx);
if get(handles.checkbox_fixe_att,'Value')==0
	mF = [h.lrai(ind) ones(size(h.lrai(ind)))]\data';
	h.f0(h.noTx) = mF(2);
	h.att_app = mF(1);
else
	h.f0(h.noTx) = ones(size(h.lrai(ind)))\(data'-h.lrai(ind)*h.att_app);
end
setappdata(handles.fig_hybrid, 'h', h)
calculeTau(handles)
update_plots(handles)
update_infos(handles)

function pushbutton_fitA0_1_Callback(hObject, eventdata, handles)
h = getappdata(handles.fig_hybrid, 'h');
[ind,indTx] = getInd(handles);
if get(handles.checkbox_fixe_att,'Value')==0
	mA = [h.lrai(ind) ones(size(h.lrai(ind)))]\h.Acorr(ind)';
	h.A0(h.noTx) = mA(2);
	h.att_app = mA(1);
else
	h.A0(h.noTx) = ones(size(h.lrai(ind)))\(h.Acorr(ind)'-h.lrai(ind)*h.att_app);
end
setappdata(handles.fig_hybrid, 'h', h)
calculeTau(handles)
update_plots(handles)
update_infos(handles)

function pushbutton_fitA0_a_Callback(hObject, eventdata, handles)
h = getappdata(handles.fig_hybrid, 'h');
if get(handles.checkbox_fixe_att,'Value')==0
	mA = [h.lrai ones(size(h.lrai))]\h.Acorr';
	h.A0(:) = mA(2);
	h.att_app = mA(1);
else
	h.A0(:) = ones(size(h.lrai))\(h.Acorr'-h.lrai*h.att_app);
end
setappdata(handles.fig_hybrid, 'h', h)
calculeTau(handles)
update_plots(handles)
update_infos(handles)


function pushbutton_fitF0_a_Callback(hObject, eventdata, handles)
h = getappdata(handles.fig_hybrid, 'h');
h.f(:) = h.ffac*str2num(get(handles.edit_f,'String'));
%for n=1:size(h.uTx,1)
%	ind = ( h.Tx(:,3) == h.uTx(n,3) & ...
%		h.Tx(:,2) == h.uTx(n,2) & h.Tx(:,1) == h.uTx(n,1));
%	h.data(ind) = h.f(n) * h.fcentroid(ind) / h.sigma_s2(n);
%end
sigma_s2 = median(h.scentroid);
h.data = h.f(1) * h.fcentroid / sigma_s2;
if get(handles.checkbox_fixe_att,'Value')==0
	mF = [h.lrai ones(size(h.lrai))]\h.data';
	h.f0(:) = mF(2);
	h.att_app = mF(1);
else
	h.f0(:) = ones(size(h.lrai))\(h.data'-h.lrai*h.att_app);
end
setappdata(handles.fig_hybrid, 'h', h)
calculeTau(handles)
update_plots(handles)
update_infos(handles)


function checkbox_fixe_att_Callback(hObject, eventdata, handles)
