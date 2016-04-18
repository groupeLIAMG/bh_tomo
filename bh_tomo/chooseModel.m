function [modelNo,file,rep] = chooseModel(varargin)

modelNo=[];
rep='';
file='';
rep2='';
file2=0;
if nargin>=2
    rep=varargin{1};
    file=varargin{2};
end

fs = 11;
if nargin>=3
    fs = varargin{3};
elseif ispc
    fs = 9;
end
vScale = 1;
if ispc
    vScale = 0.81;
end

width = 300;
height = 175*vScale;
% get screen size
su = get(groot,'Units');
set(groot,'Units','points')
scnsize = get(groot,'ScreenSize');
pos = [scnsize(3)/2-width/2 scnsize(4)/2-height/2 width height];
set(groot,'Units',su)       % Restore default root screen units

f = figure('Visible','off',...
    'Units','points',...
    'Position',pos,...
    'Tag','fig_chooseModel',...
    'Name','Choose Model',...
    'NumberTitle','off',...
    'ToolBar','none',...,
    'MenuBar','None',...
    'Resize','off',...
    'CloseRequestFcn',@closeWindow);

hdbName = uicontrol('Style','text',...
    'BackgroundColor','white',...
    'FontSize',fs+1,...
    'HorizontalAlignment','center',...
    'Units','normalized',...
    'Position',[0.1 0.77 0.8 0.142],...
    'Parent',f);
if ~isempty(file)
    hdbName.String = file;
end

uicontrol('Style','pushbutton',...
    'String','Choose DB',...
    'FontSize',fs,...
    'Units','normalized',...
    'Position',[0.2 0.54 0.6 0.13],...
    'Callback',@chooseDB,...
    'Parent',f);

hlist = uicontrol('Style','popupmenu',...
    'String',' ',...
    'FontSize',fs,...
    'Units','normalized',...
    'Position',[0.2 0.31 0.6 0.142],...
    'Parent',f);
if ~isempty(file)
    load([rep,file],'models')
    names = cell(length(models));
    for n=1:length(models)
        names{n} = models(n).name;
    end
    hlist.String = names;
end


uicontrol('Style','pushbutton',...
    'String','Cancel',...
    'FontSize',fs,...
    'Units','normalized',...
    'Position',[0.1 0.085 0.35 0.13],...
    'Callback',@cancel,...
    'Parent',f);

uicontrol('Style','pushbutton',...
    'String','Done',...
    'FontSize',fs,...
    'Units','normalized',...
    'Position',[0.55 0.085 0.35 0.13],...
    'Callback',@validate,...
    'Parent',f);

uiwait(f)

    function closeWindow(varargin)
        delete(f)
    end
    function chooseDB(varargin)
        [file2, rep2] = uigetfile('*.mat','Open DB');
        if isequal(file2,0)
            return
        end
        load([rep2,file2],'models')
        names = cell(length(models));
        for nn=1:length(models)
            names{nn} = models(nn).name;
        end
        hlist.String = names;
    end
    function cancel(varargin)
        modelNo = [];
        closeWindow()
    end
    function validate(varargin)
        if isempty(file) && isequal(file2,0)
            cancel()
            return
        end
        if ~isequal(file2,0)
            file=file2;
            rep=rep2;
        end
        modelNo = str2double(hlist.Value);
        closeWindow()
    end
end

