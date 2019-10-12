function bh_tomo2_db(varargin)
%BH_TOMO2_DB

rep='';
file='';
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

auto_pick = [];
saved = true;

width = 1000;
height = 800*vScale;
% get screen size
su = get(groot,'Units');
set(groot,'Units','points')
scnsize = get(groot,'ScreenSize');
pos = [scnsize(3)/2-width/2 scnsize(4)/2-height/2 width height];
set(groot,'Units',su)       % Restore default root screen units

f = figure('Visible','off',...
    'Units','points',...
    'Position',pos,...
    'Tag','fig_bh_tomo2_db',...
    'Name','bh_tomo_db',...
    'NumberTitle','off',...
    'ToolBar','none',...,
    'MenuBar','None',...
    'SizeChangedFcn',@resizeUI,...
    'CloseRequestFcn',@closeWindow);


%
% Menu
%
hmenu = uimenu(f,'Label','File');
uimenu(hmenu,'Label','Open ...',...
    'Accelerator','O',...
    'Callback',@openFile);
uimenu(hmenu,'Label','Save',...
    'Accelerator','S',...
    'Callback',@saveFile);
uimenu(hmenu,'Label','Save As ...',...
    'Accelerator','A',...
    'Callback',@saveFileAs);
uimenu(hmenu,'Label','Close',...
    'Separator','on',...
    'Accelerator','W',...
    'Callback',@closeWindow);

%
% Create Borehole panel
%
hborehole = BoreholeUI(f,'Units','points');
hborehole.FontSize = fs;

%
% Create MOG panel
%
hmog = MogUI(hborehole,f,'Units','points');
hmog.FontSize = fs;

%
% Create Model panel
%
hmodel= ModelUI(hborehole,hmog,f,'Units','points');
hmodel.FontSize = fs;

%
% Create info panel
%
hinfo = uipanel(f,'Title','Infos','Units','points');
hinfo.FontSize = fs+1;
htextinfo = uicontrol('Style','text',...
    'Parent',hinfo,...
    'BackgroundColor',[1 1 1],...
    'Units','normalized',...
    'HorizontalAlignment','center',...
    'FontSize',fs+1,...
    'Position',[0.05 0.05 0.9 0.9]);
htextinfo.String = {'','','Database: ','',...
    '0 Borehole(s)','0 MOG(s)','0 Model(s)','','0 traces'};


%
% Add listeners
%
addlistener(hborehole,'boreholeAdded',@updateBHinfo);
addlistener(hborehole,'boreholeDeleted',@updateBHinfo);
addlistener(hborehole,'boreholeAdded',@hmog.updateBHlist);
addlistener(hborehole,'boreholeDeleted',@hmog.updateBHlist);
addlistener(hborehole,'boreholeEdited',@dbEdited);

addlistener(hmog,'mogAdded',@updateMogInfo);
addlistener(hmog,'mogDeleted',@mogDeleted);
addlistener(hmog,'mogEdited',@dbEdited);

addlistener(hmodel,'modelAdded',@updateModelInfo);
addlistener(hmodel,'modelDeleted',@updateModelInfo);
addlistener(hmodel,'modelEdited',@dbEdited);

if ~isempty(file)
    loadDB();
end
f.Visible = 'on';



    function resizeUI(varargin)
        f.Visible = 'off';
        hBorder = 15;
        width = f.Position(3);
        height = f.Position(4);
        vBorderTop = 20*vScale;
        
        vSize1 = 4*height/7 - vBorderTop;
        vSize2 = 3*height/7 - vBorderTop;
        hsize = width/3 - hBorder;
        
        hborehole.Position = [hBorder vSize2+2*hBorder hsize vSize1];
        hmog.Position = [2*hBorder+hsize vSize2+2*hBorder 2*hsize vSize1];
        
        hmodel.Position = [hBorder hBorder 2*hsize vSize2];
        hinfo.Position = [width-hsize-hBorder hBorder hsize vSize2];
        
        f.Visible = 'on';
    end

    function dbEdited(varargin)
        saved = false;
    end
    function updateBHinfo(varargin)
        htextinfo.String{5} = [num2str(numel(hborehole.boreholes)),' Borehole(s)'];
        saved = false;
    end
    function updateMogInfo(varargin)
        htextinfo.String{6} = [num2str(numel(hmog.mogs)),' MOG(s)'];
        if numel(hmog.mogs)>0
            nt = 0;
            for n=1:numel(hmog.mogs)
                nt = nt + hmog.mogs(n).data.ntrace;
            end
            htextinfo.String{9} = [num2str(nt),' traces'];
        else
            htextinfo.String{9} = '0 traces';
        end
        saved = false;
    end
    function updateModelInfo(varargin)
        htextinfo.String{7} = [num2str(numel(hmodel.models)),' Model(s)'];
        saved = false;
    end

    function openFile(varargin)
        [file, rep] = uigetfile('*.mat','Open Database');
        if isequal(file,0)
            return
        end
        loadDB();
    end
    function loadDB(varargin)
        tmp = load([rep,file]);
        if isfield(tmp,'panels')
            warndlg(['File ',file,' is version 1 database, please convert first'])
            return
        end
        
        if ~isfield(tmp,'boreholes')
            errordlg(['File ',file,' not a bh_tomo database'],'File error')
            return
        end
        hborehole.boreholes = tmp.boreholes;
        hmog.air = tmp.air;
        hmog.mogs = tmp.mogs;
        hmodel.models = tmp.models;
        
        auto_pick = tmp.auto_pick; 
        
        htextinfo.String{3} = ['Database: ',file];
        updateBHinfo()
        updateMogInfo()
        updateModelInfo()
        hmog.updateBHlist()
        saved = true;
    end
    function saveFile(varargin)
        if isempty(file)
            [file, rep] = uiputfile('*.mat','Save Database');
            if isequal(file,0)
                return
            end
            htextinfo.String{3} = ['Database: ',file];
        end
        hw = warndlg(['Saving ',file]);
        db_file = [rep,file];
        names_mog = cell(1,length(hmog.mogs));
        for n=1:length(hmog.mogs)
            names_mog{n} = hmog.mogs(n).name;
        end
        models = hmodel.models; %#ok<NASGU>
        boreholes = hborehole.boreholes; %#ok<NASGU>
        mogs = hmog.mogs; %#ok<NASGU>
        air = hmog.air; %#ok<NASGU>
        save(db_file,'names_mog','mogs','air','boreholes','models','auto_pick','-v7.3')
        saved = true;
        delete(hw)
    end
    function saveFileAs(varargin)
        [file, rep] = uiputfile('*.mat','Save Database');
        if isequal(file,0)
            return
        end
        hw = warndlg(['Saving ',file]);
        db_file = [rep,file];
        htextinfo.String{3} = ['Database: ',file];
        names_mog = cell(1,length(hmog.mogs));
        for n=1:length(hmog.mogs)
            names_mog{n} = hmog.mogs(n).name;
        end
        models = hmodel.models; %#ok<NASGU>
        boreholes = hborehole.boreholes; %#ok<NASGU>
        mogs = hmog.mogs; %#ok<NASGU>
        air = hmog.air; %#ok<NASGU>
        save(db_file,'names_mog','mogs','air','boreholes','models','auto_pick','-v7.3')
        saved = true;
        delete(hw)
    end
    function closeWindow(varargin)
        if saved == false
            choice = questdlg('Database not saved, quit anyway?',...
                'bh_tomo_db',...
                'Don''t save','Cancel','Save','Save');
            switch choice
                case 'Don''t save'
                case 'Cancel'
                    return
                case 'Save'
                    saveFile()
            end
        end
        quitUI()
    end
    function quitUI(varargin)
        delete(f)
    end
    function mogDeleted(varargin)
        eventData = varargin{2};
        deletedMog = eventData.number;
        for n=1:length(hmodel.models)
            ind = hmodel.models(n).mogs~=deletedMog;
            hmodel.models(n).mogs = hmodel.models(n).mogs(ind);
            for nn=1:length(hmodel.models(n).mogs)
                if hmodel.models(n).mogs(nn) > deletedMog
                    hmodel.models(n).mogs(nn) = hmodel.models(n).mogs(nn)-1;
                end
            end
            if ~all(ind)
                % model contained deleted mog, we should rebuild de grid
                warndlg(['Mog removed from model ',hmodel.models(n).name,...
                    ', reinitializing grid'])
                hmodel.models(n).grid = Grid.empty;
            end
        end
        
        if numel(auto_pick) >= deletedMog
            found = -1;
            for n=1:numel(auto_pick)
                if auto_pick(n).no_mog == deletedMog
                    found = n;
                    break;
                end
            end
            if found ~= -1
                ind=1:numel(auto_pick);
                ind = ind~=found;
                auto_pick = auto_pick(ind);
            end
        end
        
        hmodel.updateListMog()
        updateMogInfo()
        updateModelInfo()
    end
end