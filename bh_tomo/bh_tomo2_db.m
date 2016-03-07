function bh_tomo2_db(varargin)
%BH_TOMO2_DB

db_file = '';
auto_pick = [];
saved = true;

width = 1000;
height = 800;
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
hborehole.FontSize = 11;

%
% Create MOG panel
%
hmog = MogUI(hborehole,f,'Units','points');
hmog.FontSize = 11;

%
% Create Model panel
%
hmodel= ModelUI(hborehole,hmog,f,'Units','points');
hmodel.FontSize = 11;

%
% Create info panel
%
hinfo = uipanel(f,'Title','Infos','Units','points');
hinfo.FontSize = 12;
htextinfo = uicontrol('Style','text',...
    'Parent',hinfo,...
    'BackgroundColor',[1 1 1],...
    'Units','normalized',...
    'HorizontalAlignment','center',...
    'FontSize',12,...
    'Position',[0.05 0.05 0.9 0.9]);
htextinfo.String = {'','','Database: ','',...
    '0 Borehole(s)','0 MOG(s)','0 Model(s)','','0 traces'};

resizeUI()

%
% Add listeners
%
addlistener(hborehole,'boreholeAdded',@updateBHinfo);
addlistener(hborehole,'boreholeDeleted',@updateBHinfo);
addlistener(hborehole,'boreholeAdded',@hmog.updateBHlist);
addlistener(hborehole,'boreholeDeleted',@hmog.updateBHlist);
addlistener(hborehole,'boreholeEdited',@dbEdited);

addlistener(hmog,'mogAdded',@updateMogInfo);
addlistener(hmog,'mogDeleted',@updateMogInfo);
addlistener(hmog,'mogEdited',@dbEdited);

addlistener(hmodel,'modelAdded',@updateModelInfo);
addlistener(hmodel,'modelDeleted',@updateModelInfo);
addlistener(hmodel,'modelEdited',@dbEdited);


    function resizeUI(varargin)
        f.Visible = 'off';
        hBorder = 15;
        width = f.Position(3);
        height = f.Position(4);
        vBorderTop = 20;
        
        vSize = height/2 - vBorderTop;
        hsize = width/3 - hBorder;
        
        hborehole.Position = [hBorder height/2+hBorder/2 hsize vSize];
        hmog.Position = [2*hBorder+hsize height/2+hBorder/2 2*hsize vSize];
        
        hmodel.Position = [hBorder hBorder 2*hsize vSize];
        hinfo.Position = [width-hsize-hBorder hBorder hsize vSize];
        
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
        db_file = [rep,file];
        tmp = load(db_file);
        if isfield(tmp,'panels')
            warndlg(['File ',file,' is version 1 database, please convert first'])
            return
        end
        
        if ~isfield(tmp,'boreholes')
            errordlg(['File ',file,' not a bh_tomo database'],'File error')
            return
        end
        hborehole.boreholes = tmp.boreholes;
        hmog.mogs = tmp.mogs;
        hmodel.models = tmp.models;
        
        hmog.air = tmp.air;
        auto_pick = tmp.auto_pick; %#ok<SETNU>
        
        htextinfo.String{3} = ['Database: ',file];
        updateBHinfo()
        updateMogInfo()
        updateModelInfo()
        hmog.updateBHlist()
        saved = true;
    end
    function saveFile(varargin)
        if isempty(db_file)
            [file, rep] = uiputfile('*.mat','Save Database');
            if isequal(file,0)
                return
            end
            db_file = [rep,file];
            htextinfo.String{3} = ['Database: ',file];
        end
        names_mog = cell(1,length(hmog.mogs));
        for n=1:length(hmog.mogs)
            names_mog{n} = hmog.mogs(n).name;
        end
        models = hmodel.models; %#ok<NASGU>
        boreholes = hborehole.boreholes; %#ok<NASGU>
        mogs = hmog.mogs; %#ok<NASGU>
        air = hmog.air; %#ok<NASGU>
        save(db_file,'names_mog','mogs','air','boreholes','models','auto_pick')
        saved = true;
    end
    function saveFileAs(varargin)
        [file, rep] = uiputfile('*.mat','Save Database');
        if isequal(file,0)
            return
        end
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
        save(db_file,'names_mog','mogs','air','boreholes','models','auto_pick')
        saved = true;
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
end