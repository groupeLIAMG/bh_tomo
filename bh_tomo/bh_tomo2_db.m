function bh_tomo2_db(varargin)
    %BH_TOMO2_DB
    
    db_file = '';
    air = AirShots.empty;
    auto_pick = [];
    
    width = 1000;
    height = 800;

    f = figure('Visible','off',...
        'Units','points',...
        'Position',[360 200 width height],...
        'Tag','fig_bh_tomo2_db',...
        'Name','bh_tomo_db',...
        'NumberTitle','off',...
        'ToolBar','none',...,
        'MenuBar','None',...
        'SizeChangedFcn',@resizeUI,...
        'CloseRequestFcn',@quitUI);
    
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
        '0 Borehole(s)','0 MOG(s)','0 Model(s)'};
    
    resizeUI()

    %
    % Add listeners
    %
    addlistener(hborehole,'boreholeAdded',@updateBHinfo);
    addlistener(hborehole,'boreholeDeleted',@updateBHinfo);
    addlistener(hborehole,'boreholeAdded',@hmog.updateBHlist);
    addlistener(hborehole,'boreholeDeleted',@hmog.updateBHlist);
    
    addlistener(hmodel,'modelAdded',@updateModelInfo);
    addlistener(hmodel,'modelDeleted',@updateModelInfo);
    
    addlistener(hmog,'mogAdded',@updateMogInfo);
    addlistener(hmog,'mogDeleted',@updateMogInfo);
  
    
  
    function resizeUI(varargin)
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

    function updateBHinfo(varargin)
        htextinfo.String{5} = [num2str(numel(hborehole.boreholes)),' Borehole(s)'];
    end
    function updateMogInfo(varargin)
        htextinfo.String{6} = [num2str(numel(hmog.mogs)),' MOG(s)'];
    end
    function updateModelInfo(varargin)
        htextinfo.String{7} = [num2str(numel(hmodel.models)),' Model(s)'];
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
        
        air = tmp.air; %#ok<SETNU>
        auto_pick = tmp.auto_pick; %#ok<SETNU>

        htextinfo.String{3} = ['Database: ',file];
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
        save(db_file,'names_mog','mogs','air','boreholes','models','auto_pick')
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
        save(db_file,'names_mog','mogs','air','boreholes','models','auto_pick')
    end
    function closeWindow(varargin)
        quitUI()
    end
    function quitUI(varargin)
        % TODO check for unsaved data
        delete(f)
    end
end