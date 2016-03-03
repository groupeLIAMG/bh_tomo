function bh_tomo2_db
    %BH_TOMO2_DB
    
    width = 1000;
    height = 800;

    f = figure('Visible','off',...
        'Units','points',...
        'Position',[360 200 width height],...
        'Tag','fig_bh_tomo2_db',...
        'Name','bh_tomo_db',...
        'NumberTitle','off',...
        'ToolBar','none',...
        'SizeChangedFcn',@resizeUI);
    
    %
    % Create Borehole panel
    %
    hborehole = BoreholeUI(f,'Units','points');
    hborehole.FontSize = 11;
    
    %
    % Create MOG panel
    %
    hmog = uipanel(f,'Title','MOGs','Units','points');
    hmog.FontSize = 11;
  
    %
    % Create Model panel
    %
    hpanel= ModelUI(f,'Units','points');
    hpanel.FontSize = 11;
    
    %
    % Create info panel
    %
    hinfo = uipanel(f,'Title','Infos','Units','points');
    hinfo.FontSize = 11;
    htextinfo = uicontrol('Style','text',...
        'Parent',hinfo,...
        'BackgroundColor',[1 1 1],...
        'Units','normalized',...
        'HorizontalAlignment','center',...
        'FontSize',12,...
        'Position',[0.05 0.05 0.9 0.9]);
    htextinfo.String = {'','Database: ','','0 Borehole(s)','0 MOG(s)'};
    addlistener(hborehole,'boreholeAdded',@updateBHinfo);
    addlistener(hborehole,'boreholeDeleted',@updateBHinfo);
    
    resizeUI()
    
  
    function resizeUI(varargin)
        hBorder = 15;
        width = f.Position(3);
        height = f.Position(4);
        vBorderTop = 20;
        
        vSize = height/2 - vBorderTop;
        hsize = width/3 - hBorder;
        
        hborehole.Position = [hBorder height/2+hBorder/2 hsize vSize];
        hmog.Position = [2*hBorder+hsize height/2+hBorder/2 2*hsize vSize];
        
        hpanel.Position = [hBorder hBorder 2*hsize vSize];
        hinfo.Position = [width-hsize-hBorder hBorder hsize vSize];
        
        f.Visible = 'on';
    end

    function updateBHinfo(varargin)
        htextinfo.String{4} = [num2str(numel(hborehole.boreholes)),' Borehole(s)'];
    end
  
end