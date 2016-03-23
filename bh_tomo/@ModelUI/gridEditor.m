function g = gridEditor(obj,varargin)

g = Grid.empty;
no = varargin{1};
if nargin>=3
    if isa(varargin{2},'Grid')
        g = varargin{2};
    end
end

fs = 11;
if ispc
    fs = 9;
end
vScale = 1;
if ispc
    vScale = 0.81;
end

width = 1000;
height = 700*vScale;

% get screen size
su = get(groot,'Units');
set(groot,'Units','points')
scnsize = get(groot,'ScreenSize');
pos = [scnsize(3)/2-width/2 scnsize(4)/2-height/2 width height];
set(groot,'Units',su)       % Restore default root screen units

f = figure('Visible','off',...
    'Units','points',...
    'Position',pos,...
    'Tag','grid_editor',...
    'Name','Grid Editor',...
    'NumberTitle','off',...
    'ToolBar','none',...,
    'MenuBar','None',...
    'SizeChangedFcn',@resizeUI,...
    'CloseRequestFcn',@closeWindow);

hboreholes = uipanel(f,'Title','Boreholes',...
    'Units','points',...
    'FontSize',fs+1);
haxes1 = axes('Units','points','Parent',hboreholes);
hbhView = uicontrol('Style','popupmenu',...
    'String',{'XY Plane','XZ Plane','YZ Plane','3D View'},...
    'Value',4,...
    'Units','points',...
    'FontSize',fs,...
    'Callback',@bhView,...
    'Parent',hboreholes);

hdone = uicontrol('Style','pushbutton',...
    'String','Done',...
    'Units','points',...
    'FontSize',fs,...
    'Callback',@saveAndQuit,...
    'Parent',f);

hcancel = uicontrol('Style','pushbutton',...
    'String','Cancel',...
    'Units','points',...
    'FontSize',fs,...
    'Callback',@closeWindow,...
    'Parent',f);
if ~isempty(g)
    hcancel.Visible = 'off';  % we cannot cancel editing grid
end

hconstr = uicontrol('Style','pushbutton',...
    'String','Add/Edit Constraints',...
    'Units','points',...
    'FontSize',fs,...
    'Callback',@constraints,...
    'Parent',f);
hinfo = uipanel(f,'Title','Infos',...
    'Units','points',...
    'FontSize',fs+1);
hinfoTxt = uicontrol('Style','text',...
    'String',{'','Number of Cells','','','Number of Data','','','',''},...
    'Units','normalized',...
    'FontSize',fs,...
    'HorizontalAlignment','center',...
    'BackgroundColor',[1 1 1],...
    'Position', [0.05 0.05 0.9 0.9],...
    'Parent',hinfo);

% Find the number of boreholes in model
model = obj.models(no);
nBH = length(model.boreholes);
% plot boreholes
boreholes = obj.bhUI.boreholes;
plot3(haxes1,boreholes(model.boreholes(1)).fdata(:,1),...
    boreholes(model.boreholes(1)).fdata(:,2),...
    boreholes(model.boreholes(1)).fdata(:,3),'r','LineWidth',2)
hold(haxes1, 'on')
plot3(haxes1,boreholes(model.boreholes(1)).X,...
    boreholes(model.boreholes(1)).Y,...
    boreholes(model.boreholes(1)).Z_surf,'ro')
text(boreholes(model.boreholes(1)).X,...
    boreholes(model.boreholes(1)).Y,...
    boreholes(model.boreholes(1)).Z_surf,...
    boreholes(model.boreholes(1)).name,...
    'Parent',haxes1,...
    'FontSize',fs)
for n=2:length(model.boreholes)
    plot3(haxes1,boreholes(model.boreholes(n)).fdata(:,1),...
        boreholes(model.boreholes(n)).fdata(:,2),...
        boreholes(model.boreholes(n)).fdata(:,3),'r','LineWidth',2)
    plot3(haxes1,boreholes(model.boreholes(n)).X,...
        boreholes(model.boreholes(n)).Y,...
        boreholes(model.boreholes(n)).Z_surf,'ro')
    text(boreholes(model.boreholes(n)).X,...
        boreholes(model.boreholes(n)).Y,...
        boreholes(model.boreholes(n)).Z_surf,...
        boreholes(model.boreholes(n)).name,...
        'Parent',haxes1,...
        'FontSize',fs)
end
mogs = obj.mogUI.mogs;
nData = 0;
nTtDataPicked = 0;
nAmpDataPicked = 0;
for n=1:length(model.mogs)
    tmp = unique([mogs(model.mogs(n)).data.Tx_x(mogs(model.mogs(n)).in)' ...
        mogs(model.mogs(n)).data.Tx_y(mogs(model.mogs(n)).in)' ...
        mogs(model.mogs(n)).data.Tx_z(mogs(model.mogs(n)).in)'],'rows');
    plot3(haxes1,tmp(:,1),tmp(:,2),tmp(:,3),'*')
    tmp = unique([mogs(model.mogs(n)).data.Rx_x(mogs(model.mogs(n)).in)' ...
        mogs(model.mogs(n)).data.Rx_y(mogs(model.mogs(n)).in)' ...
        mogs(model.mogs(n)).data.Rx_z(mogs(model.mogs(n)).in)'],'rows');
    plot3(haxes1,tmp(:,1),tmp(:,2),tmp(:,3),'g*')
    
    nTtDataPicked = nTtDataPicked + sum( mogs(model.mogs(n)).tt_done(mogs(model.mogs(n)).in)~=0 );
    nAmpDataPicked = nAmpDataPicked + sum( mogs(model.mogs(n)).amp_done(mogs(model.mogs(n)).in)~=0 );
    nData = nData + mogs(model.mogs(n)).data.ntrace;
end
grid(haxes1,'on')
hold(haxes1, 'off')
xlabel(haxes1,'X [m]','FontSize',fs)
ylabel(haxes1,'Y [m]','FontSize',fs)
zlabel(haxes1,'Elevation [m]','FontSize',fs)
set(haxes1,'DataAspectRatio',[1 1 1])

haxesl = axes('Units','points','Parent',hboreholes);
plot(haxesl,1,1,'*')
hold(haxesl, 'on')
plot(haxesl,1,2,'g*')
text(1.5,1,'Tx','Parent',haxesl,'FontSize',fs)
text(1.5,2,'Rx','Parent',haxesl,'FontSize',fs)
hold(haxesl, 'off')
haxesl.XTickLabel='';
haxesl.XTick=[];
haxesl.YTickLabel='';
haxesl.YTick=[];
haxesl.YLim = [0 3];
haxesl.XLim = [0.5 2];

hinfoTxt.String{6} = num2str(nData);
hinfoTxt.String{8} = [num2str(nTtDataPicked),' Traveltimes Picked'];
hinfoTxt.String{9} = [num2str(nAmpDataPicked),' Amplitudes Picked'];
f.Visible = 'on';

data = prepareGridData();

if isempty(g)
    if nBH==2
        % 2D grid
        type = '2D';
        [g,data] = Grid2DUI.buildGrid(data);
    else
        % find number of unique Tx-Rx borehole pairs
        pairs = zeros(length(model.mogs),2);
        for nm=1:length(model.mogs)
            pairs(nm,:) = sort( [mogs(model.mogs(nm)).Tx mogs(model.mogs(nm)).Rx] );
        end
        pairs = unique(pairs,'rows');
        % do we have a borehole in more than 2 pairs
        moreThanTwo = false;
        for nb=1:length(model.boreholes)
            if sum( model.boreholes(nb) == pairs(:) )>2
                moreThanTwo = true;
            end
        end
        
        if moreThanTwo
            % we have to build a 3D grid
            type = '3D';
            g = Grid3DUI.buildGrid(data);
        else
            % we can make a "super panel", check if OK with user
            button = questdlg('Do you want to build a 3D Model or project boreholes in a 2D plane?',...
                'Choose Type of Grid','3D Model','2D Projection','2D Projection');
            switch button
                case '3D Model'
                    type = '3D';
                    g = Grid3DUI.buildGrid(data);
                case '2D Projection'
                    type = '2D+';
                    [g,data] = Grid2DplusUI.buildGrid(data);
            end
        end
    end
else
    type = g.type;
end

switch type
    case '2D'
        gUI = Grid2DUI(f,type,fs,g,data);
    case '2D+'
        gUI = Grid2DplusUI(f,type,fs,g,data);
    case '3D'
        gUI = Grid3DUI(f,type,fs,g,data);
end

hinfoTxt.String{3} = num2str(g.getNumberOfCells());

resizeUI()
addlistener(gUI,'gridEdited',@gridEdited);

uiwait(f)


    function resizeUI(varargin)
        f.Visible = 'off';
        
        width = f.Position(3);
        height = f.Position(4);
        hBorder = 15;
        vBorderTop = 20*vScale;
        
        vFac = 0.8*height/(700*vScale);
        if vFac<1
            vFac = 1;
        end
        if ispc
            vFac = vFac*vScale;
        end
        vSize = 22*vFac;
        vSpace = 5*vFac;
        
        vSize1 = 3*height/5 - vBorderTop;
        vSize2 = 2*height/5 - vBorderTop;
        hSize = width/2 - hBorder;
        
        hboreholes.Position = [hBorder hBorder hSize vSize1];
        haxesl.Position = [hSize-3.5*vSize vSize1-2.7*vSize 3*vSize 2*vSize];
        
        hSize1 = hSize/3 - hBorder;
        hdone.Position = [hBorder 2*hBorder+vSize1 hSize1 vSize];
        hcancel.Position = [hBorder 2*hBorder+vSize1+vSpace+vSize hSize1 vSize];
        hconstr.Position = [hBorder 2*hBorder+vSize1+2*vSpace+2*vSize hSize1 vSize];
        hinfo.Position = [hBorder 2*hBorder+vSize1+3*vSpace+3*vSize hSize1 vSize2-3*vSpace-3*vSize];
        
        axBorder = haxes1.Position(2);
        haxes1.Position = [axBorder axBorder hSize-1.5*axBorder vSize1-2*axBorder];
        ext = hbhView.Extent;
        hbhView.Position = [axBorder vSize1-2*vSize 2*ext(3) vSize];
        
        gUI.Position = [2*hBorder+hSize1 2*hBorder+vSize1 hSize-hSize1-hBorder vSize2];
        
        gUI.axPanel.Position = [hSize+2*hBorder hBorder hSize-hBorder height-2*vBorderTop+hBorder];
        
        f.Visible = 'on';
    end
    function closeWindow(varargin)
        delete(f)
    end
    function bhView(varargin)
        switch hbhView.Value
            case 1
                view(haxes1,2)
            case 2
                view(haxes1,0,0)
            case 3
                view(haxes1,90,0)
            case 4
                view(haxes1,3)
        end
    end
    function saveAndQuit(varargin)
        
        closeWindow()
    end
    function constraints(varargin)
    end
    function gridEdited(varargin)
        hinfoTxt.String{3} = num2str(g.getNumberOfCells());
    end
    
    function data = prepareGridData()
        mogs = obj.mogUI.mogs;
        data.boreholes = [];
        for n=1:length(model.boreholes)
            data.boreholes = [data.boreholes boreholes(model.boreholes(n))];
        end
        data.in = [];
        data.Tx = [];
        data.Rx = [];
        data.TxCosDir = [];
        data.RxCosDir = [];
        data.Tx_Z_water = [];
        data.Rx_Z_water = [];
        for nm=1:length(model.mogs)
            data.in = [data.in; mogs(model.mogs(nm)).in'];
            data.Tx = [data.Tx; [mogs(model.mogs(nm)).data.Tx_x' ...
                mogs(model.mogs(nm)).data.Tx_y' ...
                mogs(model.mogs(nm)).data.Tx_z'] ];
            data.Rx = [data.Rx; [mogs(model.mogs(nm)).data.Rx_x' ...
                mogs(model.mogs(nm)).data.Rx_y' ...
                mogs(model.mogs(nm)).data.Rx_z'] ];
            data.TxCosDir = [data.TxCosDir; mogs(model.mogs(nm)).TxCosDir];
            data.RxCosDir = [data.RxCosDir; mogs(model.mogs(nm)).RxCosDir];
            
            bh = boreholes(mogs(model.mogs(nm)).Tx);
            if ~isnan(bh.Z_water)
                x = interp1(bh.fdata(:,3), bh.fdata(:,1), bh.Z_water);
                y = interp1(bh.fdata(:,3), bh.fdata(:,2), bh.Z_water);
            else
                x = NaN;
                y = NaN;
            end
            data.Tx_Z_water = [data.Tx_Z_water; [x y bh.Z_water] ];

            bh = boreholes(mogs(model.mogs(nm)).Rx);
            if ~isnan(bh.Z_water)
                x = interp1(bh.fdata(:,3), bh.fdata(:,1), bh.Z_water);
                y = interp1(bh.fdata(:,3), bh.fdata(:,2), bh.Z_water);
            else
                x = NaN;
                y = NaN;
            end
            data.Rx_Z_water = [data.Rx_Z_water; [x y bh.Z_water] ];
        end
        data.in = logical(data.in);
    end
end