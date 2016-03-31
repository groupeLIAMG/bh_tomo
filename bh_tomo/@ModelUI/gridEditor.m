function g = gridEditor(obj,no,varargin)

g = Grid.empty;
if nargin>=3
    if isa(varargin{1},'Grid')
        g = varargin{1};
    end
end
previousGrid = g.copy;

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
    'Callback',@closeWindow,...
    'Parent',f);

hcancel = uicontrol('Style','pushbutton',...
    'String','Cancel',...
    'Units','points',...
    'FontSize',fs,...
    'Callback',@cancel,...
    'Parent',f);

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
for n=model.mogs
    tmp = unique([mogs(n).data.Tx_x(mogs(n).in)' ...
        mogs(n).data.Tx_y(mogs(n).in)' ...
        mogs(n).data.Tx_z(mogs(n).in)'],'rows');
    plot3(haxes1,tmp(:,1),tmp(:,2),tmp(:,3),'b*')
    tmp = unique([mogs(n).data.Rx_x(mogs(n).in)' ...
        mogs(n).data.Rx_y(mogs(n).in)' ...
        mogs(n).data.Rx_z(mogs(n).in)'],'rows');
    plot3(haxes1,tmp(:,1),tmp(:,2),tmp(:,3),'g*')
    
    nTtDataPicked = nTtDataPicked + sum( mogs(n).tt_done(mogs(n).in)~=0 );
    nAmpDataPicked = nAmpDataPicked + sum( mogs(n).amp_done(mogs(n).in)~=0 );
    nData = nData + mogs(n).data.ntrace;
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
    switch type
        case '2D'
            [~,data] = Grid2DUI.buildGrid(data);
        case '2D+'
            [~,data] = Grid2DplusUI.buildGrid(data);
    end
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
    function cancel(varargin)
        g = previousGrid.copy;
        closeWindow()
    end
    function constraints(varargin)
        
        % add constraints defined at borehole if present
        
        for nn=model.boreholes
            % slowness
            if ~isempty(boreholes(nn).scont)
                switch type
                    case '2D'
                        coord = [boreholes(nn).scont.x boreholes(nn).scont.y boreholes(nn).scont.z];
                        [coord,~] = gUI.project(coord);
                        g.cont.slowness.data = [g.cont.slowness.data;
                            [coord(:,3) coord(:,1) boreholes(nn).scont.valeur boreholes(nn).scont.variance]];
                        g.cont.slowness.data = unique(g.cont.slowness.data,'rows');
                    case '2D+'
                        coord = [boreholes(nn).scont.x boreholes(nn).scont.y boreholes(nn).scont.z];
                        [coord,~] = gUI.project(coord);
                        g.cont.slowness.data = [g.cont.slowness.data;
                            [coord(:,3) coord(:,1) boreholes(nn).scont.valeur boreholes(nn).scont.variance]];
                        g.cont.slowness.data = unique(g.cont.slowness.data,'rows');
                    case '3D'
                        g.cont.slowness.data = [g.cont.slowness.data;
                            [boreholes(nn).scont.x ...
                            boreholes(nn).scont.y boreholes(nn).scont.z ...
                            boreholes(nn).scont.valeur boreholes(nn).scont.variance]];
                        g.cont.slowness.data = unique(g.cont.slowness.data,'rows');
                end
            end
            % attenuation
            if ~isempty(boreholes(nn).acont)         
                switch type
                    case '2D'
                        coord = [boreholes(nn).acont.x boreholes(nn).acont.y boreholes(nn).acont.z];
                        [coord,~] = gUI.project(coord);
                        g.cont.attenuation.data = [g.cont.attenuation.data;
                            [coord(:,3) coord(:,1) boreholes(nn).acont.valeur boreholes(nn).acont.variance]];
                        g.cont.attenuation.data = unique(g.cont.attenuation.data,'rows');
                    case '2D+'
                        coord = [boreholes(nn).acont.x boreholes(nn).acont.y boreholes(nn).acont.z];
                        [coord,~] = gUI.project(coord);
                        g.cont.attenuation.data = [g.cont.attenuation.data;
                            [coord(:,3) coord(:,1) boreholes(nn).acont.valeur boreholes(nn).acont.variance]];
                        g.cont.attenuation.data = unique(g.cont.attenuation.data,'rows');
                    case '3D'
                        g.cont.attenuation.data = [g.cont.attenuation.data;
                            [boreholes(nn).acont.x ...
                            boreholes(nn).acont.y boreholes(nn).acont.z ...
                            boreholes(nn).acont.valeur boreholes(nn).acont.variance]];
                        g.cont.attenuation.data = unique(g.cont.attenuation.data,'rows');
                end
            end
        end
        obj.constraintsEditor(g,gUI);
    end
    function gridEdited(varargin)
        hinfoTxt.String{3} = num2str(g.getNumberOfCells());
        obj.notify('modelEdited')
    end
    function data = prepareGridData()
        mogs = obj.mogUI.mogs;
        data.boreholes = [];
        for nn=model.boreholes
            data.boreholes = [data.boreholes boreholes(nn)];
        end
        data.in = [];
        data.Tx = [];
        data.Rx = [];
        data.TxCosDir = [];
        data.RxCosDir = [];
        data.Tx_Z_water = [];
        data.Rx_Z_water = [];
        for nn=model.mogs
            data.in = [data.in; mogs(nn).in'];
            data.Tx = [data.Tx; [mogs(nn).data.Tx_x' ...
                mogs(nn).data.Tx_y' ...
                mogs(nn).data.Tx_z'] ];
            data.Rx = [data.Rx; [mogs(nn).data.Rx_x' ...
                mogs(nn).data.Rx_y' ...
                mogs(nn).data.Rx_z'] ];
            data.TxCosDir = [data.TxCosDir; mogs(nn).TxCosDir];
            data.RxCosDir = [data.RxCosDir; mogs(nn).RxCosDir];
            
            bh = boreholes(mogs(nn).Tx);
            if ~isnan(bh.Z_water)
                x = interp1(bh.fdata(:,3), bh.fdata(:,1), bh.Z_water);
                y = interp1(bh.fdata(:,3), bh.fdata(:,2), bh.Z_water);
            else
                x = NaN;
                y = NaN;
            end
            data.Tx_Z_water = [data.Tx_Z_water; [x y bh.Z_water] ];

            bh = boreholes(mogs(nn).Rx);
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