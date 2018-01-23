function bh_tomo2_tlinv( varargin )
%BH_TOMO2_TLINV Perform time-lapse tomographic inversion

rep='';
file='';
if nargin>=2
    rep=varargin{1};
    file=varargin{2};
end

modelNo = [];
model = [];
saved = true;
previousTypeData = 1;
previousTypeInv = 1;
tomo = [];
param = [];
gridViewer = [];

cminAmp = 1;
cmaxAmp = 3;
cminTT = 0.06;
cmaxTT = 0.12;

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

width = 1400;
height = 875*vScale;
% get screen size
su = get(groot,'Units');
set(groot,'Units','points')
scnsize = get(groot,'ScreenSize');
pos = [scnsize(3)/2-width/2 scnsize(4)/2-height/2 width height];
set(groot,'Units',su)       % Restore default root screen units

f = figure('Visible','off',...
    'Units','points',...
    'Position',pos,...
    'Tag','fig_bh_tomo2_tlinv',...
    'Name','bh_tomo_tlinv',...
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
uimenu(hmenu,'Label','Close',...
    'Separator','on',...
    'Accelerator','W',...
    'Callback',@closeWindow);
hresultsMenu = uimenu(f,'Label','Results');
uimenu(hresultsMenu,'Label','Export ...',...
    'Accelerator','E',...
    'Callback',@exportTomo);
uimenu(hresultsMenu,'Label','Tomograms',...
    'Accelerator','T',...
    'Callback',@showTomo);
uimenu(hresultsMenu,'Label','Difference Tomogram',...
    'Accelerator','G',...
    'Callback',@showDiff);
uimenu(hresultsMenu,'Label','Departure From Ref. Tomog.',...
    'Accelerator','H',...
    'Callback',@showDep);
uimenu(hresultsMenu,'Label','Rays',...
    'Accelerator','R',...
    'Callback',@showRays);
uimenu(hresultsMenu,'Label','Ray Density',...
    'Accelerator','D',...
    'Callback',@showRayDensity);
uimenu(hresultsMenu,'Label','Residuals',...
    'Accelerator','L',...
    'Callback',@showResiduals);


%
% Main Panels
%

pdata = uipanel(f,'Title','Data',...
    'Units','points',...
    'FontSize',fs+1);
pprevious = uipanel(f,'Title','Previous Inversion',...
    'Units','points',...
    'FontSize',fs+1);
pinv = uipanel(f,'Title','Inversion Parameters',...
    'Units','points',...
    'FontSize',fs+1);
pfig = uipanel(f,'Title','Figures',...
    'Units','points',...
    'FontSize',fs+1);

hmessage = uicontrol('Style','text',...
    'ForegroundColor','red',...
    'BackgroundColor','white',...
    'FontSize',fs+1,...
    'Units','points',...
    'Parent',f);

nLines=8;
vSizeTot = nLines*22 + (nLines+1)*5;
vSize = 22/vSizeTot;
vSpace = 5/vSizeTot;

% Data

hmodelName = uicontrol('Style','text',...
    'String','Model: ',...
    'ForegroundColor','red',...
    'FontSize',fs+1,...
    'HorizontalAlignment','left',...
    'Units','normalized',...
    'Position',[0.05 7*vSize+8*vSpace 0.9 vSize],...
    'Parent',pdata);
htypeData = uicontrol('Style','popupmenu',...
    'String',{'Traveltime','Amplitude - Peak-to-Peak','Amplitude - Centroid Frequency'},...
    'FontSize',fs,...
    'Units','normalized',...
    'Position',[0.05 5.5*vSize+6.5*vSpace 0.5 vSize],...
    'Callback',@changeTypeData,...
    'Parent',pdata);
uicontrol('Style','pushbutton',...
    'String',['Show ',char(916),' data'],...
    'FontSize',fs,...
    'Units','normalized',...
    'Position',[0.6 5.5*vSize+6.5*vSpace 0.3, vSize],...
    'Callback',@showDelta,...
    'Parent',pdata);
uicontrol('Style','text',...
    'String','Baseline Data',...
    'FontSize',fs,...
    'Units','normalized',...
    'Position',[0.03 4*vSize+5*vSpace 0.45 vSize],...
    'Parent',pdata);
uicontrol('Style','text',...
    'String','Repeat Data',...
    'FontSize',fs,...
    'Units','normalized',...
    'Position',[0.51 4*vSize+5*vSpace 0.45 vSize],...
    'Parent',pdata);
hlistBaseline = uicontrol('Style','listbox',...
    'FontSize',fs,...
    'Units','normalized',...
    'HorizontalAlignment','center',...
    'Position',[0.03 vSpace 0.45 5*vSize+5*vSpace],...
    'Parent',pdata);
hlistRepeat = uicontrol('Style','listbox',...
    'FontSize',fs,...
    'Units','normalized',...
    'HorizontalAlignment','center',...
    'Position',[0.51 vSpace 0.45 5*vSize+5*vSpace],...
    'Parent',pdata);

% Previous Inversions

nLines=1;
vSizeTot = nLines*22 + (nLines+1)*5;
vSize = 22/vSizeTot;
vSpace = 5/vSizeTot;

hpreviousInv = uicontrol('Style','popupmenu',...
    'String','-',...
    'FontSize',fs,...
    'Units','normalized',...
    'Position',[0.025 vSpace 0.6 vSize],...
    'Parent',pprevious);
uicontrol('Style','pushbutton',...
    'String','Load',...
    'FontSize',fs,...
    'Units','normalized',...
    'Position',[0.65 vSpace 0.15 vSize],...
    'Callback',@loadPrevious,...
    'Parent',pprevious);
uicontrol('Style','pushbutton',...
    'String','Delete',...
    'FontSize',fs,...
    'Units','normalized',...
    'Position',[0.825 vSpace 0.15 vSize],...
    'Callback',@deletePrevious,...
    'Parent',pprevious);

% Inversion params

nLines=16;
vSizeTot = nLines*22 + (nLines+1)*5;
vSize = 22/vSizeTot;
vSpace = 5/vSizeTot;

uicontrol('Style','pushbutton',...
    'String','GO',...
    'FontSize',fs+1,...
    'FontWeight','bold',...
    'Units','normalized',...
    'Position',[0.3 1.5*vSpace 0.4 vSize],...
    'Callback',@doInv,...
    'Parent',pinv);

psimult = uipanel(pinv,'Title','Parameters - Simultaneous Inversion',...
    'Units','normalized',...
    'Position',[0.025 2*vSpace+1.25*vSize 0.95 8*vSize+8*vSpace],...
    'FontSize',fs,...
    'Visible','on');
pdiff = uipanel(pinv,'Title','Parameters - Difference Inversion',...
    'Units','normalized',...
    'Position',[0.025 2*vSpace+1.25*vSize 0.95 8*vSize+8*vSpace],...
    'FontSize',fs,...
    'Visible','off');

pweight = uipanel(pinv,'Title','Weight - Reservoir Cells',...
    'Units','normalized',...
    'Position',[0.025 11*vSpace+9.7*vSize 0.95 1.5*vSize+2*vSpace],...
    'FontSize',fs,...
    'Visible','on');
pref = uipanel(pinv,'Title','Reference Tomogram',...
    'Units','normalized',...
    'Position',[0.025 13*vSpace+11.75*vSize 0.95 1.5*vSize+2*vSpace],...
    'FontSize',fs,...
    'Visible','on');


uicontrol('Style','text',...
    'String','Algorithm',...
    'FontSize',fs,...
    'Units','normalized',...
    'HorizontalAlignment','right',...
    'Position',[0.1 16*vSpace+14.75*vSize 0.3 vSize],...
    'Parent',pinv);
htypeInv = uicontrol('Style','popupmenu',...
    'String',{'Simultaneous Inversion',...
    'Difference Inversion',...
    ['Straight-Ray ',char(916),'t Inversion']},...
    'FontSize',fs,...
    'Units','normalized',...
    'Position',[0.45 16*vSpace+14.75*vSize 0.4 vSize],...
    'Callback',@changeTypeInv,...
    'Parent',pinv);
uicontrol('Style','text',...
    'String','Number of Iterations',...
    'FontSize',fs,...
    'Units','normalized',...
    'HorizontalAlignment','right',...
    'Position',[0.1 15*vSpace+13.75*vSize 0.3 vSize],...
    'Parent',pinv);
hnumIt = uicontrol('Style','edit',...
    'String',5,...
    'FontSize',fs,...
    'Units','normalized',...
    'Position',[0.45 15*vSpace+13.75*vSize 0.2 vSize],...
    'Callback',@changeTypeInv,...
    'Parent',pinv);

nLines=1;
vSizeTot = nLines*22 + (nLines+1)*5;
vSize = 22/vSizeTot;
vSpace = 5/vSizeTot;

hrefTomo = uicontrol('Style','popupmenu',...
    'String',{'-'},...
    'FontSize',fs,...
    'Units','normalized',...
    'Position',[0.05 vSpace 0.5 vSize],...
    'Parent',pref);
uicontrol('Style','pushbutton',...
    'String','View',...
    'FontSize',fs,...
    'Units','normalized',...
    'HorizontalAlignment','right',...
    'Position',[0.65 vSpace 0.2 vSize],...
    'Callback',@showRef,...
    'Parent',pref);



nLines=8;
vSizeTot = nLines*22 + (nLines+1)*5;
vSize = 22/vSizeTot;
vSpace = 5/vSizeTot;

uicontrol('Style','text',...
    'String',[char(945),' : '],...
    'FontSize',fs,...
    'HorizontalAlignment','right',...
    'Units','normalized',...
    'Position',[0.05 8*vSpace+7*vSize 0.4 vSize],...
    'Parent',psimult);
halphaS = uicontrol('Style','edit',...
    'String','1',...
    'FontSize',fs,...
    'Units','normalized',...
    'Position',[0.5 8*vSpace+7*vSize 0.2 vSize],...
    'Parent',psimult);
uicontrol('Style','text',...
    'String',[char(946),' : '],...
    'FontSize',fs,...
    'HorizontalAlignment','right',...
    'Units','normalized',...
    'Position',[0.05 7*vSpace+6*vSize 0.4 vSize],...
    'Parent',psimult);
hbetaS = uicontrol('Style','edit',...
    'String','75',...
    'FontSize',fs,...
    'Units','normalized',...
    'Position',[0.5 7*vSpace+6*vSize 0.2 vSize],...
    'Parent',psimult);
uicontrol('Style','text',...
    'String',[char(955),' : '],...
    'FontSize',fs,...
    'HorizontalAlignment','right',...
    'Units','normalized',...
    'Position',[0.05 6*vSpace+5*vSize 0.4 vSize],...
    'Parent',psimult);
hlambdaS = uicontrol('Style','edit',...
    'String','0.2',...
    'FontSize',fs,...
    'Units','normalized',...
    'Position',[0.5 6*vSpace+5*vSize 0.2 vSize],...
    'Parent',psimult);
uicontrol('Style','text',...
    'String','Smoothing Weight x : ',...
    'FontSize',fs,...
    'HorizontalAlignment','right',...
    'Units','normalized',...
    'Position',[0.05 5*vSpace+4*vSize 0.4 vSize],...
    'Parent',psimult);
hmuxS = uicontrol('Style','edit',...
    'String','1',...
    'FontSize',fs,...
    'Units','normalized',...
    'Position',[0.5 5*vSpace+4*vSize 0.2 vSize],...
    'Parent',psimult);
uicontrol('Style','text',...
    'String','Smoothing Weight y : ',...
    'FontSize',fs,...
    'HorizontalAlignment','right',...
    'Units','normalized',...
    'Position',[0.05 4*vSpace+3*vSize 0.4 vSize],...
    'Parent',psimult);
hmuyS = uicontrol('Style','edit',...
    'String','1',...
    'FontSize',fs,...
    'Units','normalized',...
    'Position',[0.5 4*vSpace+3*vSize 0.2 vSize],...
    'Parent',psimult);
uicontrol('Style','text',...
    'String','Smoothing Weight z : ',...
    'FontSize',fs,...
    'HorizontalAlignment','right',...
    'Units','normalized',...
    'Position',[0.05 3*vSpace+2*vSize 0.4 vSize],...
    'Parent',psimult);
hmuzS = uicontrol('Style','edit',...
    'String','1',...
    'FontSize',fs,...
    'Units','normalized',...
    'Position',[0.5 3*vSpace+2*vSize 0.2 vSize],...
    'Parent',psimult);
uicontrol('Style','text',...
    'String','Step Damping Factor : ',...
    'FontSize',fs,...
    'HorizontalAlignment','right',...
    'Units','normalized',...
    'Position',[0.05 2*vSpace+vSize 0.4 vSize],...
    'Parent',psimult);
hdampS = uicontrol('Style','edit',...
    'String','0.2',...
    'FontSize',fs,...
    'Units','normalized',...
    'Position',[0.5 2*vSpace+vSize 0.2 vSize],...
    'Parent',psimult);
hraysStext = uicontrol('Style','text',...
    'String','Rays : ',...
    'FontSize',fs,...
    'HorizontalAlignment','right',...
    'Units','normalized',...
    'Position',[0.05 vSpace 0.4 vSize],...
    'Visible','off',...
    'Parent',psimult);
hraysS = uicontrol('Style','popupmenu',...
    'String',{'-'},...
    'FontSize',fs,...
    'Units','normalized',...
    'Position',[0.5 vSpace 0.45 vSize],...
    'Visible','off',...
    'Parent',psimult);


uicontrol('Style','text',...
    'String',[char(946),' : '],...
    'FontSize',fs,...
    'HorizontalAlignment','right',...
    'Units','normalized',...
    'Position',[0.05 8*vSpace+7*vSize 0.4 vSize],...
    'Parent',pdiff);
hbetaD = uicontrol('Style','edit',...
    'String','75',...
    'FontSize',fs,...
    'Units','normalized',...
    'Position',[0.5 8*vSpace+7*vSize 0.2 vSize],...
    'Parent',pdiff);
uicontrol('Style','text',...
    'String',[char(955),' : '],...
    'FontSize',fs,...
    'HorizontalAlignment','right',...
    'Units','normalized',...
    'Position',[0.05 7*vSpace+6*vSize 0.4 vSize],...
    'Parent',pdiff);
hlambdaD = uicontrol('Style','edit',...
    'String','0.2',...
    'FontSize',fs,...
    'Units','normalized',...
    'Position',[0.5 7*vSpace+6*vSize 0.2 vSize],...
    'Parent',pdiff);
uicontrol('Style','text',...
    'String','Smoothing Weight x : ',...
    'FontSize',fs,...
    'HorizontalAlignment','right',...
    'Units','normalized',...
    'Position',[0.05 6*vSpace+5*vSize 0.4 vSize],...
    'Parent',pdiff);
hmuxD = uicontrol('Style','edit',...
    'String','1',...
    'FontSize',fs,...
    'Units','normalized',...
    'Position',[0.5 6*vSpace+5*vSize 0.2 vSize],...
    'Parent',pdiff);
uicontrol('Style','text',...
    'String','Smoothing Weight y : ',...
    'FontSize',fs,...
    'HorizontalAlignment','right',...
    'Units','normalized',...
    'Position',[0.05 5*vSpace+4*vSize 0.4 vSize],...
    'Parent',pdiff);
hmuyD = uicontrol('Style','edit',...
    'String','1',...
    'FontSize',fs,...
    'Units','normalized',...
    'Position',[0.5 5*vSpace+4*vSize 0.2 vSize],...
    'Parent',pdiff);
uicontrol('Style','text',...
    'String','Smoothing Weight z : ',...
    'FontSize',fs,...
    'HorizontalAlignment','right',...
    'Units','normalized',...
    'Position',[0.05 4*vSpace+3*vSize 0.4 vSize],...
    'Parent',pdiff);
hmuzD = uicontrol('Style','edit',...
    'String','1',...
    'FontSize',fs,...
    'Units','normalized',...
    'Position',[0.5 4*vSpace+3*vSize 0.2 vSize],...
    'Parent',pdiff);
uicontrol('Style','text',...
    'String','Step Damping Factor : ',...
    'FontSize',fs,...
    'HorizontalAlignment','right',...
    'Units','normalized',...
    'Position',[0.05 3*vSpace+2*vSize 0.4 vSize],...
    'Parent',pdiff);
hdampD = uicontrol('Style','edit',...
    'String','0.2',...
    'FontSize',fs,...
    'Units','normalized',...
    'Position',[0.5 3*vSpace+2*vSize 0.2 vSize],...
    'Parent',pdiff);
hraysDtext = uicontrol('Style','text',...
    'String','Rays : ',...
    'FontSize',fs,...
    'HorizontalAlignment','right',...
    'Units','normalized',...
    'Position',[0.05 2*vSpace+vSize 0.4 vSize],...
    'Visible','off',...
    'Parent',pdiff);
hraysD = uicontrol('Style','popupmenu',...
    'String',{'-'},...
    'FontSize',fs,...
    'Units','normalized',...
    'Position',[0.5 2*vSpace+vSize 0.45 vSize],...
    'Visible','off',...
    'Parent',pdiff);


nLines=1;
vSizeTot = nLines*22 + (nLines+1)*5;
vSize = 22/vSizeTot;
vSpace = 5/vSizeTot;

hreservoir = uicontrol('Style','checkbox',...
    'String','Apply',...
    'FontSize',fs,...
    'Units','normalized',...
    'Position',[0.1 vSpace 0.15 vSize],...
    'Parent',pweight);
uicontrol('Style','text',...
    'String','Value: ',...
    'FontSize',fs,...
    'HorizontalAlignment','right',...
    'Units','normalized',...
    'Position',[0.3 vSpace 0.3 vSize],...
    'Parent',pweight);
hresWeight = uicontrol('Style','edit',...
    'String','0.0001',...
    'FontSize',fs,...
    'Units','normalized',...
    'Position',[0.65 vSpace 0.2 vSize],...
    'Parent',pweight);



hcolorbar = uicontrol('Style','checkbox',...
    'String','Set Color Limits',...
    'FontSize',fs,...
    'Units','normalized',...
    'Position',[0.02 vSpace 0.12 vSize],...
    'Callback',@doColorbar,...
    'Parent',pfig);
uicontrol('Style','text',...
    'String','Min: ',...
    'HorizontalAlignment','right',...
    'FontSize',fs,...
    'Units','normalized',...
    'Position',[0.15 vSpace 0.04 vSize],...
    'Parent',pfig)
hcmin = uicontrol('Style','edit',...
    'String','0.06',...
    'FontSize',fs,...
    'Units','normalized',...
    'Position',[0.2 vSpace 0.07 vSize],...
    'Callback',@setClim,...
    'Parent',pfig);
uicontrol('Style','text',...
    'String','Max: ',...
    'HorizontalAlignment','right',...
    'FontSize',fs,...
    'Units','normalized',...
    'Position',[0.28 vSpace 0.04 vSize],...
    'Parent',pfig)
hcmax = uicontrol('Style','edit',...
    'String','0.12',...
    'FontSize',fs,...
    'Units','normalized',...
    'Position',[0.33 vSpace 0.07 vSize],...
    'Callback',@setClim,...
    'Parent',pfig);

m = {'plasma','magma','inferno','viridis','cmr','polarmap','parula','jet',...
    'hsv','hot','cool','autumn','spring','winter',...
    'summer','gray','bone','copper','pink','prism','flag','colorcube','lines'};

hcmap = uicontrol('Style','popupmenu',...
    'String',m,...
    'FontSize',fs,...
    'Units','normalized',...
    'Position',[0.41 vSpace 0.12 vSize],...
    'Callback',@doMap,...
    'Parent',pfig);
hsliderText = uicontrol('Style','text',...
    'String','Y Plane',...
    'FontSize',fs,...
    'Units','normalized',...
    'Position',[0.54 vSpace 0.07 vSize],...
    'Visible','off',...
    'Parent',pfig);


haxes1 = axes('Units','points','Parent',f,'Visible','on');
haxes2 = axes('Units','points','Parent',f,'Visible','on');

f.Visible = 'on';


    function resizeUI(varargin)
        f.Visible = 'off';
        hBorder = 15;
        width = f.Position(3);
        height = f.Position(4);
        
        vFac = 1;
        if ispc
            vFac = 0.81*vFac;
        end
        vSize = 22*vFac;
        vSpace = 5*vFac;
        vBorder = 10*vFac;
        
        axBorder = 60;%haxes1.Position(2);
        rBorder = 2*hBorder;
        hSize = 420;
        
        vSize2 = 8*vSize+9*vSpace;
        vPos = height-2*vBorder-vSize2;
        pdata.Position = [hBorder vPos hSize vSize2];

        vSize2 = 1.5*vSize+3*vSpace;
        vPos = vPos-vBorder-vSize2;
        pprevious.Position = [hBorder vPos hSize vSize2];
        
        vSize2 = 17*vSize+18*vSpace;
        vPos = vPos-2*vBorder-vSize2;
        pinv.Position = [hBorder vPos hSize vSize2];
        
        vSize2 = height-2*axBorder-4*vSize;
        hSize2 = width-(hBorder+hSize+axBorder+rBorder);
        hSize3 = (hSize2-rBorder)/2;

        hmessage.Position = [hBorder+hSize+axBorder height-2*vSize hSize2 vSize];
        pfig.Position = [hBorder+hSize+axBorder height-4.25*vSize hSize2 1.5*vSize+2*vSpace];
        
        haxes1.Position = [hBorder+hSize+axBorder axBorder hSize3 vSize2];
        haxes2.Position = [hBorder+hSize+axBorder+hSize3+rBorder axBorder hSize3 vSize2];

        f.Visible = 'on';
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
    function saveFile(varargin)
        if isempty(model)
            return
        end
        
        if isempty(tomo)
            % save current model even if empty tomo,
            %    in case previous inversion were removed
            load([rep,file],'models')
            models(modelNo).inv_res = model.inv_res; %#ok<STRNU>
            save([rep,file],'models','-append')
            saved = true;
            return
        end
        
        switch htypeData.Value
            case 1
                dType = '-vel';
            otherwise
                dType = '-att';
        end
        switch htypeInv.Value
            case 1
                iType = 'simult';
            case 2
                iType = 'diff';
            case 3
                iType = 'dt';
                
        end
        
        nameDefault = ['tomo(',tomo.date,')',dType,iType];  
        prompt = {'Inversion name:                                                               '};
        title = 'Save inversion results';
        nblines = 1;
        answer = myinputdlg(prompt,title,nblines,{nameDefault},'on');
        if ~isempty(answer)
            name=answer{1};
        end
        
        if isempty(model.tlinv_res)
            no_tlinv_res = 1;
        else
            flag=0;
            no_tlinv_res = 1+length(model.tlinv_res);
            for n=1:length(model.tlinv_res)
                if strcmp(model.tlinv_res(n).name, name)
                    no_tlinv_res = n;
                    flag = 1;
                    break;
                end
            end
            if flag==1
                answer=questdlg(['Overwrite ',name,'?']);
                if ~strcmp(answer,'Yes')
                    return
                end
            end
        end
        model.tlinv_res(no_tlinv_res).name = name;
        model.tlinv_res(no_tlinv_res).tomo = tomo;
        model.tlinv_res(no_tlinv_res).param = param;
        
        load([rep,file],'models')
        models(modelNo) = model; %#ok<NASGU>
        save([rep,file],'models','-append')
        saved = true;
        
        names = cell(numel(model.tlinv_res),1);
        for n=1:numel(model.tlinv_res)
            names{n} = model.tlinv_res(n).name;
        end
        hpreviousInv.String = names;
    end
    function openFile(varargin)
        [modelNo,file2,rep2] = chooseModel(rep,file);
        if isempty(modelNo)
            return
        end
        rep=rep2;
        file=file2;
        try
            load([rep,file],'models','mogs')
        catch ME
            errordlg(ME.message)
            return
        end
        model = models(modelNo);
        if isempty(model.grid)
            errordlg('Grid not created, inversion cannot be computed')
            return
        end
        
        %
        % Reset UI
        %
        mname = cell(1,numel(model.mogs));
        for n = 1:numel(model.mogs)
            mname{n} = [mogs(model.mogs(n)).name,' - ',mogs(model.mogs(n)).date];
        end
        hlistBaseline.String = mname;
        hlistBaseline.Max = numel(model.mogs);
        hlistRepeat.String = mname;
        hlistRepeat.Max = numel(model.mogs);
        hmodelName.String = model.name;
        
        if ~isempty(model.tlinv_res)
            nr = cell(numel(model.tlinv_res),1);
            for n=1:numel(model.tlinv_res)
                nr{n} = model.tlinv_res(n).name;
            end
        else
            nr = cell(1);
            nr{1} = '-';
        end
        hpreviousInv.String = nr;
        hpreviousInv.Value = 1;
        
        inv_name = cell(1,length(model.inv_res));
        for n=1:length(model.inv_res)
            inv_name{n} = model.inv_res(n).name;
        end
        hrefTomo.String = inv_name;
        hraysS.String = inv_name;
        hraysD.String = inv_name;
        
        htypeData.Value = 1;
        gridViewer = GridViewer(model.grid);
        if strcmp(model.grid.type, '3D')
            hsliderText.Visible = 'on';
            gridViewer.createSliders('Parent',pfig,...
                'Units','normalized','Visible','on');
            nLines=1;
            vSizeTot = nLines*22 + 2*5;
            vSize = 22/vSizeTot;
            vSpace = 5/vSizeTot;
            gridViewer.slider1.Position = [0.64 vSpace 0.15 vSize];
            gridViewer.slider2.Position = [0.80 vSpace 0.15 vSize];
        end
        
    end

    function doInv(varargin)
        if isempty(model)
            return
        end
        hmessage.String = 'Starting ...';
        cla(haxes1);cla(haxes2);
        cbh = findobj( 0, 'tag', 'Colorbar' );
        for i = 1: length(cbh)
            colorbar(cbh(i),'off')
        end

        param = [];
        param.db_file = [rep,file];
        param.saveInvData = 1;
        if ~isempty(model.grid.cont.ind_reservoir) && hreservoir.Value == 1
            param.ind_reservoir = model.grid.cont.ind_reservoir;
        else
            param.ind_reservoir = false(model.grid.getNumberOfCells(),1);
        end
        param.weight_reservoir = str2double(hresWeight.String);
        param.max_it = str2double(hnumIt.String);
        param.ref_inv_no = hrefTomo.Value;

        cmap = hcmap.String{hcmap.Value};
        clim = [];
        if htypeData.Value == 1
            param.typeData = 'tt';
            param.tomoAtt = 0;
            
            if model.inv_res(param.ref_inv_no).param.tomoAtt == 1
                errordlg('Cannot use attenuation model as reference for traveltime inversion')
                hmessage.String = '';
                return
            end
            if hcolorbar.Value==1
                clim = [cminTT cmaxTT];
            end
        elseif htypeData.Value == 2
            param.typeData = 'amp';
            param.tomoAtt = 1;
            
            if model.inv_res(param.ref_inv_no).param.tomoAtt == 0
                errordlg('Cannot use velocity model as reference for amplitude inversion')
                hmessage.String = '';
                return
            end
            if hcolorbar.Value==1
                clim = [cminAmp cmaxAmp];
            end
        else
            param.typeData = 'fce';
            param.tomoAtt = 1;
            
            if model.inv_res(param.ref_inv_no).param.tomoAtt == 0
                errordlg('Cannot use velocity model as reference for amplitude inversion')
                hmessage.String = '';
                return
            end
            if hcolorbar.Value==1
                clim = [cminAmp cmaxAmp];
            end
        end
        
        gh = {clim; cmap; haxes1; haxes2};

        if htypeInv.Value == 1
            % Simultaneous
            param.alpha = str2double(halphaS.String);
            param.beta = str2double(hbetaS.String);
            param.lambda = str2double(hlambdaS.String);
            param.mux = str2double(hmuxS.String);
            param.muy = str2double(hmuyS.String);
            param.muz = str2double(hmuzS.String);
            param.damping = str2double(hdampS.String);
            
            param.mog_no0 = hlistBaseline.Value;
            param.mog_no1 = hlistRepeat.Value;
            param.L_tomo_no = hraysS.Value;
            
            if param.tomoAtt == 1 && ~isfield(model.inv_res(param.L_tomo_no).tomo,'no_trace0')
                errordlg('Selected rays not obtained with simultaneous inversion')
                hmessage.String = '';
                return
            end
            
            tomo = invSimultaneous(model,param,hmessage,gh,gridViewer);
        elseif htypeInv.Value == 2
            % Difference
            param.beta = str2double(hbetaD.String);
            param.lambda = str2double(hlambdaD.String);
            param.mux = str2double(hmuxD.String);
            param.muy = str2double(hmuyD.String);
            param.muz = str2double(hmuzD.String);
            param.damping = str2double(hdampD.String);

            param.mog_no = hlistRepeat.Value;
            param.L_tomo_no = hraysD.Value;
            
            % Show Reference model
            if param.tomoAtt == 1
                gridViewer.plotTomo(model.inv_res(param.ref_inv_no).tomo.s,...
                    'Baseline Survey: Attenuation','Distance [m]','Elevation [m]',haxes1)
            else
                gridViewer.plotTomo(1./model.inv_res(param.ref_inv_no).tomo.s,...
                    'Baseline Survey: Velocity','Distance [m]','Elevation [m]',haxes1)
            end
            if ~isempty(clim), caxis(haxes1,clim), end
            colorbar('peer', haxes1)
            colormap(haxes1,cmap)
            drawnow
            
            tomo = invDifference(model,param,hmessage,gh,gridViewer);
        else
            warndlg('Not implemented yet')
        end
        saved = false;
        hmessage.String = '';
        
    end
    function changeTypeData(varargin)
        if isempty(model)
            warndlg('Data not Loaded')
            htypeData.Value = previousTypeData;
            return
        end
        if htypeData.Value ~= previousTypeData
            % type was really changed
            if previousTypeData == 1
                hraysS.Visible = 'on';
                hraysD.Visible = 'on';
                hraysStext.Visible = 'on';
                hraysDtext.Visible = 'on';
                
            elseif htypeData.Value == 1
                hraysS.Visible = 'off';
                hraysD.Visible = 'off';
                hraysStext.Visible = 'off';
                hraysDtext.Visible = 'off';
            end                
            previousTypeData = htypeData.Value;
        end
    end
    function changeTypeInv(varargin)
        if isempty(model)
            warndlg('Data not Loaded')
            htypeInv.Value = previousTypeInv;
            return
        end
        if htypeInv.Value ~= previousTypeInv
            % type was really changed
            if htypeInv.Value==1
                pdiff.Visible = 'off';
                psimult.Visible = 'on';
            else
                pdiff.Visible = 'on';
                psimult.Visible = 'off';
            end
            previousTypeInv = htypeInv.Value;
        end
    end

    function doColorbar(varargin)
        if hcolorbar.Value==1
            cl = [str2double(hcmin.String) str2double(hcmax.String)];
            haxes1.CLim = cl;
            haxes2.CLim = cl;
        else
            haxes1.CLimMode = 'auto';
            haxes2.CLimMode = 'auto';
        end
    end
    function setClim(src,varargin)
        if htypeData.Value==1
            if src==hcmin
                if str2double(src.String)>cmaxTT
                    warndlg({'Value should be smaller than max value',...
                        'Resetting value'})
                    src.String = num2str(cminTT);
                    return
                else
                    cminTT = str2double(src.String);
                end
            else
                if str2double(src.String)<cminTT
                    warndlg({'Value should be greater than min value',...
                        'Resetting value'})
                    src.String = num2str(cmaxTT);
                    return
                else
                    cmaxTT = str2double(src.String);
                end
            end
        else
            if src==hcmin
                if str2double(src.String)>cmaxAmp
                    warndlg({'Value should be smaller than max value',...
                        'Resetting value'})
                    src.String = num2str(cminAmp);
                    return
                else
                    cminAmp = str2double(src.String);
                end
            else
                if str2double(src.String)<cminAmp
                    warndlg({'Value should be greater than min value',...
                        'Resetting value'})
                    src.String = num2str(cmaxAmp);
                    return
                else
                    cmaxAmp = str2double(src.String);
                end
            end
        end
        if hcolorbar.Value==1
            doColorbar()
        end
    end
    function doMap(varargin)
        colormap(haxes1, hcmap.String{hcmap.Value})
        colormap(haxes2, hcmap.String{hcmap.Value})
    end
    function showRef(varargin)
        if isempty(model)
            warndlg('Data not Loaded')
            return
        end
        if isempty(model.inv_res)
            warndlg('No Inversion Results Found For This Model')
            return
        end
        nf=figure;
        ax=axes('Parent',nf);
        if model.inv_res(hrefTomo.Value).param.tomoAtt == 1
            gridViewer.plotTomo(model.inv_res(hrefTomo.Value).tomo.s,...
                'Reference Model: Attenuation','Distance [m]','Elevation [m]',ax)
        else
            gridViewer.plotTomo(1./model.inv_res(hrefTomo.Value).tomo.s,...
                'Reference Model: Velocity','Distance [m]','Elevation [m]',ax)
        end
        colorbar('peer', ax)
        colormap(ax,hcmap.String{hcmap.Value})
    end

    function exportTomo(varargin)
        
    end
    function showTomo(varargin)
        if isempty(tomo)
            return
        end
        if param.tomoAtt==0
            t1 = 1./tomo.s;
            t0 = 1./tomo.s0;
        else
            t1 = tomo.s;
            t0 = tomo.s0;
        end
        nf=figure;
        ax0=subplot(1,2,1,'Parent',nf);
        ax1=subplot(1,2,2,'Parent',nf);
        
        if strcmp(model.grid.type,'3D')==1
            gv = GridViewer(model.grid);
            gv.createSliders('Parent',nf);
            gv.plotTomo(t0,'Baseline','Distance [m]','Elevation [m]',ax0)
            gv.plotTomo(t1,'Repeat','Distance [m]','Elevation [m]',ax1)
        else
            gridViewer.plotTomo(t0,'Baseline','Distance [m]','Elevation [m]',ax0)
            gridViewer.plotTomo(t1,'Repeat','Distance [m]','Elevation [m]',ax1)
        end
        
        if hcolorbar.Value==1
            cmin = str2double(hcmin.String);
            cmax = str2double(hcmax.String);
            caxis(ax0,[cmin cmax])
            caxis(ax1,[cmin cmax])
        end
        hb0=colorbar('peer',ax0);
        hb1=colorbar('peer',ax1);
        colormap(nf,hcmap.String{hcmap.Value})
        
        load([rep,file],'mogs')
        d = mogs(hlistBaseline.Value).data;
        if param.tomoAtt==0
            units = [d.cunits,'/',d.tunits];
        else
            units = ['Np/',d.cunits];
        end
        set(get(hb0,'Title'),'String',units,'FontSize',12)
        set(get(hb1,'Title'),'String',units,'FontSize',12)
    end

    function showDiff(varargin)
        if isempty(tomo)
            return
        end
        if param.tomoAtt==0
            t = 1./tomo.s - 1./tomo.s0;
        else
            t = tomo.s - tomo.s0;
        end
        nf=figure;
        ax=axes('Parent',nf);
        if strcmp(model.grid.type,'3D')==1
            gv = GridViewer(model.grid);
            gv.createSliders('Parent',nf);
            gv.slider2.Visible = 'off';
            gv.plotTomo(t,'Difference','Distance [m]','Elevation [m]',ax)
        else
            gridViewer.plotTomo(t,'Difference','Distance [m]','Elevation [m]',ax)
        end
        hb=colorbar('peer',ax);
        colormap(nf,hcmap.String{hcmap.Value})
        load([rep,file],'mogs')
        d = mogs(hlistBaseline.Value).data;
        if param.tomoAtt==0
            units = [d.cunits,'/',d.tunits];
        else
            units = ['Np/',d.cunits];
        end
        set(get(hb,'Title'),'String',units,'FontSize',12)
    end

    function showDep(varargin)
        if isempty(tomo)
            return
        end
        s_ref = model.inv_res(param.ref_inv_no).tomo.s;
        if param.tomoAtt==0
            t1 = 1./tomo.s - 1./s_ref;
            t0 = 1./tomo.s0 - 1./s_ref;
        else
            t1 = tomo.s - s_ref;
            t0 = tomo.s0 - s_ref;
        end
        
        cmin = min([t1(:); t0(:)]);
        cmax = max([t1(:); t0(:)]);
        
        
        nf=figure;
        ax0=subplot(1,2,1,'Parent',nf);
        ax1=subplot(1,2,2,'Parent',nf);
        
        if strcmp(model.grid.type,'3D')==1
            gv = GridViewer(model.grid);
            gv.createSliders('Parent',nf);
            gv.plotTomo(t0,'Baseline','Distance [m]','Elevation [m]',ax0)
            gv.plotTomo(t1,'Repeat','Distance [m]','Elevation [m]',ax1)
        else
            gridViewer.plotTomo(t0,'Baseline','Distance [m]','Elevation [m]',ax0)
            gridViewer.plotTomo(t1,'Repeat','Distance [m]','Elevation [m]',ax1)
        end
        caxis(ax0,[cmin cmax])
        caxis(ax1,[cmin cmax])
        hb0=colorbar('peer',ax0);
        hb1=colorbar('peer',ax1);
        colormap(nf,hcmap.String{hcmap.Value})
        
        load([rep,file],'mogs')
        d = mogs(hlistBaseline.Value).data;
        if param.tomoAtt==0
            units = [d.cunits,'/',d.tunits];
        else
            units = ['Np/',d.cunits];
        end
        set(get(hb0,'Title'),'String',units,'FontSize',12)
        set(get(hb1,'Title'),'String',units,'FontSize',12)
    end

    function showRays(varargin)
        if isempty(tomo)
            return
        end
        if isempty(tomo.rays)
            return
        end
        nf=figure;
        ax=axes('Parent',nf);
        
        rmin = 1.001*min(tomo.invData(end).res);
        rmax = 1.001*max(tomo.invData(end).res);
        rmax = max(abs([rmin rmax]));
        rmin = -rmax;
        c = [0 0 1;0.8 0.8 0.8;1 0 0];
        c = interp1((-1:1)',c,(-1:0.02:1)');
        
        m = (size(c,1)-1)/(rmax-rmin);
        b = 1-rmin*m;
        p = m*tomo.invData(end).res(1)+b;
        
        if strcmp(model.grid.type,'3D')
            couleur = interp1(c,p);
            plot3(ax,tomo.rays{1}(:,1),tomo.rays{1}(:,2),tomo.rays{1}(:,3),...
                'Color',couleur)
            hold(ax,'on')
            for n=2:length(tomo.rays)
                p = m*tomo.invData(end).res(n)+b;
                couleur = interp1(c,p);
                plot3(ax,tomo.rays{n}(:,1),tomo.rays{n}(:,2),tomo.rays{n}(:,3),...
                    'Color',couleur)
            end
            hold(ax,'off')
            xlabel(ax,'X','FontSize',12)
            ylabel(ax,'Y','FontSize',12)
            zlabel('Elevation [m]','FontSize',12)
        else
            couleur = interp1(c,p);
            plot(ax,tomo.rays{1}(:,1),tomo.rays{1}(:,end),'Color',couleur)
            hold(ax,'on')
            for n=2:length(tomo.rays)
                p = m*tomo.invData(end).res(n)+b;
                if p>size(c,1)
                    p=size(c,1);
                elseif p<1
                    p=1;
                end
                    
                couleur = interp1(c,p);
                plot(ax,tomo.rays{n}(:,1),tomo.rays{n}(:,end),'Color',couleur)
            end
            hold(ax,'off')
            xlabel('Distance [m]','FontSize',12)
            ylabel('Elevation [m]','FontSize',12)
        end
        set(ax,'DataAspectRatio',[1 1 1])
        axis(ax,'tight')
        colormap(c)%jet)
        hb=colorbar('peer',ax);
        caxis(ax,[rmin rmax])
        set(get(hb,'Title'),'String','Residuals','FontSize',12)
    end

    function showRayDensity(varargin)
        if isempty(tomo)
            return
        end
        nf=figure;
        ax0=subplot(1,2,1,'Parent',nf);
        ax1=subplot(1,2,2,'Parent',nf);

        if isfield(tomo,'xi')
            nCells=size(tomo.L0,2)/2;
            Lx = tomo.L0(:,1:nCells);
            Lz = tomo.L0(:,(1+nCells):end);
            rd0 = full(sum(sqrt(Lx.^2+Lz.^2)));
            nCells=size(tomo.L,2)/2;
            Lx = tomo.L(:,1:nCells);
            Lz = tomo.L(:,(1+nCells):end);
            rd1 = full(sum(sqrt(Lx.^2+Lz.^2)));
        else
            rd0 = full(sum(tomo.L0));
            rd1 = full(sum(tomo.L));
        end
        if strcmp(model.grid.type,'3D')==1
            gv = GridViewer(model.grid);
            gv.createSliders('Parent',nf);
            gv.plotTomo(rd0,'Repeat','Distance [m]','Elevation [m]',ax0)
            gv.plotTomo(rd1,'Repeat','Distance [m]','Elevation [m]',ax1)
        else
            gridViewer.plotTomo(rd0,'Repeat','Distance [m]','Elevation [m]',ax0)
            gridViewer.plotTomo(rd1,'Repeat','Distance [m]','Elevation [m]',ax1)
        end
        hb0=colorbar('peer',ax0);
        hb1=colorbar('peer',ax1);
        colormap(nf,hcmap.String{hcmap.Value})
        set(get(hb0,'Title'),'String','Ray Density','FontSize',12)
        set(get(hb1,'Title'),'String','Ray Density','FontSize',12)
    end

    function showResiduals(varargin)
        if isempty(tomo)
            return
        end
        
        
        if param.tomoAtt == 0
            [data,idata] = Model.getModelData(model,[rep,file],'tt',param.mog_no0);
            [depth,~] = Model.getModelData(model,[rep,file],'depth',param.mog_no0,[],'tt');
            data = [model.grid.Tx(idata,:) model.grid.Rx(idata,:) data];
        else
            switch htypeData.Value
                case 2
                    [data,idata] = Model.getModelData(model,[rep,file],'amp',param.mog_no0);
                    [depth,~] = Model.getModelData(model,[rep,file],'depth',param.mog_no0,[],'tt');
                case 3
                    [data,idata] = Model.getModelData(model,[rep,file],'fce',param.mog_no0);
                    [depth,~] = Model.getModelData(model,[rep,file],'depth',param.mog_no0,[],'fce');
                case 4
                    [data,idata] = Model.getModelData(model,[rep,file],'hyb',param.mog_no0);
                    [depth,~] = Model.getModelData(model,[rep,file],'depth',param.mog_no0,[],'hyb');
            end
            data = [model.grid.Tx(idata,:) model.grid.Rx(idata,:) data];
        end
        hyp = sqrt( sum((data(:,1:3)-data(:,4:6)).^2, 2) );
        dz = data(:,6)-data(:,3);
        theta = 180/pi*asin(dz./hyp);

        nIt = length(tomo.invData);
        rms = zeros(nIt,1);
        for n=1:nIt
            rms(n) = rmsv(tomo.invData(n).res0);
        end
        
        
        figure
        subplot(2,2,1)
        plot(1:nIt, rms,'o')
        ylabel('||res||')
        xlabel('Iteration')
        
        res = tomo.invData(nIt).res0;
        subplot(2,2,2)
        plot(theta, res,'o')
        xlabel('Angle w/r to horizontal [deg]')
        ylabel('Residuals')
        
        vres = var(res);
        h1=subplot(2,2,3);
        hist(res,30)
        xlabel('Residuals')
        ylabel('Count')
        title(['\sigma^2 = ', num2str(vres)])
        
        dTx = sort(unique(depth(:,1)));
        dRx = sort(unique(depth(:,2)));
        imdata = nan(length(dTx),length(dRx));
        for i=1:length(dTx)
            for j=1:length(dRx)
                ind = dTx(i)==depth(:,1) & dRx(j)==depth(:,2);
                if sum(ind)==1
                    imdata(i,j) = res(ind);
                end
            end
        end
        
        p = [0 0 1;1 1 1;1 0 0];
        p = interp1((-1:1)',p,(-1:0.02:1)');
        
        z=imdata;
        z(isnan(imdata))=0;
        z(~isnan(imdata))=1;
        
        h2=subplot(2,2,4);
        imagesc(dRx,dTx,imdata);
        set(gca,'color',[0.8 0.8 0.8]);
        alpha(z);
        axis image;
        
        ca = caxis;
        caxis([-max(abs(ca)) max(abs(ca))])
        
        xlabel('Rx depth')
        ylabel('Tx depth')
        colorbar
        
        colormap(h2,p)
        set(get(h1,'Children'),'FaceColor',[0 0 1])
        suptitle('Baseline')
        
        
        if param.tomoAtt == 0
            [data,idata] = Model.getModelData(model,[rep,file],'tt',param.mog_no1);
            [depth,~] = Model.getModelData(model,[rep,file],'depth',param.mog_no1,[],'tt');
            data = [model.grid.Tx(idata,:) model.grid.Rx(idata,:) data];
        else
            switch htypeData.Value
                case 2
                    [data,idata] = Model.getModelData(model,[rep,file],'amp',param.mog_no1);
                    [depth,~] = Model.getModelData(model,[rep,file],'depth',param.mog_no1,[],'tt');
                case 3
                    [data,idata] = Model.getModelData(model,[rep,file],'fce',param.mog_no1);
                    [depth,~] = Model.getModelData(model,[rep,file],'depth',param.mog_no1,[],'fce');
                case 4
                    [data,idata] = Model.getModelData(model,[rep,file],'hyb',param.mog_no1);
                    [depth,~] = Model.getModelData(model,[rep,file],'depth',param.mog_no1,[],'hyb');
            end
            data = [model.grid.Tx(idata,:) model.grid.Rx(idata,:) data];
        end
        hyp = sqrt( sum((data(:,1:3)-data(:,4:6)).^2, 2) );
        dz = data(:,6)-data(:,3);
        theta = 180/pi*asin(dz./hyp);

        nIt = length(tomo.invData);
        rms = zeros(nIt,1);
        for n=1:nIt
            rms(n) = rmsv(tomo.invData(n).res);
        end
        
        
        figure
        subplot(2,2,1)
        plot(1:nIt, rms,'o')
        ylabel('||res||')
        xlabel('Iteration')
        
        res = tomo.invData(nIt).res;
        subplot(2,2,2)
        plot(theta, res,'o')
        xlabel('Angle w/r to horizontal [deg]')
        ylabel('Residuals')
        
        vres = var(res);
        h1=subplot(2,2,3);
        hist(res,30)
        xlabel('Residuals')
        ylabel('Count')
        title(['\sigma^2 = ', num2str(vres)])
        
        dTx = sort(unique(depth(:,1)));
        dRx = sort(unique(depth(:,2)));
        imdata = nan(length(dTx),length(dRx));
        for i=1:length(dTx)
            for j=1:length(dRx)
                ind = dTx(i)==depth(:,1) & dRx(j)==depth(:,2);
                if sum(ind)==1
                    imdata(i,j) = res(ind);
                end
            end
        end
        
        p = [0 0 1;1 1 1;1 0 0];
        p = interp1((-1:1)',p,(-1:0.02:1)');
        
        z=imdata;
        z(isnan(imdata))=0;
        z(~isnan(imdata))=1;
        
        h2=subplot(2,2,4);
        imagesc(dRx,dTx,imdata);
        set(gca,'color',[0.8 0.8 0.8]);
        alpha(z);
        axis image;
        
        ca = caxis;
        caxis([-max(abs(ca)) max(abs(ca))])
        
        xlabel('Rx depth')
        ylabel('Tx depth')
        colorbar
        
        colormap(h2,p)
        set(get(h1,'Children'),'FaceColor',[0 0 1])
        suptitle('Repeat')
        
    end

    function loadPrevious(varargin)
        if isempty(model)
            return
        end
        if isempty(model.tlinv_res)
            return
        end
        no = hpreviousInv.Value;
        tomo = model.tlinv_res(no).tomo;
        param = model.tlinv_res(no).param;
    end
    function deletePrevious(varargin)
        if isempty(model)
            return
        end
        if isempty(model.tlinv_res)
            return
        end
        no = hpreviousInv.Value;
        nos=1:numel(model.tlinv_res);
        ind = nos~=no;
        model.tlinv_res = model.tlinv_res(ind);
        
        if isempty(model.tlinv_res)
            hpreviousInv.String = {'-'};
        else
            names = cell(numel(model.tlinv_res),1);
            for n=1:numel(model.tlinv_res)
                names{n} = model.tlinv_res(n).name;
            end
            hpreviousInv.String = names;
        end
        tomo = [];
        saved = false;
    end

    function showDelta(varargin)
        
        if htypeData.Value == 1
            tomoAtt = 0;
            titre = '\Delta t (t_1 - t_0)';
        else
            tomoAtt = 1;
            titre = '\Delta \tau (\tau_1 - \tau_0)';
        end
        
        mog_no0 = hlistBaseline.Value;
        mog_no1 = hlistRepeat.Value;
            
        if tomoAtt == 0
            [data0,idata] = Model.getModelData(model,[rep,file],'tt',mog_no0);
            [depth0,~] = Model.getModelData(model,[rep,file],'depth',mog_no0,[],'tt');
            data0 = [model.grid.Tx(idata,:) model.grid.Rx(idata,:) data0];
        else
            switch htypeData.Value
                case 2
                    [data0,idata] = Model.getModelData(model,[rep,file],'amp',mog_no0);
                    [depth0,~] = Model.getModelData(model,[rep,file],'depth',mog_no0,[],'tt');
                case 3
                    [data0,idata] = Model.getModelData(model,[rep,file],'fce',mog_no0);
                    [depth0,~] = Model.getModelData(model,[rep,file],'depth',mog_no0,[],'fce');
                case 4
                    [data0,idata] = Model.getModelData(model,[rep,file],'hyb',mog_no0);
                    [depth0,~] = Model.getModelData(model,[rep,file],'depth',mog_no0,[],'hyb');
            end
            data0 = [model.grid.Tx(idata,:) model.grid.Rx(idata,:) data0];
        end

        if tomoAtt == 0
            [data1,idata] = Model.getModelData(model,[rep,file],'tt',mog_no1);
%            [depth1,~] = Model.getModelData(model,[rep,file],'depth',mog_no1,[],'tt');
            data1 = [model.grid.Tx(idata,:) model.grid.Rx(idata,:) data1];
        else
            switch htypeData.Value
                case 2
                    [data1,idata] = Model.getModelData(model,[rep,file],'amp',mog_no1);
%                    [depth1,~] = Model.getModelData(model,[rep,file],'depth',mog_no1,[],'tt');
                case 3
                    [data1,idata] = Model.getModelData(model,[rep,file],'fce',mog_no1);
%                    [depth1,~] = Model.getModelData(model,[rep,file],'depth',mog_no1,[],'fce');
                case 4
                    [data1,idata] = Model.getModelData(model,[rep,file],'hyb',mog_no1);
%                    [depth1,~] = Model.getModelData(model,[rep,file],'depth',mog_no1,[],'hyb');
            end
            data1 = [model.grid.Tx(idata,:) model.grid.Rx(idata,:) data1];
        end

        dx = model.grid.grx(2) - model.grid.grx(1);
        dz = model.grid.grz(2) - model.grid.grz(1);
        tol = 0.5*(dx+dz);
        
        dTx = sort(unique(depth0(:,1)));
        dRx = sort(unique(depth0(:,2)));
        if length(dTx) > 2*length(model.grid.grz)
            o = floor(log10(dz)-1);
            fac = 10^o;
            depth0(:,1) = round(depth0(:,1)/fac) * fac;
            dTx = sort(unique(depth0(:,1)));
        end
        if length(dRx) > 2*length(model.grid.grz)
            o = floor(log10(dz)-1);
            fac = 10^o;
            depth0(:,2) = round(depth0(:,2)/fac) * fac;
            dRx = sort(unique(depth0(:,2)));
        end
            
            
        ddata = nan(length(dTx),length(dRx));
        nval = zeros(length(dTx),length(dRx));
        
        for n0=1:size(data0,1)
            for n1=1:size(data1,1)
                
                if sqrt(sum((data0(n0,1:3)-data1(n1,1:3)).^2)) < tol && ...
                        sqrt(sum((data0(n0,4:6)-data1(n1,4:6)).^2)) < tol
                    iTx = (dTx == depth0(n0,1));
                    iRx = (dRx == depth0(n0,2));
                    if isnan(ddata(iTx,iRx))
                        ddata(iTx,iRx) = data1(n1,7) - data0(n0,7);
                    else
                        ddata(iTx,iRx) = ddata(iTx,iRx) + data1(n1,7) - data0(n0,7);
                    end
                    nval(iTx,iRx) = nval(iTx,iRx) + 1;
                    
                end
            end
        end
        
        ind = nval>0;
        ddata(ind) = ddata(ind) ./ nval(ind);
        
        figure
        
        p = [0 0 1;1 1 1;1 0 0];
        p = interp1((-1:1)',p,(-1:0.02:1)');
        
        z=ddata;
        z(isnan(ddata))=0;
        z(~isnan(ddata))=1;
        
        imagesc(dRx,dTx,ddata);
        set(gca,'color',[0.8 0.8 0.8]);
        alpha(z);
        axis image;
        
        ca = caxis;
        caxis([-max(abs(ca)) max(abs(ca))])
        
        xlabel('Rx depth')
        ylabel('Tx depth')
        colorbar
        colormap(p)
        
        title(titre, 'FontSize',16)
        
    end

end

