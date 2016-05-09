function bh_tomo2_interp( varargin )
%BH_TOMO2_INTERP Compute physical properties from tomograms

rep='';
file='';
if nargin>=2
    rep=varargin{1};
    file=varargin{2};
end

% Main variables

saved = true;
model = Model.empty();
islo = [];
iatt = [];
gv = GridViewer.empty();

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

width = 1300;
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
    'Tag','fig_bh_tomo2_interp',...
    'Name','bh_tomo_interp',...
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
uimenu(hmenu,'Label','Save Figure',...
    'Accelerator','S',...
    'Callback',@saveFigure);
uimenu(hmenu,'Label','Close',...
    'Separator','on',...
    'Accelerator','W',...
    'Callback',@closeWindow);




%
% Main Panels
%
ptomo = uipanel(f,'Title','Tomograms',...
    'Units','points',...
    'FontSize',fs+1);
pprop = uipanel(f,'Title','Physical Property',...
    'Units','points',...
    'FontSize',fs+1);
ppetro = uipanel(f,'Title','Petrophysical Model',...
    'Units','points',...
    'FontSize',fs+1);
pfig = uipanel(f,'Title','Figures',...
    'Units','points',...
    'FontSize',fs+1);


nLines=1;
vSizeTot = nLines*22 + 2*5;
vSize = 22/vSizeTot;
vSpace = 5/vSizeTot;

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

m = {'cmr','polarmap','parula','jet','hsv','hot','cool','autumn','spring','winter',...
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

nLines=6;
vSizeTot = nLines*22 + (nLines+1)*5;
vSize = 22/vSizeTot;
vSpace = 5/vSizeTot;

uicontrol('Style','text',...
    'String','Slowness',...
    'FontSize',fs,...
    'HorizontalAlignment','left',...
    'Units','normalized',...
    'Position',[0.1 5*vSize+6*vSpace 0.8 vSize],...
    'Parent',ptomo);
hslo = uicontrol('Style','popupmenu',...
    'String',{'-'},...
    'FontSize',fs,...
    'Units','normalized',...
    'Position',[0.1 4*vSize+5*vSpace 0.8 vSize],...
    'Parent',ptomo);
uicontrol('Style','text',...
    'String','Attenuation',...
    'FontSize',fs,...
    'HorizontalAlignment','left',...
    'Units','normalized',...
    'Position',[0.1 3*vSize+4*vSpace 0.8 vSize],...
    'Parent',ptomo);
hatt = uicontrol('Style','popupmenu',...
    'String',{'-'},...
    'FontSize',fs,...
    'Units','normalized',...
    'Position',[0.1 2*vSize+3*vSpace 0.8 vSize],...
    'Parent',ptomo);
uicontrol('Style','text',...
    'String','Type',...
    'FontSize',fs,...
    'HorizontalAlignment','left',...
    'Units','normalized',...
    'Position',[0.1 vSize+2*vSpace 0.8 vSize],...
    'Parent',ptomo);
htype = uicontrol('Style','popupmenu',...
    'String',{'Cokriging','Cosimulation'},...
    'FontSize',fs,...
    'Units','normalized',...
    'Position',[0.1 vSpace 0.8 vSize],...
    'Callback',@plotFig,...
    'Parent',ptomo);

nLines=2;
vSizeTot = nLines*22 + (nLines+1)*5;
vSize = 22/vSizeTot;
vSpace = 5/vSizeTot;

pp={'Velocity',...
    'Slowness',...
    ['Attenuation ',char(945)],...
    ['Volumetric Water Content ',char(952)]...
    ['Dielectric Permittivity ',char(949)],...
    ['Dielectric Constant ',char(954)],...
    ['Electric Conductivity ',char(963)],...
    ['Electric Resistivity ',char(961)],...
    ['Loss Tangent ',char(948)]};

hprop = uicontrol('Style','popupmenu',...
    'String',pp,...
    'FontSize',fs,...
    'Units','normalized',...
    'Position',[0.1 vSize+2*vSpace 0.8 vSize],...
    'Callback',@plotFig,...
    'Parent',pprop);
uicontrol('Style','text',...
    'String','Frequency [MHz]',...
    'FontSize',fs,...
    'HorizontalAlignment','right',...
    'Units','normalized',...
    'Position',[0.1 vSpace 0.38 vSize],...
    'Parent',pprop);
hfreq = uicontrol('Style','edit',...
    'String','100',...
    'FontSize',fs,...
    'HorizontalAlignment','center',...
    'Units','normalized',...
    'Position',[0.5 vSpace 0.2 vSize],...
    'Parent',pprop);

hpetro = uicontrol('Style','popupmenu',...
    'String',{'Topp','CRIM','Hanai-Bruggeman'},...
    'FontSize',fs,...
    'Units','normalized',...
    'Position',[0.1 vSize+2*vSpace 0.8 vSize],...
    'Callback',@modelPetro,...
    'Parent',ppetro);
htwater = uicontrol('Style','text',...
    'String',[char(954),' Water'],...
    'FontSize',fs,...
    'HorizontalAlignment','right',...
    'Units','normalized',...
    'Position',[0.05 vSpace 0.18 vSize],...
    'Visible','off',...
    'Parent',ppetro);
hkwater = uicontrol('Style','edit',...
    'String','80',...
    'FontSize',fs,...
    'Units','normalized',...
    'Position',[0.25 vSpace 0.18 vSize],...
    'Visible','off',...
    'Callback',@plotFig,...
    'Parent',ppetro);
htmatrix = uicontrol('Style','text',...
    'String',[char(954),' Matrix'],...
    'FontSize',fs,...
    'HorizontalAlignment','right',...
    'Units','normalized',...
    'Position',[0.5 vSpace 0.18 vSize],...
    'Visible','off',...
    'Parent',ppetro);
hkmatrix = uicontrol('Style','edit',...
    'String','5',...
    'FontSize',fs,...
    'Units','normalized',...
    'Position',[0.7 vSpace 0.18 vSize],...
    'Visible','off',...
    'Callback',@plotFig,...
    'Parent',ppetro);

haxes1 = axes('Units','points','Parent',f,'Visible','on');


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
        hSize = 325;

        vSize2 = 6.5*vSize+7*vSpace;
        vPos = height-2*vBorder-vSize2;
        ptomo.Position = [hBorder vPos hSize vSize2];
        
        vSize2 = 2.5*vSize+3*vSpace;
        vPos = vPos-vBorder-vSize2;
        pprop.Position = [hBorder vPos hSize vSize2];
        
        vSize2 = 2.5*vSize+3*vSpace;
        vPos = vPos-vBorder-vSize2;
        ppetro.Position = [hBorder vPos hSize vSize2];
        
        vSize2 = 1.5*vSize+2*vSpace;
        hSize2 = width-(hBorder+hSize+axBorder+rBorder);
        vPos = height-2*vBorder-vSize2;
        pfig.Position = [hBorder+hSize+axBorder vPos hSize2 vSize2];
        
        
        vSize2 = height-2*axBorder-2*vBorder-vSize2;
        haxes1.Position = [hBorder+hSize+axBorder axBorder hSize2 vSize2];
        
        
        f.Visible = 'on';
    end

    function closeWindow(varargin)
        if saved == false
            choice = questdlg('Database not saved, quit anyway?',...
                'bh_tomo_db',...
                'Don''t save','Cancel','Save Figure','Save Figure');
            switch choice
                case 'Don''t save'
                case 'Cancel'
                    return
                case 'Save'
                    saveFigure()
            end
        end
        quitUI()
    end
    function quitUI(varargin)
        delete(f)
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
        if isempty(model.inv_res)
            errordlg('No Inversion Results For This Model')
            return
        end
        
        %
        % Reset UI
        %
        cla(haxes1);
        hsliderText.Visible = 'off';
        cbh = findobj( 0, 'tag', 'Colorbar' );
        for i = 1:length(cbh)
            colorbar(cbh(i),'off')
        end

        gv = GridViewer(model.grid);
        if strcmp(model.grid.type, '3D')
            hsliderText.Visible = 'on';
            gv.createSliders('Parent',pfig,...
                'Units','normalized','Visible','off');
            nLines=1;
            vSizeTot = nLines*22 + 2*5;
            vSize = 22/vSizeTot;
            vSpace = 5/vSizeTot;
            gv.slider1.Position = [0.64 vSpace 0.15 vSize];
            gv.slider1.Visible = 'on';
            gv.slider2.Position = [0.80 vSpace 0.15 vSize];
        end
        
        nslo = cell(length(model.inv_res),1);
        natt = cell(length(model.inv_res),1);
        islo = zeros(length(model.inv_res),1);
        iatt = zeros(length(model.inv_res),1);
        ns = 0;
        na = 0;
        for n=1:length(model.inv_res)
            if model.inv_res(n).param.tomoAtt==0
                ns = ns+1;
                nslo{ns} = model.inv_res(n).name;
                islo(ns) = n;
            else
                na = na+1;
                natt{na} = model.inv_res(n).name;
                iatt(na) = n;
            end
        end
        if ns>0
            nslo = nslo{1:ns};
            islo = islo(1:ns);
            hslo.String = nslo;
        else
            islo = [];
        end
        if na>0
            natt = natt{1:na};
            iatt = iatt(1:na);
            hatt.String = natt;
        else
            iatt = [];
        end
        
        plotFig()
    end
    function saveFigure(varargin)
        [filename,pathname,filterindex] = uiputfile({'Portable Document Format *.pdf';...
            'Encapsulated Postscript *.eps';...
            'Portable Network Graphics file *.png';...
            'TIFF image *.tif';...
            'TIFF no compression image *.tif';...
            'Scalable Vector Graphics file *.svg';...
            'Bitmap file *.bmp'});
        if isequal(filename,0) || isequal(pathname,0)
            return
        end
        f.PaperPositionMode = 'auto';
        ptomo.Visible = 'off';
        pprop.Visible = 'off';
        ppetro.Visible = 'off';
        pfig.Visible = 'off';
        
        switch filterindex
            case 1
                print(f,[pathname,filename],'-dpdf','-noui')
            case 2
                print(f,[pathname,filename],'-depsc2','-noui')
            case 3
                print(f,[pathname,filename],'-dpng','-noui')
            case 4
                print(f,[pathname,filename],'-dtiff','-noui')
            case 5
                print(f,[pathname,filename],'-dtiffnocompression','-noui')
            case 6
                print(f,[pathname,filename],'-dsvg','-noui')
            case 7
                print(f,[pathname,filename],'-dbmp16m','-noui')
        end
        ptomo.Visible = 'on';
        pprop.Visible = 'on';
        ppetro.Visible = 'on';
        pfig.Visible = 'on';

    end
    function plotFig(varargin)
        cla(haxes1)
        cbh = findobj( 0, 'tag', 'Colorbar' );
        for i = 1:length(cbh)
            colorbar(cbh(i),'off')
        end
        switch hprop.Value
            case 1
                try
                    data = 1./getSlowness();
                catch ME
                    warndlg(ME.msgtext)
                    return
                end
            case 2
                try
                    data = getSlowness();
                catch ME
                    warndlg(ME.message)
                    return
                end
            case 3
                try
                    data = getAttenuation();
                catch ME
                    warndlg(ME.message)
                    return
                end
            case 4
                try
                    data = getWaterContent();
                catch ME
                    warndlg(ME.message)
                    return
                end
            case 5
                try
                    data = getK()*8.8541878e-12;
                catch ME
                    warndlg(ME.message)
                    return
                end
            case 6
                try
                    data = getK();
                catch ME
                    warndlg(ME.message)
                    return
                end
            case 7
                try
                    data = getSigma();
                catch ME
                    warndlg(ME.message)
                    return
                end
            case 8
                try
                    data = 1./getSigma();
                catch ME
                    warndlg(ME.message)
                    return
                end
            case 9
                try
                    sig = getSigma();
                    k = getK();
                catch ME
                    warndlg(ME.message)
                    return
                end
                f = str2double(hfreq.String)*1e6;
                w = 2*pi*f;
                eps0 = 8.8541878e-12;
                e=k*eps0;
                data = sig./(w*e);
        end
        gv.plotTomo(data,hprop.String{hprop.Value},'Distance [m]','Elevation [m]',haxes1)
        colorbar('peer',haxes1)
        colormap(hcmap.String{hcmap.Value})

    end
    function doColorbar(varargin)
        if hcolorbar.Value==1
            cl = [str2double(hcmin.String) str2double(hcmax.String)];
            haxes1.CLim = cl;
        else
            haxes1.CLimMode = 'auto';
        end
    end
    function doMap(varargin)
        colormap(f, hcmap.String{hcmap.Value})
    end
    function setClim(varargin)
        % TODO: store values
        
        if hcolorbar.Value==1
            doColorbar()
        end
    end
    function modelPetro(varargin)
        if hpetro.Value==1
            htwater.Visible = 'off';
            hkwater.Visible = 'off';
            htmatrix.Visible = 'off';
            hkmatrix.Visible = 'off';
        else
            htwater.Visible = 'on';
            hkwater.Visible = 'on';
            htmatrix.Visible = 'on';
            hkmatrix.Visible = 'on';
        end
        if hprop.Value == 4
            plotFig()
        end
    end

    function data = getSlowness()
        if isempty(islo)
            ME = MException('bh_tomo_interp:noSuchVariable',...
                'No Slowness Results For This Model');
            throw(ME)
        end
        if htype.Value==2
            if ~isfield(model.inv_res(islo(hslo.Value)).tomo,'sgr')
                ME = MException('bh_tomo_interp:noSuchVariable',...
                    'Simulations Were Not Performed');
                throw(ME)
            end
            data = model.inv_res(islo(hslo.Value)).tomo.sgr;
        else
            data = model.inv_res(islo(hslo.Value)).tomo.s;
        end
    end
    function data = getAttenuation()
        if isempty(iatt)
            ME = MException('bh_tomo_interp:noSuchVariable',...
                'No Attenuation Results For This Model');
            throw(ME)
        end
        if htype.Value==2
            if ~isfield(model.inv_res(iatt(hatt.Value)).tomo,'sgr')
                ME = MException('bh_tomo_interp:noSuchVariable',...
                    'Simulations Were Not Performed');
                throw(ME)
            end
            data = model.inv_res(iatt(hatt.Value)).tomo.sgr;
        else
            data = model.inv_res(iatt(hatt.Value)).tomo.s;
        end
    end
    function k = getK()
        try
            s = getSlowness();
        catch
            rethrow(ME)
        end
        try
            a = getAttenuation();
        catch
            a = [];
        end
        if isempty(a)
            k = 0.2298^2 * s.^2;
        else
            f = str2double(hfreq.String)*1e6;
            s = s*1e-9;
            mu0 = 4*pi*1e-7;
            eps0 = 8.8541878e-12;
            w = 2*pi*f;
            k = (1/(mu0*eps0))*(s.^2 - (a./w).^2);
        end
    end
    function s = getSigma()
        try
            s = getSlowness();
        catch ME
            rethrow(ME)
        end
        try
            a = getAttenuation();
        catch ME
            rethrow(ME)
        end
        s = s*1e-9;
        mu0 = 4*pi*1e-7;
        s = 2*a.*s./mu0;
    end
    function theta = getWaterContent()
        try
            k = getK();
        catch ME
            rethrow(ME)
        end
        switch hpetro.Value
            case 1
                theta = topp(k);
            case 2
                k_w = str2double(hkwater.String);
                k_m = str2double(hkmatrix.String);
                theta = (sqrt(k) - sqrt(k_m))./(sqrt(k_w) - sqrt(k_m));
            case 3
                k_w = str2double(hkwater.String);
                k_m = str2double(hkmatrix.String);
                W = 1/3;
                theta = 1-((k_w-k)./(k_w-k_m)).*(k_m./k).^W;
        end
    end
end

