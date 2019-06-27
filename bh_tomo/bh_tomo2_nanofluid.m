function bh_tomo2_nanofluid(varargin)


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

saved = true;
model = Model.empty();
modelNo = [];
gridViewer = [];

nf_sat.epsilon_m = 8;
nf_sat.a = 0.62;
nf_sat.m = 1.72;
nf_sat.sigma_nf = 0.01;
nf_sat.epsilon_nf = 75;
nf_sat.mu_nf = 3;
nf_sat.epsilon_f0 = 2.25;
nf_sat.sigma_f0 = 0.001;
nf_sat.omega = 2*pi*1e6*50;
nf_sat.S0 = 0.5;
sat_map = [];

width = 1100;
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
    'Tag','fig_bh_tomo2_nanofluid',...
    'Name','bh_tomo_nanofluid',...
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
uimenu(hresultsMenu,'Label','Plot',...
    'Accelerator','P',...
    'Callback',@plot_res);



%
% Main Panels
%
pinput = uipanel(f,'Title','Input Data',...
    'Units','points',...
    'FontSize',fs+2);

pfig = uipanel(f,'Title','Figures',...
    'Units','points',...
    'FontSize',fs+1);

psat_map = uipanel(f,'Title','Saturation Maps',...
    'Units','points',...
    'FontSize',fs+1);
hsat_map = uicontrol('Style','popupmenu',...
    'String','-',...
    'FontSize',fs,...
    'Units','normalized',...
    'Position',[0.04 0.1 0.92 0.8],...
    'Parent',psat_map);

bcompute = uicontrol('Style','pushbutton',...
    'String','Compute',...
    'FontSize',fs+1,...
    'FontWeight','bold',...
    'Units','points',...
    'Callback',@compute,...
    'Parent',f);


nLines=1;
vSizeTot = nLines*22 + 2*5;
vSize = 22/vSizeTot;
vSpace = 5/vSizeTot;

hcolorbar = uicontrol('Style','checkbox',...
    'String','Set Color Limits',...
    'FontSize',fs,...
    'Units','normalized',...
    'Position',[0.02 vSpace 0.19 vSize],...
    'Callback',@doColorbar,...
    'Parent',pfig);
uicontrol('Style','text',...
    'String','Min: ',...
    'HorizontalAlignment','right',...
    'FontSize',fs,...
    'Units','normalized',...
    'Position',[0.21 vSpace 0.06 vSize],...
    'Parent',pfig)
hcmin = uicontrol('Style','edit',...
    'String','0.0',...
    'FontSize',fs,...
    'Units','normalized',...
    'Position',[0.28 vSpace 0.07 vSize],...
    'Callback',@setClim,...
    'Parent',pfig);
uicontrol('Style','text',...
    'String','Max: ',...
    'HorizontalAlignment','right',...
    'FontSize',fs,...
    'Units','normalized',...
    'Position',[0.36 vSpace 0.06 vSize],...
    'Parent',pfig)
hcmax = uicontrol('Style','edit',...
    'String','1',...
    'FontSize',fs,...
    'Units','normalized',...
    'Position',[0.43 vSpace 0.07 vSize],...
    'Callback',@setClim,...
    'Parent',pfig);

m = {'plasma','magma','inferno','viridis','cmr','polarmap','parula','jet',...
    'hsv','hot','cool','autumn','spring','winter',...
    'summer','gray','bone','copper','pink','prism','flag','colorcube','lines'};

hcmap = uicontrol('Style','popupmenu',...
    'String',m,...
    'FontSize',fs,...
    'Units','normalized',...
    'Position',[0.52 vSpace 0.2 vSize],...
    'Callback',@doMap,...
    'Parent',pfig);
hsliderText = uicontrol('Style','text',...
    'String','Y Plane',...
    'FontSize',fs,...
    'Units','normalized',...
    'Position',[0.75 vSpace 0.24 vSize],...
    'Visible','off',...
    'Parent',pfig);


vPos = 0.02;
vSize = 0.94 * 3/5 * 3/7;
pfluid = uipanel(pinput,'Title','Properties - Fluid in Place',...
    'Units','normalized',...
    'Position',[0.02 vPos 0.47 vSize],...
    'FontSize',fs+1);

uicontrol('Style','text',...
    'String','Relative permittivity',...
    'FontSize',fs,...
    'Units','normalized',...
    'HorizontalAlignment','right',...
    'Position',[0.02 0.5 0.65 0.35],...
    'Parent',pfluid);
hepsilon_f0 = uicontrol('Style','edit',...
    'String',num2str(nf_sat.epsilon_f0),...
    'FontSize',fs,...
    'Units','normalized',...
    'Position',[0.7 0.53 0.28 0.37],...
    'Parent',pfluid);
uicontrol('Style','text',...
    'String','Conductivity (S/m)',...
    'FontSize',fs,...
    'Units','normalized',...
    'HorizontalAlignment','right',...
    'Position',[0.02 0.02 0.65 0.35],...
    'Parent',pfluid);
hsigma_f0 = uicontrol('Style','edit',...
    'String',num2str(nf_sat.sigma_f0),...
    'FontSize',fs,...
    'Units','normalized',...
    'Position',[0.7 0.05 0.28 0.37],...
    'Parent',pfluid);

pparams = uipanel(pinput,'Title','Parameters',...
    'Units','normalized',...
    'Position',[0.51 vPos 0.47 vSize],...
    'FontSize',fs+1);

uicontrol('Style','text',...
    'String','Frequency (MHz)',...
    'FontSize',fs,...
    'Units','normalized',...
    'HorizontalAlignment','right',...
    'Position',[0.02 0.5 0.65 0.35],...
    'Parent',pparams);
hfreq = uicontrol('Style','edit',...
    'String',num2str(nf_sat.omega/2/pi * 1e-6),...
    'FontSize',fs,...
    'Units','normalized',...
    'Position',[0.7 0.53 0.28 0.37],...
    'Parent',pparams);
uicontrol('Style','text',...
    'String','Starting value of S',...
    'FontSize',fs,...
    'Units','normalized',...
    'HorizontalAlignment','right',...
    'Position',[0.02 0.02 0.65 0.35],...
    'Parent',pparams);
hS0 = uicontrol('Style','edit',...
    'String',num2str(nf_sat.S0),...
    'FontSize',fs,...
    'Units','normalized',...
    'Position',[0.7 0.05 0.28 0.37],...
    'Parent',pparams);

vPos = vPos + 0.01 + vSize;
vSize = 0.94 * 3/5 * 4/7;
pmatrix = uipanel(pinput,'Title','Properties - Rock Matrix',...
    'Units','normalized',...
    'Position',[0.02 vPos 0.47 vSize],...
    'FontSize',fs+1);

uicontrol('Style','text',...
    'String','Relative permittivity',...
    'FontSize',fs,...
    'Units','normalized',...
    'HorizontalAlignment','right',...
    'Position',[0.02 0.66 0.65 0.24],...
    'Parent',pmatrix);
hepsilon_m = uicontrol('Style','edit',...
    'String',num2str(nf_sat.epsilon_m),...
    'FontSize',fs,...
    'Units','normalized',...
    'Position',[0.7 0.69 0.28 0.25],...
    'Parent',pmatrix);
uicontrol('Style','text',...
    'String','Archie parameter a',...
    'FontSize',fs,...
    'Units','normalized',...
    'HorizontalAlignment','right',...
    'Position',[0.02 0.34 0.65 0.24],...
    'Parent',pmatrix);
ha = uicontrol('Style','edit',...
    'String',num2str(nf_sat.a),...
    'FontSize',fs,...
    'Units','normalized',...
    'Position',[0.7 0.37 0.28 0.25],...
    'Parent',pmatrix);
uicontrol('Style','text',...
    'String','Archie parameter m',...
    'FontSize',fs,...
    'Units','normalized',...
    'HorizontalAlignment','right',...
    'Position',[0.02 0.02 0.65 0.24],...
    'Parent',pmatrix);
hm = uicontrol('Style','edit',...
    'String',num2str(nf_sat.m),...
    'FontSize',fs,...
    'Units','normalized',...
    'Position',[0.7 0.05 0.28 0.25],...
    'Parent',pmatrix);

pnano = uipanel(pinput,'Title','Properties - Nano Fluid',...
    'Units','normalized',...
    'Position',[0.51 vPos 0.47 vSize],...
    'FontSize',fs+1);

uicontrol('Style','text',...
    'String','Relative permittivity',...
    'FontSize',fs,...
    'Units','normalized',...
    'HorizontalAlignment','right',...
    'Position',[0.02 0.66 0.65 0.24],...
    'Parent',pnano);
hepsilon_nf = uicontrol('Style','edit',...
    'String',num2str(nf_sat.epsilon_nf),...
    'FontSize',fs,...
    'Units','normalized',...
    'Position',[0.7 0.69 0.28 0.25],...
    'Parent',pnano);
uicontrol('Style','text',...
    'String','Relative permeability',...
    'FontSize',fs,...
    'Units','normalized',...
    'HorizontalAlignment','right',...
    'Position',[0.02 0.34 0.65 0.24],...
    'Parent',pnano);
hmu_nf = uicontrol('Style','edit',...
    'String',num2str(nf_sat.mu_nf),...
    'FontSize',fs,...
    'Units','normalized',...
    'Position',[0.7 0.37 0.28 0.25],...
    'Parent',pnano);
uicontrol('Style','text',...
    'String','Conductivity (S/m)',...
    'FontSize',fs,...
    'Units','normalized',...
    'HorizontalAlignment','right',...
    'Position',[0.02 0.02 0.65 0.24],...
    'Parent',pnano);
hsigma_nf = uicontrol('Style','edit',...
    'String',num2str(nf_sat.sigma_nf),...
    'FontSize',fs,...
    'Units','normalized',...
    'Position',[0.7 0.05 0.28 0.25],...
    'Parent',pnano);

vPos = vPos + 0.01 + vSize;
vSize = 0.94 * 2/5;
ptlinv = uipanel(pinput,'Title','Time-lapse inversion data',...
    'Units','normalized',...
    'Position',[0.02 vPos 0.96 vSize],...
    'FontSize',fs+1);

hdata_att = uicontrol('Style','popupmenu',...
    'String','-',...
    'FontSize',fs,...
    'Units','normalized',...
    'Position',[0.04 0.02 0.92 0.24],...
    'Parent',ptlinv);
uicontrol('Style','text',...
    'String','Attenuation',...
    'FontSize',fs,...
    'Units','normalized',...
    'HorizontalAlignment','left',...
    'Position',[0.04 0.28 0.92 0.2],...
    'Parent',ptlinv);
hdata_vel = uicontrol('Style','popupmenu',...
    'String','-',...
    'FontSize',fs,...
    'Units','normalized',...
    'Position',[0.04 0.50 0.92 0.24],...
    'Parent',ptlinv);
uicontrol('Style','text',...
    'String','Velocity',...
    'FontSize',fs,...
    'Units','normalized',...
    'HorizontalAlignment','left',...
    'Position',[0.04 0.76 0.92 0.2],...
    'Parent',ptlinv);


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
        hSize = 400;

        vSize2 = 2.2*vSize;
        vPos = height-12*vBorder-vSize2;
        psat_map.Position = [hBorder vPos hSize vSize2];
        
        vSize2 = 425*vScale-4*vBorder-vSpace-vSize;
        vPos = vPos-4*vBorder-vSize2;
        pinput.Position = [hBorder vPos hSize vSize2];
        
        vPos = vPos - vSize - 4*vBorder;
        bcompute.Position = [2*hBorder vPos hSize-2*hBorder vSize];
        
        vSize2 = 2.2*vSize;
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
        if isempty(model.tlinv_res)
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

        names = cell(1,length(model.tlinv_res));
        for n=1:length(model.tlinv_res)
            names{n} = model.tlinv_res(n).name;
        end
        hdata_vel.String = names;
        hdata_att.String = names;
        
        if ~isempty(model.nf_sat)
            names = cell(1,length(model.nf_sat));
            for n=1:length(model.nf_sat)
                names{n} = model.nf_sat(n).name;
            end
            hsat_map.String = names;
        else
            hsat_map.String = '-';
        end
        
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

    function saveFile(varargin)
        if isempty(model)
            return
        end
        if saved
            warndlg('Results already saved')
            return
        end
        
        nameDefault = ['saturation ',nf_sat.vel_data, ' & ',nf_sat.att_data];  
        prompt = {'Saturation result name:                                                               '};
        title = 'Save saturation estimation results';
        nblines = 1;
        answer = myinputdlg(prompt,title,nblines,{nameDefault},'on');
        if ~isempty(answer)
            name=answer{1};
        else
            return
        end
        if isempty(model.nf_sat)
            no_nf_sat = 1;
        else
            flag=0;
            no_nf_sat = 1+length(model.nf_sat);
            for n=1:length(model.nf_sat)
                if strcmp(model.nf_sat(n).name, name)
                    no_nf_sat = n;
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
        model.nf_sat(no_nf_sat).name = name;
        model.nf_sat(no_nf_sat).S = sat_map;
        model.nf_sat(no_nf_sat).param = nf_sat;

        load([rep,file],'models')
        models(modelNo) = model; %#ok<NASGU>
        save([rep,file],'models','-append')
        saved = true;

        names = cell(1,length(model.nf_sat));
        for n=1:length(model.nf_sat)
            names{n} = model.nf_sat(n).name;
        end
        hsat_map.String = names;

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
        
    end
    function compute(varargin)
        if isempty(model)
            return
        end
        
        no_vel = hdata_vel.Value;
        no_att = hdata_att.Value;
        if ~strcmp(model.tlinv_res(no_vel).param.typeData,'tt')
            errordlg([model.tlinv_res(no_vel).name,' should contain velocity models'])
            return
        end
        if ~strcmp(model.tlinv_res(no_att).param.typeData,'amp')
            errordlg([model.tlinv_res(no_vel).name,' should contain attenuation models'])
            return
        end
        
        data_v0 = 1./model.tlinv_res(no_vel).tomo.s0;
        data_v1 = 1./model.tlinv_res(no_vel).tomo.s;
        data_a0 = model.tlinv_res(no_att).tomo.s0;
        data_a1 = model.tlinv_res(no_att).tomo.s;
        
        nf_sat.vel_data = model.tlinv_res(no_vel).name;
        nf_sat.att_data = model.tlinv_res(no_att).name;

        nf_sat.epsilon_m = str2double(hepsilon_m.String);
        nf_sat.a = str2double(ha.String);
        nf_sat.m = str2double(hm.String);
        nf_sat.sigma_nf = str2double(hsigma_nf.String);
        nf_sat.epsilon_nf = str2double(hepsilon_nf.String);
        nf_sat.mu_nf = str2double(hmu_nf.String);
        nf_sat.sigma_f0 = str2double(hsigma_f0.String);
        nf_sat.epsilon_f0 = str2double(hepsilon_f0.String);
        nf_sat.omega = 2*pi*1e6*str2double(hfreq.String);
        nf_sat.S0 = str2double(hS0.String);
        
        sat_map = nan(size(data_v1));
        
        hh = waitbar(0,'1','Name','Computing saturation ...',...
            'CreateCancelBtn','setappdata(gcf,''canceling'',1)');
        setappdata(gcf,'canceling',0)
        for n=125:130%1:length(sat_map)
            % Check for Cancel button press
            drawnow
            if getappdata(gcf,'canceling')
                break
            end
            
            waitbar(n/length(sat_map), hh, sprintf('Node %d of %d',n, length(sat_map)))
            sat_map(n) = findS(data_v0(n),data_a0(n),data_v1(n),data_a1(n),...
                nf_sat.epsilon_m,nf_sat.a,nf_sat.m,nf_sat.sigma_nf,nf_sat.epsilon_nf,...
                nf_sat.mu_nf,nf_sat.epsilon_f0,nf_sat.sigma_f0,nf_sat.omega,nf_sat.S0);
        end
        delete (hh)

        gridViewer.plotTomo(sat_map,'Nano fluid saturation',...
            'Distance [m]','Elevation [m]',haxes1)
        if hcolorbar.Value==1
            cmin = str2double(hcmin.String);
            cmax = str2double(hcmax.String);
            caxis(haxes1,[cmin cmax])
        end
        colorbar('peer',haxes1)
        colormap(haxes1,hcmap.String{hcmap.Value})
        saved = false;
    end
    function plot_res(varargin)
        
        no_nf_sat = hsat_map.Value;
        figure()
        ax = gca();
        model.nf_sat(no_nf_sat).S
        gridViewer.plotTomo(model.nf_sat(no_nf_sat).S,'Nano fluid saturation',...
            'Distance [m]','Elevation [m]',ax)
        if hcolorbar.Value==1
            cmin = str2double(hcmin.String);
            cmax = str2double(hcmax.String);
            caxis(ax,[cmin cmax])
        end
        colorbar('peer',ax)
        colormap(ax,hcmap.String{hcmap.Value})
        
    end



end