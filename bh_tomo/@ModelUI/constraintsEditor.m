function constraintsEditor(obj,grid,gUI,varargin)

c = grid.cont.copy();
contEdited = false;

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
    'Name','Constraints Editor',...
    'NumberTitle','off',...
    'ToolBar','none',...,
    'MenuBar','None',...
    'SizeChangedFcn',@resizeUI,...
    'CloseRequestFcn',@closeWindow);

enableY = false;
switch grid.type
    case '2D'
        prop = {'Velocity','Attenuation','Reservoir','xi','Tilt Angle'};
    case '2D+'
        prop = {'Velocity','Attenuation','Reservoir','xi','Tilt Angle'};
    case '3D'
        prop = {'Velocity','Attenuation'};  % no time-lapse or anisotropy in 3D
        enableY = true;
end

haxes = axes('Units','points','Parent',f);

pprop = uipanel(f,'Title','Property',...
    'Units','points',...
    'FontSize',fs+1);

hdone =  uicontrol('Style','pushbutton',...
    'String','Done',...
    'Units','points',...
    'FontSize',fs,...
    'Callback',@saveCont,...
    'Parent',f);

hcancel = uicontrol('Style','pushbutton',...
    'String','Cancel',...
    'Units','points',...
    'FontSize',fs,...
    'Callback',@cancel,...
    'Parent',f);

hprop = uicontrol('Style','popupmenu',...
    'String',prop,...
    'Value',1,...
    'Units','points',...
    'FontSize',fs,...
    'Callback',@changeProperty,...
    'Parent',pprop);

hvalue = uicontrol('Style','text',...
    'String','Value',...
    'Units','points',...
    'HorizontalAlignment','right',...
    'FontSize',fs,...
    'Parent',pprop);

hvalueEdit = uicontrol('Style','edit',...
    'String','',...
    'Units','points',...
    'FontSize',fs,...
    'Parent',pprop);

hedit = uicontrol('Style','pushbutton',...
    'String','Edit',...
    'Units','points',...
    'FontSize',fs,...
    'Callback',@edit,...
    'Parent',pprop);

himport = uicontrol('Style','pushbutton',...
    'String','Import',...
    'Units','points',...
    'FontSize',fs,...
    'Callback',@import,...
    'Parent',pprop);

hexport = uicontrol('Style','pushbutton',...
    'String','Export',...
    'Units','points',...
    'FontSize',fs,...
    'Callback',@export,...
    'Parent',pprop);

hreinit = uicontrol('Style','pushbutton',...
    'String','Reinitialize',...
    'Units','points',...
    'FontSize',fs,...
    'Callback',@reinit,...
    'Parent',pprop);

pvariance = uipanel(pprop,'Title','Variance',...
    'Units','points',...
    'FontSize',fs);

uicontrol('Style','text',...
    'String','Value',...
    'Units','normalized',...
    'FontSize',fs,...
    'HorizontalAlignment','right',...
    'Position',[0.05 0.5 0.45 0.32],...
    'Parent',pvariance);

hvarEdit = uicontrol('Style','edit',...
    'String','0',...
    'Units','normalized',...
    'FontSize',fs,...
    'Position',[0.5 0.5 0.45 0.32],...
    'Parent',pvariance);

hvarToggle = uicontrol('Style','togglebutton',...
    'String','Display',...
    'Units','normalized',...
    'FontSize',fs,...
    'Position',[0.05 0.07 0.9 0.32],...
    'Callback',@toggleVariance,...
    'Parent',pvariance);

hcmin = uicontrol('Style','text',...
    'String','cmin',...
    'Units','points',...
    'FontSize',fs,...
    'FontName','FixedWidth',...
    'FontWeight','bold',...
    'Parent',f);

hcminEdit = uicontrol('Style','edit',...
    'String','',...
    'Units','points',...
    'FontSize',fs,...
    'Callback',@clim,...
    'Parent',f);

hcmax = uicontrol('Style','text',...
    'String','cmax',...
    'Units','points',...
    'FontSize',fs,...
    'FontName','FixedWidth',...
    'FontWeight','bold',...
    'Parent',f);

hcmaxEdit = uicontrol('Style','edit',...
    'String','',...
    'Units','points',...
    'FontSize',fs,...
    'Callback',@clim,...
    'Parent',f);

if enableY
    hy = uicontrol('Style','text',...
        'String','Y Coordinate',...
        'Units','points',...
        'FontSize',fs,...
        'Parent',f);
    hypopup = uicontrol('Style','popupmenu',...
        'String',num2str(grid.gry(:)),...
        'Units','points',...
        'FontSize',fs,...
        'Callback',@yaxis,...
        'Parent',f);
end

% plotting variables


gx = [grid.grx(1) grid.grx(2)-grid.grx(1) grid.grx(end)];
gz = [grid.grz(1) grid.grz(2)-grid.grz(1) grid.grz(end)];
dx = gx(2);
dz = gz(2);
xx = (gx(1)+dx/2):dx:(gx(3)-dx/2);
zz = (gz(1)+dz/2):dz:(gz(3)-dz/2);
x=grid.grx(:);
z=grid.grz(:);
aa = [min(x)*ones(length(z),1) max(x)*ones(length(z),1);x x]';
bb = [z z; min(z)*ones(length(x),1) max(z)*ones(length(x),1)]';
if enableY
    gy = [grid.gry(1) grid.gry(2)-grid.gry(1) grid.gry(end)];
    dy = gy(2);
    yy = (gy(1)+dy/2):dy:(gy(3)-dy/2);
end
values = [];
valuesVar = [];
fillValues();
updateFig();

uiwait(f)


    function resizeUI(varargin)
        f.Visible = 'off';
        width = f.Position(3);
        height= f.Position(4);
        
        hBorder = 15;
        
        vFac = 1;
        if ispc
            vFac = 0.81*vFac;
        end
        vSize = 22*vFac;
        vSpace = 5*vFac;
        vBorderTop = 45*vFac;
        vBorder = 15*vFac;
        
        axBorder = haxes.Position(2);
        hSize = 180;
        hSize2 = hSize-14;
        
        haxes.Position = [axBorder axBorder width-3*axBorder-hSize-hBorder height-2*axBorder];
        hcmin.Position = [width-axBorder-hSize-2*hBorder axBorder-vSpace-vSize axBorder vSize];
        hcminEdit.Position = [width-axBorder-hSize-2*hBorder axBorder axBorder vSize];
        hcmax.Position = [width-axBorder-hSize-2*hBorder height-axBorder+vSpace axBorder vSize];
        hcmaxEdit.Position = [width-axBorder-hSize-2*hBorder height-axBorder-vSize axBorder vSize];
        
        hdone.Position = [width-hBorder-hSize vBorder hSize vSize];
        hcancel.Position = [width-hBorder-hSize vBorder+vSize+vSpace hSize vSize];
        
        vSizeVar = 2*vSize+7*vSpace;
        vSizeProp = 11*vSpace+vSizeVar+6*vSize;
        pprop.Position = [width-hBorder-hSize height-vSizeProp-vBorderTop hSize vSizeProp];
        
        pvariance.Position  = [5 vSpace hSize2 vSizeVar];
        hreinit.Position    = [5 2*vSpace+vSizeVar hSize2 vSize];
        hexport.Position    = [5 3*vSpace+vSizeVar+vSize hSize2 vSize];
        himport.Position    = [5 4*vSpace+vSizeVar+2*vSize hSize2 vSize];
        hedit.Position      = [5 5*vSpace+vSizeVar+3*vSize hSize2 vSize];
        hvalue.Position     = [5 6*vSpace+vSizeVar+4*vSize hSize2/2 vSize];
        hvalueEdit.Position = [hSize/2 6*vSpace+vSizeVar+4*vSize hSize2/2 vSize];
        hprop.Position      = [5 7*vSpace+vSizeVar+5*vSize hSize2 vSize];
        
        if enableY
            hy.Position = [width-hBorder-hSize height-vSizeProp-vBorderTop-vBorder-vSize hSize vSize];
            hypopup.Position = [width-hBorder-hSize height-vSizeProp-vBorderTop-vBorder-2*vSize hSize vSize];
        end
        
        f.Visible = 'on';
    end
    function closeWindow(varargin)
        if contEdited
            choice = questdlg('Do you want to save the constraints?',...
                'Constraints',...
                'Yes','no','Cancel','Cancel');
            switch choice
                case 'Yes'
                    grid.cont = c.copy();
                    obj.notify('modelEdited')
                case 'Cancel'
                    return
            end
        end
        delete(f)
    end
    function saveCont(varargin)
        if contEdited
            grid.cont = c.copy();
            obj.notify('modelEdited')
        end
        delete(f)
    end
    function cancel(varargin)
        delete(f)
    end
    function changeProperty(varargin)
        fillValues()
        updateFig()
    end
    function toggleVariance(varargin)
        updateFig()
    end
    function edit(varargin)
        if isempty(hvalueEdit.String)
            warndlg('Define first value of property')
            return
        end
        val = str2double(hvalueEdit.String);
        if hprop.Value == 3 && val == 0
            % for reservoir
            val = NaN;
        end
        if ~isnumeric(val)
            warndlg('Property value should be numeric')
            return            
        end
        if isempty(hvarEdit.String)
            warndlg('Define value of variance')
            return
        end
        var = str2double(hvarEdit.String);
        if ~isnumeric(val)
            warndlg('Variance value should be numeric')
            return            
        end
        
        axes(haxes)
        b=1;
        while b==1
            [x,z,b] = ginput(1);
            if x<grid.grx(1) || x>grid.grx(end) || z<grid.grz(1) || z>grid.grz(end)
                break
            end
            ix = findnear(x,xx);
            iz = findnear(z,zz);
            values(iz,ix) = val;
            valuesVar(iz,ix) = var;
            updateFig();
        end
        fillCont();        
        contEdited = true;
    end
    function import(varargin)
        switch hprop.Value
            case 1
                str = 'File with Velocity Constraints';
            case 2
                str = 'File with Attenuation Constraints';
            case 3
                str = 'File with Reservoir Cells';
            case 4
                str = 'File with xi Constraints';
            case 5
                str = 'File with Tilt Angle Constraints';
        end
        
        [file, rep] = uigetfile('*.con',str);
        if file==0
            return
        end
        icont = load([rep,file]);
        if size(icont,2)~=4 && size(icont,2)~=5
            errordlg('Constaints file should contain X Y Z value [variance]')
            return
        end
        if size(icont,2)==4
            % add zero variance if not defined
            icont = [icont zeros(size(icont,1))];
        end
        switch grid.type
            case '2D'
                % Project on plane
                
                [coord,dist] = gUI.project(icont(:,1:3));
                fd=figure;
                plot(dist,'*')
                xlabel('Point Number')
                ylabel('Distance')
                title('Distance of projected point from grid plane')
                
                cutoff = myinputdlg('Enter cutoff distance');
                close(fd)
                if ~isempty(cutoff)
                    ind = dist<str2double(cutoff);
                else
                    ind = 1:length(dist);
                end
                icont = [coord(ind,1) coord(ind,3) icont(ind,4:5)];
                
                switch hprop.Value
                    case 1
                        c.slowness.data = [c.slowness.data; ...
                            [icont(:,1:2) 1./icont(:,3) ...
                            icont(:,4)./icont(:,3).^4]];
                        c.slowness.data = unique(c.slowness.data,'rows');
                    case 2
                        c.attenuation.data = [c.attenuation.data; icont];
                        c.attenuation.data = unique(c.attenuation.data,'rows');
                    case 3
                        ix = zeros(size(icont,1),1);
                        iz = zeros(size(icont,1),1);
                        for n=1:size(icont,1)
                            ix(n) = findnear(icont(n,1),xx);
                            iz(n) = findnear(icont(n,2),zz);
                        end
                        ix = unique(ix);
                        iz = unique(iz);
                        nz=length(zz);
                        c.ind_reservoir = (ix-1)*nz+iz;
                    case 4
                        c.slowness.data_xi = [c.slowness.data_xi; icont];
                        c.slowness.data_xi = unique(c.slowness.data_xi,'rows');
                    case 5
                        c.slowness.data_theta = [c.slowness.data_theta; icont];
                        c.slowness.data_theta = unique(c.slowness.data_theta,'rows');
                end
                
            case '2D+'
                % Project on planes
                [coord,dist] = gUI.project(icont(:,1:3));
                fd=figure;
                plot(dist,'*')
                xlabel('Point Number')
                ylabel('Distance')
                title('Distance of projected point from grid plane')
                
                cutoff = myinputdlg('Enter cutoff distance');
                close(fd)
                if ~isempty(cutoff)
                    ind = dist<str2double(cutoff);
                else
                    ind = 1:length(dist);
                end
                icont = [coord(ind,1) coord(ind,2) icont(ind,4:5)];
                
                switch hprop.Value
                    case 1
                        c.slowness.data = [c.slowness.data;...
                            [icont(:,1:2) 1./icont(:,3) ...
                            icont(:,4)./icont(:,3).^4]];
                        c.slowness.data = unique(c.slowness.data,'rows');
                    case 2
                        c.attenuation.data = [c.attenuation.data; icont];
                        c.attenuation.data = unique(c.attenuation.data,'rows');
                    case 3
                        ix = zeros(size(icont,1),1);
                        iz = zeros(size(icont,1),1);
                        for n=1:size(icont,1)
                            ix(n) = findnear(icont(n,1),xx);
                            iz(n) = findnear(icont(n,2),zz);
                        end
                        ix = unique(ix);
                        iz = unique(iz);
                        nz=length(zz);
                        c.ind_reservoir = (ix-1)*nz+iz;
                    case 4
                        c.slowness.data_xi = [c.slowness.data_xi; icont];
                        c.slowness.data_xi = unique(c.slowness.data_xi,'rows');
                    case 5
                        c.slowness.data_theta = [c.slowness.data_theta; icont];
                        c.slowness.data_theta = unique(c.slowness.data_theta,'rows');
                end
                
            case '3D'
                switch hprop.Value
                    case 1
                        c.slowness.data = [c.slowness.data;...
                            [icont(:,1:3) 1./icont(:,4) ...
                            icont(:,5)./icont(:,4).^4]];
                        c.slowness.data = unique(c.slowness.data,'rows');
                    case 2
                        c.attenuation.data = [c.attenuation.data; icont];
                        c.attenuation.data = unique(c.attenuation.data,'rows');
                    case 3
                        % TODO
                    case 4
                        c.slowness.data_xi = [c.slowness.data_xi; icont];
                        c.slowness.data_xi = unique(c.slowness.data_xi,'rows');
                    case 5
                        c.slowness.data_theta = [c.slowness.data_theta; icont];
                        c.slowness.data_theta = unique(c.slowness.data_theta,'rows');
                end
        end
        
        fillValues()
        updateFig()
        
        contEdited = true;
    end
    function export(varargin)
        
        switch grid.type
            case '2D'
                % TODO we must "reverse project" to get coordonates in 3D
                warndlg('Not yet working in 2D')
            case '2D+'
                % TODO we must "reverse project" to get coordonates in 3D
                warndlg('Not yet working in 2D')
            case '3D'
                
                switch hprop.Value
                    case 1
                        str = 'File with Velocity Constraints';
                    case 2
                        str = 'File with Attenuation Constraints';
                    case 3
                        str = 'File with Reservoir Cells';
                    case 4
                        str = 'File with xi Constraints';
                    case 5
                        str = 'File with Tilt Angle Constraints';
                end
                
                [file, rep] = uiputfile('*.con',str);
                if file==0
                    return
                end
                
                if hprop.Value==3
                    [iz,ix]=find(~isnan(values));
                    xt = xx(ix);
                    zt = zz(iz);
                    cont = [xt(:) zt(:) ones(numel(ix),1) zeros(numel(ix),1)]; %#ok<NASGU>
                    save([rep,file],'cont','-ascii')
                else
                    fillCont()
                    switch hprop.Value
                        case 1
                            cont = [c.slowness.data(:,1:2) 1./c.slowness.data(:,3) ...
                                c.slowness.data(:,4)./c.slowness.data(:,3).^4]; %#ok<NASGU>
                            save([rep,file],'cont','-ascii')
                        case 2
                            save([rep,file],'c.attenuation.data','-ascii')
                        case 4
                            save([rep,file],'c.slowness.data_xi','-ascii')
                        case 5
                            save([rep,file],'c.slowness.data_theta','-ascii')
                    end
                end
        end
    end
    function reinit(varargin)
        c = grid.cont;
        fillValues();
        updateFig();
    end
    function fillValues()
        values = nan(length(zz), length(xx));
        valuesVar = nan(length(zz), length(xx));
        
        if enableY==false  % 2D
            if hprop.Value == 3  % reservoir (we don't care about variance)
                if ~isempty(c.ind_reservoir)
                    nz=length(zz);
                    ix = floor(0.0001+c.ind_reservoir/nz)+1;
                    iz = c.ind_reservoir - (ix-1)*nz;
                    for n=1:length(ix)
                        values(iz(n),ix(n)) = 1;
                    end
                end
                return
            end
            
            switch hprop.Value
                case 1  % velocity
                    if ~isempty(c.slowness.data)
                        for n=1:size(c.slowness.data,1)
                            if c.slowness.data(n,1)>=xx(1) && ...
                                    c.slowness.data(n,1)<=xx(end) && ...
                                    c.slowness.data(n,2)>=zz(1) && ...
                                    c.slowness.data(n,2)<=zz(end)
                                ix = findnear(c.slowness.data(n,1), xx);
                                iz = findnear(c.slowness.data(n,2), zz);
                                values(iz,ix) = 1/c.slowness.data(n,3);
                                % http://math.stackexchange.com/questions/269216/inverse-of-random-variable
                                valuesVar(iz,ix) = c.slowness.data(n,4)/c.slowness.data(n,3)^4;
                            end
                        end
                    end
                case 2  % attenuation
                    if ~isempty(c.attenuation.data)
                        for n=1:size(c.attenuation.data,1)
                            if c.attenuation.data(n,1)>=xx(1) && ...
                                    c.attenuation.data(n,1)<=xx(end) && ...
                                    c.attenuation.data(n,2)>=zz(1) && ...
                                    c.attenuation.data(n,2)<=zz(end)
                                ix = findnear(c.attenuation.data(n,1), xx);
                                iz = findnear(c.attenuation.data(n,2), zz);
                                values(iz,ix) = c.attenuation.data(n,3);
                                valuesVar(iz,ix) = c.attenuation.data(n,4);
                            end
                        end
                    end
                case 4  % xi
                    if ~isempty(c.slowness.data_xi)
                        for n=1:size(c.slowness.data_xi,1)
                            if c.slowness.data_xi(n,1)>=xx(1) && ...
                                    c.slowness.data_xi(n,1)<=xx(end) && ...
                                    c.slowness.data_xi(n,2)>=zz(1) && ...
                                    c.slowness.data_xi(n,2)<=zz(end)
                                ix = findnear(c.slowness.data_xi(n,1), xx);
                                iz = findnear(c.slowness.data_xi(n,2), zz);
                                values(iz,ix) = c.slowness.data_xi(n,3);
                                valuesVar(iz,ix) = c.slowness.data_xi(n,4);
                            end
                        end
                    end
                case 5  % tilt angle
                    if ~isempty(c.slowness.data_theta)
                        for n=1:size(c.slowness.data_theta,1)
                            if c.slowness.data_theta(n,1)>=xx(1) && ...
                                    c.slowness.data_theta(n,1)<=xx(end) && ...
                                    c.slowness.data_theta(n,2)>=zz(1) && ...
                                    c.slowness.data_theta(n,2)<=zz(end)
                                ix = findnear(c.slowness.data_theta(n,1), xx);
                                iz = findnear(c.slowness.data_theta(n,2), zz);
                                values(iz,ix) = c.slowness.data_theta(n,3);
                                valuesVar(iz,ix) = c.slowness.data_theta(n,4);
                            end
                        end
                    end
            end
        else % 3D
            if hprop.Value == 3  % reservoir (we don't care about variance)
                if ~isempty(c.ind_reservoir)
                    for n=1:size(c.ind_reservoir,1)
                        if c.ind_reservoir(n,2)>=yy(1) && c.ind_reservoir(n,2)<=yy(end)
                            iy = findnear(c.ind_reservoir(n,2), yy);
                            if iy==hypopup.Value
                                if c.ind_reservoir(n,1)>=xx(1) && ...
                                        c.ind_reservoir(n,1)<=xx(end) && ...
                                        c.ind_reservoir(n,3)>=zz(1) && ...
                                        c.ind_reservoir(n,3)<=zz(end)
                                    ix = findnear(c.ind_reservoir(n,1), xx);
                                    iz = findnear(c.ind_reservoir(n,3), zz);
                                    values(iz,ix) = c.ind_reservoir(n,4);
                                end
                            end
                        end
                    end
                end
                return
            end
            
            switch hprop.Value
                case 1  % velocity
                    if ~isempty(c.slowness.data)
                        for n=1:size(c.slowness.data,1)
                            if c.slowness.data(n,2)>=yy(1) && c.slowness.data(n,2)<=yy(end)
                                iy = findnear(c.slowness.data(n,2), yy);
                                if iy==hypopup.Value
                                    if c.slowness.data(n,1)>=xx(1) && ...
                                            c.slowness.data(n,1)<=xx(end) && ...
                                            c.slowness.data(n,3)>=zz(1) && ...
                                            c.slowness.data(n,3)<=zz(end)
                                        ix = findnear(c.slowness.data(n,1), xx);
                                        iz = findnear(c.slowness.data(n,3), zz);
                                        values(iz,ix) = 1/c.slowness.data(n,4);
                                        % http://math.stackexchange.com/questions/269216/inverse-of-random-variable
                                        valuesVar(iz,ix) = c.slowness.data(n,5)/c.slowness.data(n,4)^4;
                                    end
                                end
                            end
                        end
                    end
                case 2  % attenuation
                    if ~isempty(c.attenuation.data)
                        for n=1:size(c.attenuation.data,1)
                            if c.attenuation.data(n,2)>=yy(1) && c.attenuation.data(n,2)<=yy(end)
                                iy = findnear(c.attenuation.data(n,2), yy);
                                if iy==hypopup.Value
                                    if c.attenuation.data(n,1)>=xx(1) && ...
                                            c.attenuation.data(n,1)<=xx(end) && ...
                                            c.attenuation.data(n,3)>=zz(1) && ...
                                            c.attenuation.data(n,3)<=zz(end)
                                        ix = findnear(c.attenuation.data(n,1), xx);
                                        iz = findnear(c.attenuation.data(n,3), zz);
                                        values(iz,ix) = c.attenuation.data(n,4);
                                        valuesVar(iz,ix) = c.attenuation.data(n,5);
                                    end
                                end
                            end
                        end
                    end
                case 4  % xi
                    if ~isempty(c.slowness.data_xi)
                        for n=1:size(c.slowness.data_xi,1)
                            if c.slowness.data_xi(n,2)>=yy(1) && c.slowness.data_xi(n,2)<=yy(end)
                                iy = findnear(c.slowness.data_xi(n,2), yy);
                                if iy==hypopup.Value
                                    if c.slowness.data_xi(n,1)>=xx(1) && ...
                                            c.slowness.data_xi(n,1)<=xx(end) && ...
                                            c.slowness.data_xi(n,3)>=zz(1) && ...
                                            c.slowness.data_xi(n,3)<=zz(end)
                                        ix = findnear(c.slowness.data_xi(n,1), xx);
                                        iz = findnear(c.slowness.data_xi(n,3), zz);
                                        values(iz,ix) = c.slowness.data_xi(n,4);
                                        valuesVar(iz,ix) = c.slowness.data_xi(n,5);
                                    end
                                end
                            end
                        end
                    end
                case 5  % tilt angle
                    if ~isempty(c.slowness.data_theta)
                        for n=1:size(c.slowness.data_theta,1)
                            if c.slowness.data_theta(n,2)>=yy(1) && c.slowness.data_theta(n,2)<=yy(end)
                                iy = findnear(c.slowness.data_theta(n,2), yy);
                                if iy==hypopup.Value
                                    if c.slowness.data_theta(n,1)>=xx(1) && ...
                                            c.slowness.data_theta(n,1)<=xx(end) && ...
                                            c.slowness.data_theta(n,3)>=zz(1) && ...
                                            c.slowness.data_theta(n,3)<=zz(end)
                                        ix = findnear(c.slowness.data_theta(n,1), xx);
                                        iz = findnear(c.slowness.data_theta(n,3), zz);
                                        values(iz,ix) = c.slowness.data_theta(n,4);
                                        valuesVar(iz,ix) = c.slowness.data_theta(n,5);
                                    end
                                end
                            end
                        end
                    end
            end
        end
    end
    function fillCont()
        [iz,ix]=find(~isnan(values));
        v=zeros(numel(ix),1);
        vv=v;
        for n=1:numel(ix)
            v(n) = values(iz(n),ix(n));
            vv(n)= valuesVar(iz(n),ix(n));
        end
        if enableY==false  % 2D
            switch hprop.Value
                case 1
                    c.slowness.data = [xx(ix)' zz(iz)' 1./v vv./v.^4];
                case 2
                    c.attenuation.data = [xx(ix)' zz(iz)' v vv];
                case 3
                    nz=length(zz);
                    c.ind_reservoir = (ix-1)*nz+iz;
                case 4
                    c.slowness.data_xi = [xx(ix)' zz(iz)' v vv];
                case 5
                    c.slowness.data_theta = [xx(ix)' zz(iz)' v vv];
            end
        else  %3D
            y = grid.gry(hypopup.Value);
            switch hprop.Value
                case 1
                    for n=1:numel(ix)
                        ind = find(sum(c.slowness.data(:,1:3)==kron([xx(ix(n)) y zz(iz(n))],...
                            ones(size(c.slowness.data,1),1)),2)==3,1);
                        if isempty(ind)
                            c.slowness.data = [c.slowness.data;
                                [xx(ix(n)) y zz(iz(n)) 1/values(iz(n),ix(n)) ...
                                valuesVar(iz(n),ix(n))./values(iz(n),ix(n))^4]];
                        else
                            c.slowness.data(ind,3) = 1/values(iz(n),ix(n));
                            c.slowness.data(ind,4) = valuesVar(iz(n),ix(n))./values(iz(n),ix(n))^4;
                        end
                    end
                        
                case 2
                    for n=1:numel(ix)
                        ind = find(sum(c.attenuation.data(:,1:3)==kron([xx(ix(n)) y zz(iz(n))],...
                            ones(size(c.attenuation.data,1),1)),2)==3,1);
                        if isempty(ind)
                            c.attenuation.data = [c.attenuation.data;
                                [xx(ix(n)) y zz(iz(n)) values(iz(n),ix(n)) ...
                                valuesVar(iz(n),ix(n))]];
                        else
                            c.attenuation.data(ind,3) = values(iz(n),ix(n));
                            c.attenuation.data(ind,4) = valuesVar(iz(n),ix(n));
                        end
                    end
                case 3
                    for n=1:numel(ix)
                        ind = find(sum(c.ind_reservoir(:,1:3)==kron([xx(ix(n)) y zz(iz(n))],...
                            ones(size(c.ind_reservoir,1),1)),2)==3,1);
                        if isempty(ind)
                            c.ind_reservoir = [c.ind_reservoir;
                                [xx(ix(n)) y zz(iz(n)) values(iz(n),ix(n)) ...
                                valuesVar(iz(n),ix(n))]];
                        else
                            c.ind_reservoir(ind,3) = values(iz(n),ix(n));
                            c.ind_reservoir(ind,4) = valuesVar(iz(n),ix(n));
                        end
                    end
                case 4
                    for n=1:numel(ix)
                        ind = find(sum(c.slowness.data_xi(:,1:3)==kron([xx(ix(n)) y zz(iz(n))],...
                            ones(size(c.slowness.data_xi,1),1)),2)==3,1);
                        if isempty(ind)
                            c.slowness.data_xi = [c.slowness.data_xi;
                                [xx(ix(n)) y zz(iz(n)) values(iz(n),ix(n)) ...
                                valuesVar(iz(n),ix(n))]];
                        else
                            c.slowness.data_xi(ind,3) = values(iz(n),ix(n));
                            c.slowness.data_xi(ind,4) = valuesVar(iz(n),ix(n));
                        end
                    end
                case 5
                    for n=1:numel(ix)
                        ind = find(sum(c.slowness.data_theta(:,1:3)==kron([xx(ix(n)) y zz(iz(n))],...
                            ones(size(c.slowness.data_theta,1),1)),2)==3,1);
                        if isempty(ind)
                            c.slowness.data_theta = [c.slowness.data_theta;
                                [xx(ix(n)) y zz(iz(n)) values(iz(n),ix(n)) ...
                                valuesVar(iz(n),ix(n))]];
                        else
                            c.slowness.data_theta(ind,3) = values(iz(n),ix(n));
                            c.slowness.data_theta(ind,4) = valuesVar(iz(n),ix(n));
                        end
                    end
            end
        end
    end
    function updateFig()
        
        if hvarToggle.Value == 0 || hprop.Value == 3
            imagesc(xx,zz,values,'Parent',haxes)
        else
            imagesc(xx,zz,valuesVar,'Parent',haxes)
        end
        hold(haxes,'on')
        plot(haxes,aa,bb,'Color',[0.5 0.5 0.5])
        set(haxes,'YDir','normal')
        hold(haxes,'off')
        axis(haxes,'equal')
        axis(haxes,'tight')
        xlabel(haxes,'X')
        ylabel(haxes,'Z')
        colorbar('peer',haxes);
        hcminEdit.String = num2str(haxes.CLim(1));
        hcmaxEdit.String = num2str(haxes.CLim(2));
    end
    function yaxis(varargin)
        fillValues();
        updateFig();
    end

    function clim(varargin)
        if ~isempty(hcmaxEdit.String)
            haxes.CLim(2) = str2double(hcmaxEdit.String);
        end
        if ~isempty(hcminEdit.String)
            haxes.CLim(1) = str2double(hcminEdit.String);
        end
    end
end