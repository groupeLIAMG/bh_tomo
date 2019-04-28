function prune(obj,varargin)

no = obj.handles.listMogs.Value;
if no<=0 || no>length(obj.mogs)
    warndlg('No MOG selected')
    return
end
fs = 11;
if isa(varargin{1},'matlab.ui.control.UIControl')
    fs = varargin{1}.FontSize;
elseif ispc
    fs = 9;
end
vScale = 1;
if ispc
    vScale = 0.81;
end
mog = obj.mogs(no);
SNR = [];
width = 700;
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
    'Tag','fig_bh_tomo2_prune',...
    'Name','Prune MOG',...
    'NumberTitle','off',...
    'ToolBar','none',...,
    'MenuBar','None',...
    'SizeChangedFcn',@resizeUI,...
    'CloseRequestFcn',@closeWindow);

hskipTx = uicontrol('Style','text',...
    'Units','points',...
    'String','Number of Stations to Skip - Tx',...
    'Fontsize',fs,...
    'HorizontalAlignment','center',...
    'Parent',f);
hskipTxEdit = uicontrol('Style','edit',...
    'Units','points',...
    'String',num2str(mog.pruneParams.sautTx),...
    'Fontsize',fs,...
    'Callback',@update,...
    'Parent',f);

hskipRx = uicontrol('Style','text',...
    'Units','points',...
    'String','Number of Stations to Skip - Rx',...
    'Fontsize',fs,...
    'HorizontalAlignment','center',...
    'Parent',f);
hskipRxEdit = uicontrol('Style','edit',...
    'Units','points',...
    'String',num2str(mog.pruneParams.sautRx),...
    'Fontsize',fs,...
    'Callback',@update,...
    'Parent',f);

hround = uicontrol('Style','text',...
    'Units','points',...
    'String','Rounding Factor',...
    'Fontsize',fs,...
    'HorizontalAlignment','center',...
    'Parent',f);
hroundEdit = uicontrol('Style','edit',...
    'Units','points',...
    'String',num2str(mog.pruneParams.arrondi),...
    'Fontsize',fs,...
    'Callback',@update,...
    'Parent',f);

hsnr = uicontrol('Style','checkbox',...
    'Units','points',...
    'String','Threshhold - SNR',...
    'Value',mog.pruneParams.use_SB,...
    'Fontsize',fs,...
    'Callback',@update,...
    'Parent',f);
hsnrEdit = uicontrol('Style','edit',...
    'Units','points',...
    'String',num2str(mog.pruneParams.seuil_SB),...
    'Fontsize',fs,...
    'Callback',@update,...
    'Parent',f);

htmin = uicontrol('Style','text',...
    'Units','points',...
    'String','Minimum Ray Angle',...
    'Fontsize',fs,...
    'HorizontalAlignment','center',...
    'Parent',f);
htminEdit = uicontrol('Style','edit',...
    'Units','points',...
    'String',num2str(mog.pruneParams.thetaMin),...
    'Fontsize',fs,...
    'Callback',@update,...
    'Parent',f);

htmax = uicontrol('Style','text',...
    'Units','points',...
    'String','Maximum Ray Angle',...
    'Fontsize',fs,...
    'HorizontalAlignment','center',...
    'Parent',f);
htmaxEdit = uicontrol('Style','edit',...
    'Units','points',...
    'String',num2str(mog.pruneParams.thetaMax),...
    'Fontsize',fs,...
    'Callback',@update,...
    'Parent',f);

hzmin = uicontrol('Style','text',...
    'Units','points',...
    'String','Minimum Elevation',...
    'Fontsize',fs,...
    'HorizontalAlignment','center',...
    'Parent',f);
hzminEdit = uicontrol('Style','edit',...
    'Units','points',...
    'String',num2str(mog.pruneParams.zmin),...
    'Fontsize',fs,...
    'Callback',@update,...
    'Parent',f);

hzmax = uicontrol('Style','text',...
    'Units','points',...
    'String','Maximum Elevation',...
    'Fontsize',fs,...
    'HorizontalAlignment','center',...
    'Parent',f);
hzmaxEdit = uicontrol('Style','edit',...
    'Units','points',...
    'String',num2str(mog.pruneParams.zmax),...
    'Fontsize',fs,...
    'Callback',@update,...
    'Parent',f);

hzeros = uicontrol('Style','checkbox',...
    'Units','points',...
    'String','Ignore traces filled with 0',...
    'Value',mog.pruneParams.skip_zeros,...
    'Fontsize',fs,...
    'Callback',@update,...
    'Parent',f);

hinfo = uicontrol('Style','text',...
    'Units','points',...
    'String','Infos',...
    'Fontsize',fs,...
    'BackgroundColor',[1 1 1],...
    'HorizontalAlignment','center',...
    'Parent',f);

hdone = uicontrol('Style','pushbutton',...
    'Units','points',...
    'String','Done',...
    'Fontsize',fs,...
    'Callback',@closeWindow,...
    'Parent',f);


% next two buttons are off until I find a way get ride of the out of memory
% problem
hzoom = uicontrol('Style','togglebutton',...
    'CData',zoomCData,...
    'Callback',@zoom,...
    'Visible','off',...  
    'Parent',f);
hrotate3d = uicontrol('Style','togglebutton',...
    'CData',rotateCData,...
    'Callback',@rotate3d,...
    'Visible','off',...
    'Parent',f);

haxes = axes('Units','points','Parent',f);

update()
uiwait(f)

    function resizeUI(varargin)
        f.Visible = 'off';
        
        hBorderRight = 25;
        width = f.Position(3);
        height = f.Position(4);
        
        hSize = 77;
        
        vFac = 0.8*height/750;
        if vFac<1
            vFac = 1;
        end
        if ispc
            vFac = 0.81*vFac;
        end
        vSize = 22*vFac;
        vSpace = 5*vFac;
        vBorderTop = 45*vFac;
        vBorder = 20*vFac;
        
        hdone.Position = [width-hBorderRight-2*hSize vBorder 2*hSize vSize];
        hzoom.Position = [width-hBorderRight-10/7*hSize vBorder+vSize+2*vSpace 22 22];
        hrotate3d.Position = [width-hBorderRight-6/7*hSize vBorder+vSize+2*vSpace 22 22];
        
        pos = haxes.Position;
        haxes.Position = [pos(1) pos(2) width-3*hSize-pos(1) height-vBorderTop-pos(2)];
        
        ext = hskipTx.Extent;
        hSize2 = 1.2*ext(3);
        hskipTx.Position = [width-hBorderRight-hSize-hSize2/2 height-vBorderTop-vSize hSize2 vSize];
        hskipTxEdit.Position = [width-hBorderRight-3/2*hSize height-vBorderTop-2*vSize hSize vSize];
        hskipRx.Position = [width-hBorderRight-hSize-hSize2/2 height-vBorderTop-3*vSize-vSpace hSize2 vSize];
        hskipRxEdit.Position = [width-hBorderRight-3/2*hSize height-vBorderTop-4*vSize-vSpace hSize vSize];
        hround.Position = [width-hBorderRight-hSize-hSize2/2 height-vBorderTop-5*vSize-2*vSpace hSize2 vSize];
        hroundEdit.Position = [width-hBorderRight-3/2*hSize height-vBorderTop-6*vSize-2*vSpace hSize vSize];
        
        ext = hsnr.Extent;
        hSize2 = 1.3*ext(3);
        hsnr.Position = [width-hBorderRight-hSize-hSize2/2 height-vBorderTop-7*vSize-3*vSpace hSize2 vSize];
        hsnrEdit.Position = [width-hBorderRight-3/2*hSize height-vBorderTop-8*vSize-3*vSpace hSize vSize];
        htmin.Position = [width-hBorderRight-hSize-hSize2/2 height-vBorderTop-9*vSize-4*vSpace hSize2 vSize];
        htminEdit.Position = [width-hBorderRight-3/2*hSize height-vBorderTop-10*vSize-4*vSpace hSize vSize];
        htmax.Position = [width-hBorderRight-hSize-hSize2/2 height-vBorderTop-11*vSize-5*vSpace hSize2 vSize];
        htmaxEdit.Position = [width-hBorderRight-3/2*hSize height-vBorderTop-12*vSize-5*vSpace hSize vSize];
        hzmin.Position = [width-hBorderRight-hSize-hSize2/2 height-vBorderTop-13*vSize-6*vSpace hSize2 vSize];
        hzminEdit.Position = [width-hBorderRight-3/2*hSize height-vBorderTop-14*vSize-6*vSpace hSize vSize];
        hzmax.Position = [width-hBorderRight-hSize-hSize2/2 height-vBorderTop-15*vSize-7*vSpace hSize2 vSize];
        hzmaxEdit.Position = [width-hBorderRight-3/2*hSize height-vBorderTop-16*vSize-7*vSpace hSize vSize];
        
        hzeros.Position(3) = 3*hSize;
        
        hinfo.Position = [width-hBorderRight-2*hSize vBorder+2*vSize+5*vSpace 2*hSize 7*vSize];
        
        f.Visible = 'on';
    end
    function closeWindow(varargin)
        delete(f)
    end
    function update(varargin)
        skipTx = str2double(hskipTxEdit.String);
        skipRx = str2double(hskipRxEdit.String);
        zmin = str2double(hzminEdit.String);
        zmax = str2double(hzmaxEdit.String);
        thetaMin = str2double(htminEdit.String)*pi/180;
        thetaMax = str2double(htmaxEdit.String)*pi/180;
        
        inTx = mog.data.Tx_z >= zmin & mog.data.Tx_z <= zmax;
        inRx = mog.data.Rx_z >= zmin & mog.data.Rx_z <= zmax;
        
        if ~any(inTx) || ~any(inRx)
            zmin = min([mog.data.Tx_z mog.data.Rx_z]);
            zmax = max([mog.data.Tx_z mog.data.Rx_z]);
            hzminEdit.String = num2str(zmin);
            hzmaxEdit.String = num2str(zmax);
            inTx = mog.data.Tx_z >= zmin & mog.data.Tx_z <= zmax;
            inRx = mog.data.Rx_z >= zmin & mog.data.Rx_z <= zmax;
        end
        
        Tx = [mog.data.Tx_x; mog.data.Tx_y; mog.data.Tx_z]';
        Rx = [mog.data.Rx_x; mog.data.Rx_y; mog.data.Rx_z]';
        
        dr = sqrt( (mog.data.Tx_x-mog.data.Rx_x).^2 +...
            (mog.data.Tx_y-mog.data.Rx_y).^2 );
        theta = atan2(mog.data.Rx_z-mog.data.Tx_z, dr);
        
        inTheta = theta>thetaMin & theta<thetaMax;

        arrondi = str2double(hroundEdit.String);
        if arrondi > 0
            Tx = arrondi*(round(Tx/arrondi));
            Rx = arrondi*(round(Rx/arrondi));
        end
        uTx = unique(Tx(inTx,:), 'rows');
        uRx = unique(Rx(inRx,:), 'rows');
        
        [~,n] = max(var(uTx));
        [~,n] = sort(uTx(:,n));
        uTx = uTx(n,:);
        [~,n] = max(var(uRx));
        [~,n] = sort(uRx(:,n));
        uRx = uRx(n,:);
        
        uTx = uTx(1:(1+skipTx):end,:);
        uRx = uRx(1:(1+skipRx):end,:);

        inTx = false(1,mog.data.ntrace);
        nTx = size(uTx,1);
        for n1=1:length(inTx)
            inTx(n1) = ~isempty(find(all(repmat(Tx(n1,:),nTx,1)==uTx,2),1));
        end
        
        inRx = false(1,mog.data.ntrace);
        nRx = size(uRx,1);
        for n1=1:length(inRx)
            inRx(n1) = ~isempty(find(all(repmat(Rx(n1,:),nRx,1)==uRx,2),1));
        end
        
        use_SB = hsnr.Value;
        skip_zeros = hzeros.Value;
        
        if isempty(SNR)
            SNR = computeSNR();
        end
        seuil_SB = str2double(hsnrEdit.String);
        
        if skip_zeros == 1
            in_zeros = sum(mog.data.rdata == 0, 1)~=mog.data.nptsptrc;
        else
            in_zeros = true(1,mog.data.ntrace);
        end
        mog.in = inTx & inRx & SNR>seuil_SB & inTheta & in_zeros;
        
        
        [~, a]=Grid.lsplane([uTx;uRx]);
        el = (pi-a(3))*180/pi;
        az = atan( cos(a(2))/cos(a(1)) )*180/pi;
        
        texte = cell(8,1);
        texte{1} = 'Infos';
        texte{2} = '';
        texte{3} = [num2str(size(uTx,1)), ' Tx'];
        texte{4} = [num2str(size(uRx,1)), ' Rx'];
        texte{5} = [num2str(round(100*sum(~(inTx & inRx))/length(mog.in)),'%d'), '% removed - Tx & Rx'];
        texte{6} = [num2str(round(100*sum(SNR<=seuil_SB)/length(mog.in)),'%d'), '% removed - S/M ratio'];
        texte{7} = [num2str(round(100*sum(~inTheta)/length(mog.in)),'%d'), '% removed - ray angle'];
        texte{8} = [num2str(round(100*sum(~(in_zeros))/length(mog.in)),'%d'), '% removed - zeros'];
        texte{9} = [num2str(round(100*sum(mog.in)/length(mog.in)),'%d'), '% of traces kept'];
        hinfo.String = texte;
        
        plot3(haxes,uTx(:,1), uTx(:,2), uTx(:,3),'b*')
        hold(haxes,'on')
        plot3(haxes,uRx(:,1), uRx(:,2), uRx(:,3),'gx')
        haxes.DataAspectRatio = [1 1 1];
        grid(haxes,'on')
        xlabel(haxes,'X')
        ylabel(haxes,'Y')
        zlabel(haxes,'Z')
        view(haxes,az,el)
        legend(haxes,'Tx','Rx','Location','Best')
        hold(haxes,'off')
        
        mog.pruneParams.sautTx = skipTx;
        mog.pruneParams.sautRx = skipRx;
        mog.pruneParams.arrondi = arrondi;
        mog.pruneParams.use_SB = use_SB;
        mog.pruneParams.seuil_SB = seuil_SB;
        mog.pruneParams.zmin = zmin;
        mog.pruneParams.zmax = zmax;
        mog.pruneParams.thetaMin = thetaMin*180/pi;
        mog.pruneParams.thetaMax = thetaMax*180/pi;
        mog.pruneParams.skip_zeros = skip_zeros;

    end
    function zoom(varargin)
        button_state = hzoom.Value;
        if button_state == hzoom.Max
            % toggle button is pressed
            zoom(haxes,'on')
        elseif button_state == hzoom.Min
            % toggle button is not pressed
            zoom(haxes,'off')
        end
    end
    function rotate3d(varargin)
        button_state = hrotate3d.Value;
        if button_state == hrotate3d.Max
            % toggle button is pressed
            rotate3d(haxes,'on')
        elseif button_state == hrotate3d.Min
            % toggle button is not pressed
            rotate3d(haxes,'off')
        end
    end

    function SNR = computeSNR()
        SNR = ones(size(1,mog.data.ntrace));
        [~,i] = max(abs(detrend_rad(mog.data.rdata)));
        
        largeur = 60;
        
        i1 = i-largeur/2;
        i2 = i+largeur/2;
        i1(i1<1) = 1;
        i2(i2>mog.data.nptsptrc) = mog.data.nptsptrc;
        for n=1:mog.data.ntrace
            SNR(n) = std(mog.data.rdata(i1(n):i2(n), n))/std(mog.data.rdata(1:largeur, n));
            if isnan(SNR(n))
                SNR(n) = 1e-9;
            end
        end
    end

    function z=zoomCData()
        z=nan(15,16,3);
        z(:,:,1) = [NaN   NaN   NaN   NaN   NaN   NaN   NaN   NaN     0     0     0     0   NaN   NaN   NaN   NaN
            NaN   NaN   NaN   NaN   NaN   NaN     0     0   NaN   NaN   NaN   NaN     0     0   NaN   NaN
            NaN   NaN   NaN   NaN   NaN     0   NaN   NaN   NaN   NaN   NaN   NaN   NaN   NaN     0   NaN
            NaN   NaN   NaN   NaN   NaN     0   NaN   NaN   NaN     0     0   NaN   NaN   NaN     0   NaN
            NaN   NaN   NaN   NaN     0   NaN   NaN   NaN   NaN     0     0   NaN   NaN   NaN   NaN     0
            NaN   NaN   NaN   NaN     0   NaN   NaN     0     0     0     0     0     0   NaN   NaN     0
            NaN   NaN   NaN   NaN     0   NaN   NaN     0     0     0     0     0     0   NaN   NaN     0
            NaN   NaN   NaN   NaN     0   NaN   NaN   NaN   NaN     0     0   NaN   NaN   NaN   NaN     0
            NaN   NaN   NaN   NaN   NaN     0   NaN   NaN   NaN     0     0   NaN   NaN   NaN     0   NaN
            NaN   NaN   NaN   NaN   NaN     0   NaN   NaN   NaN   NaN   NaN   NaN   NaN   NaN     0   NaN
            NaN   NaN   NaN   NaN     0     0     0     0   NaN   NaN   NaN   NaN     0     0   NaN   NaN
            NaN   NaN   NaN     0     0     0   NaN   NaN     0     0     0     0   NaN   NaN   NaN   NaN
            NaN   NaN     0     0     0   NaN   NaN   NaN   NaN   NaN   NaN   NaN   NaN   NaN   NaN   NaN
            NaN     0     0     0   NaN   NaN   NaN   NaN   NaN   NaN   NaN   NaN   NaN   NaN   NaN   NaN
            NaN   NaN     0   NaN   NaN   NaN   NaN   NaN   NaN   NaN   NaN   NaN   NaN   NaN   NaN   NaN];
        
        z(:,:,2) = [NaN   NaN   NaN   NaN   NaN   NaN   NaN   NaN     0     0     0     0   NaN   NaN   NaN   NaN
            NaN   NaN   NaN   NaN   NaN   NaN     0     0   NaN   NaN   NaN   NaN     0     0   NaN   NaN
            NaN   NaN   NaN   NaN   NaN     0   NaN   NaN   NaN   NaN   NaN   NaN   NaN   NaN     0   NaN
            NaN   NaN   NaN   NaN   NaN     0   NaN   NaN   NaN     0     0   NaN   NaN   NaN     0   NaN
            NaN   NaN   NaN   NaN     0   NaN   NaN   NaN   NaN     0     0   NaN   NaN   NaN   NaN     0
            NaN   NaN   NaN   NaN     0   NaN   NaN     0     0     0     0     0     0   NaN   NaN     0
            NaN   NaN   NaN   NaN     0   NaN   NaN     0     0     0     0     0     0   NaN   NaN     0
            NaN   NaN   NaN   NaN     0   NaN   NaN   NaN   NaN     0     0   NaN   NaN   NaN   NaN     0
            NaN   NaN   NaN   NaN   NaN     0   NaN   NaN   NaN     0     0   NaN   NaN   NaN     0   NaN
            NaN   NaN   NaN   NaN   NaN     0   NaN   NaN   NaN   NaN   NaN   NaN   NaN   NaN     0   NaN
            NaN   NaN   NaN   NaN     0     0     0     0   NaN   NaN   NaN   NaN     0     0   NaN   NaN
            NaN   NaN   NaN     0     0     0   NaN   NaN     0     0     0     0   NaN   NaN   NaN   NaN
            NaN   NaN     0     0     0   NaN   NaN   NaN   NaN   NaN   NaN   NaN   NaN   NaN   NaN   NaN
            NaN     0     0     0   NaN   NaN   NaN   NaN   NaN   NaN   NaN   NaN   NaN   NaN   NaN   NaN
            NaN   NaN     0   NaN   NaN   NaN   NaN   NaN   NaN   NaN   NaN   NaN   NaN   NaN   NaN   NaN];
        
        z(:,:,3) = [NaN   NaN   NaN   NaN   NaN   NaN   NaN   NaN     0     0     0     0   NaN   NaN   NaN   NaN
            NaN   NaN   NaN   NaN   NaN   NaN     0     0   NaN   NaN   NaN   NaN     0     0   NaN   NaN
            NaN   NaN   NaN   NaN   NaN     0   NaN   NaN   NaN   NaN   NaN   NaN   NaN   NaN     0   NaN
            NaN   NaN   NaN   NaN   NaN     0   NaN   NaN   NaN     0     0   NaN   NaN   NaN     0   NaN
            NaN   NaN   NaN   NaN     0   NaN   NaN   NaN   NaN     0     0   NaN   NaN   NaN   NaN     0
            NaN   NaN   NaN   NaN     0   NaN   NaN     0     0     0     0     0     0   NaN   NaN     0
            NaN   NaN   NaN   NaN     0   NaN   NaN     0     0     0     0     0     0   NaN   NaN     0
            NaN   NaN   NaN   NaN     0   NaN   NaN   NaN   NaN     0     0   NaN   NaN   NaN   NaN     0
            NaN   NaN   NaN   NaN   NaN     0   NaN   NaN   NaN     0     0   NaN   NaN   NaN     0   NaN
            NaN   NaN   NaN   NaN   NaN     0   NaN   NaN   NaN   NaN   NaN   NaN   NaN   NaN     0   NaN
            NaN   NaN   NaN   NaN     0     0     0     0   NaN   NaN   NaN   NaN     0     0   NaN   NaN
            NaN   NaN   NaN     0     0     0   NaN   NaN     0     0     0     0   NaN   NaN   NaN   NaN
            NaN   NaN     0     0     0   NaN   NaN   NaN   NaN   NaN   NaN   NaN   NaN   NaN   NaN   NaN
            NaN     0     0     0   NaN   NaN   NaN   NaN   NaN   NaN   NaN   NaN   NaN   NaN   NaN   NaN
            NaN   NaN     0   NaN   NaN   NaN   NaN   NaN   NaN   NaN   NaN   NaN   NaN   NaN   NaN   NaN];
    end

    function r = rotateCData()
        r=nan(15,16,3);
        r(:,:,1) = [NaN   NaN   NaN   NaN   NaN   NaN   NaN   NaN   NaN   NaN   NaN   NaN   NaN   NaN   NaN   NaN
            NaN   NaN   NaN   NaN   NaN     0   NaN     0     0     0     0   NaN   NaN   NaN   NaN   NaN
            NaN   NaN   NaN   NaN     0     0     0   NaN   NaN   NaN   NaN     0     0   NaN   NaN   NaN
            NaN   NaN   NaN     0     0     0     0   NaN   NaN   NaN   NaN   NaN   NaN     0   NaN   NaN
            NaN   NaN   NaN   NaN   NaN   NaN   NaN   NaN   NaN   NaN   NaN   NaN   NaN     0   NaN   NaN
            NaN   NaN   NaN   NaN   NaN   NaN   NaN   NaN   NaN   NaN   NaN   NaN   NaN   NaN   NaN   NaN
            NaN   NaN   NaN   NaN   NaN   NaN   NaN   NaN   NaN   NaN   NaN   NaN   NaN   NaN     0   NaN
            NaN   NaN     0   NaN   NaN   NaN   NaN   NaN   NaN   NaN   NaN   NaN   NaN   NaN   NaN   NaN
            NaN   NaN   NaN   NaN   NaN   NaN   NaN   NaN   NaN   NaN   NaN   NaN   NaN   NaN     0   NaN
            NaN   NaN   NaN   NaN   NaN   NaN   NaN   NaN   NaN   NaN   NaN   NaN   NaN   NaN   NaN   NaN
            NaN   NaN   NaN     0   NaN   NaN   NaN   NaN   NaN   NaN   NaN   NaN   NaN     0   NaN   NaN
            NaN   NaN   NaN   NaN   NaN   NaN   NaN   NaN   NaN   NaN   NaN   NaN   NaN   NaN   NaN   NaN
            NaN   NaN   NaN   NaN   NaN     0   NaN   NaN   NaN   NaN   NaN     0   NaN   NaN   NaN   NaN
            NaN   NaN   NaN   NaN   NaN   NaN   NaN   NaN     0   NaN   NaN   NaN   NaN   NaN   NaN   NaN
            NaN   NaN   NaN   NaN   NaN   NaN   NaN   NaN   NaN   NaN   NaN   NaN   NaN   NaN   NaN   NaN];
        
        r(:,:,2) = [NaN   NaN   NaN   NaN   NaN   NaN   NaN   NaN   NaN   NaN   NaN   NaN   NaN   NaN   NaN   NaN
            NaN   NaN   NaN   NaN   NaN     0   NaN     0     0     0     0   NaN   NaN   NaN   NaN   NaN
            NaN   NaN   NaN   NaN     0     0     0   NaN   NaN   NaN   NaN     0     0   NaN   NaN   NaN
            NaN   NaN   NaN     0     0     0     0   NaN   NaN   NaN   NaN   NaN   NaN     0   NaN   NaN
            NaN   NaN   NaN   NaN   NaN   NaN   NaN   NaN   NaN   NaN   NaN   NaN   NaN     0   NaN   NaN
            NaN   NaN   NaN   NaN   NaN   NaN   NaN   NaN   NaN   NaN   NaN   NaN   NaN   NaN   NaN   NaN
            NaN   NaN   NaN   NaN   NaN   NaN   NaN   NaN   NaN   NaN   NaN   NaN   NaN   NaN     0   NaN
            NaN   NaN     0   NaN   NaN   NaN   NaN   NaN   NaN   NaN   NaN   NaN   NaN   NaN   NaN   NaN
            NaN   NaN   NaN   NaN   NaN   NaN   NaN   NaN   NaN   NaN   NaN   NaN   NaN   NaN     0   NaN
            NaN   NaN   NaN   NaN   NaN   NaN   NaN   NaN   NaN   NaN   NaN   NaN   NaN   NaN   NaN   NaN
            NaN   NaN   NaN     0   NaN   NaN   NaN   NaN   NaN   NaN   NaN   NaN   NaN     0   NaN   NaN
            NaN   NaN   NaN   NaN   NaN   NaN   NaN   NaN   NaN   NaN   NaN   NaN   NaN   NaN   NaN   NaN
            NaN   NaN   NaN   NaN   NaN     0   NaN   NaN   NaN   NaN   NaN     0   NaN   NaN   NaN   NaN
            NaN   NaN   NaN   NaN   NaN   NaN   NaN   NaN     0   NaN   NaN   NaN   NaN   NaN   NaN   NaN
            NaN   NaN   NaN   NaN   NaN   NaN   NaN   NaN   NaN   NaN   NaN   NaN   NaN   NaN   NaN   NaN];
        
        r(:,:,3) = [NaN   NaN   NaN   NaN   NaN   NaN   NaN   NaN   NaN   NaN   NaN   NaN   NaN   NaN   NaN   NaN
            NaN   NaN   NaN   NaN   NaN     0   NaN     0     0     0     0   NaN   NaN   NaN   NaN   NaN
            NaN   NaN   NaN   NaN     0     0     0   NaN   NaN   NaN   NaN     0     0   NaN   NaN   NaN
            NaN   NaN   NaN     0     0     0     0   NaN   NaN   NaN   NaN   NaN   NaN     0   NaN   NaN
            NaN   NaN   NaN   NaN   NaN   NaN   NaN   NaN   NaN   NaN   NaN   NaN   NaN     0   NaN   NaN
            NaN   NaN   NaN   NaN   NaN   NaN   NaN   NaN   NaN   NaN   NaN   NaN   NaN   NaN   NaN   NaN
            NaN   NaN   NaN   NaN   NaN   NaN   NaN   NaN   NaN   NaN   NaN   NaN   NaN   NaN     0   NaN
            NaN   NaN     0   NaN   NaN   NaN   NaN   NaN   NaN   NaN   NaN   NaN   NaN   NaN   NaN   NaN
            NaN   NaN   NaN   NaN   NaN   NaN   NaN   NaN   NaN   NaN   NaN   NaN   NaN   NaN     0   NaN
            NaN   NaN   NaN   NaN   NaN   NaN   NaN   NaN   NaN   NaN   NaN   NaN   NaN   NaN   NaN   NaN
            NaN   NaN   NaN     0   NaN   NaN   NaN   NaN   NaN   NaN   NaN   NaN   NaN     0   NaN   NaN
            NaN   NaN   NaN   NaN   NaN   NaN   NaN   NaN   NaN   NaN   NaN   NaN   NaN   NaN   NaN   NaN
            NaN   NaN   NaN   NaN   NaN     0   NaN   NaN   NaN   NaN   NaN     0   NaN   NaN   NaN   NaN
            NaN   NaN   NaN   NaN   NaN   NaN   NaN   NaN     0   NaN   NaN   NaN   NaN   NaN   NaN   NaN
            NaN   NaN   NaN   NaN   NaN   NaN   NaN   NaN   NaN   NaN   NaN   NaN   NaN   NaN   NaN   NaN];
        
    end
end
