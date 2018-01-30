function zop(obj,varargin)

fs = 11;
if isa(varargin{1},'matlab.ui.control.UIControl')
    fs = varargin{1}.FontSize;
elseif ispc
    fs = 9;
end

no = obj.handles.listMogs.Value;
if no<=0 || no>length(obj.mogs)
    warndlg('No MOG selected')
    return
end


width = 1000;
height = 600;

% get screen size
su = get(groot,'Units');
set(groot,'Units','points')
scnsize = get(groot,'ScreenSize');
pos = [scnsize(3)/2-width/2 scnsize(4)/2-height/2 width height];
set(groot,'Units',su)       % Restore default root screen units


f = figure('Visible','off',...
    'Units','points',...
    'Position',pos,...
    'Tag','fig_bh_tomo2_zop',...
    'Name','Zero-Offset Profile',...
    'NumberTitle','off',...
    'ToolBar','none',...,
    'MenuBar','None',...
    'SizeChangedFcn',@resizeUI,...
    'CloseRequestFcn',@closeWindow);

haxes1 = axes('Units','points','Parent',f);
haxes2 = axes('Units','points','Parent',f);

hcontrol = uipanel(f,'Title','Control',...
    'Units','points',...
    'FontSize',fs);

hprint = uicontrol('Style','pushbutton',...
    'String','Print',...
    'Units','points',...
    'FontSize',fs,...
    'Callback',@print,...
    'Parent',hcontrol);

hrays = uicontrol('Style','pushbutton',...
    'String','Show Rays',...
    'Units','points',...
    'FontSize',fs,...
    'Callback',@rays,...
    'Parent',hcontrol);

hspreading = uicontrol('Style','checkbox',...
    'String','Corr. Geometrical Spreading',...
    'Units','points',...
    'FontSize',fs,...
    'Callback',@updateZOP,...
    'Parent',hcontrol);
hamplitude = uicontrol('Style','checkbox',...
    'String','Show Amplitude Data',...
    'Units','points',...
    'FontSize',fs,...
    'Callback',@updateFig2,...
    'Parent',hcontrol);
hcont = uicontrol('Style','checkbox',...
    'String','Show BH Velocity Constraints',...
    'Units','points',...
    'FontSize',fs,...
    'Callback',@updateFig2,...
    'Parent',hcontrol);
hvapp = uicontrol('Style','checkbox',...
    'String','Show Apparent Velocity',...
    'Units','points',...
    'FontSize',fs,...
    'Callback',@updateFig2,...
    'Parent',hcontrol);

hscale = uipanel(hcontrol,'Title','Color Scale',...
    'Units','points',...
    'FontSize',fs);

hscalePopup = uicontrol('Style','popupmenu',...
    'String',{'Low','Medium','High'},...
    'Units','points',...
    'FontSize',fs,...
    'Callback',@setCaxis,...
    'Parent',hscale);

hscaleEdit = uicontrol('Style','edit',...
    'String','',...
    'Units','points',...
    'FontSize',fs,...
    'Callback',@updateZOP,...
    'Parent',hscale);

hoffset = uicontrol('Style','text',...
    'String','Vertical Tx-Rx Offset Tolerance',...
    'Units','points',...
    'FontSize',fs,...
    'HorizontalAlignment','center',...
    'Parent',hcontrol);
hoffsetEdit = uicontrol('Style','edit',...
    'String','',...
    'Units','points',...
    'FontSize',fs,...
    'Callback',@updateZOP,...
    'Parent',hcontrol);

hzmin = uicontrol('Style','text',...
    'String','z min',...
    'Units','points',...
    'FontSize',fs,...
    'HorizontalAlignment','right',...
    'Parent',hcontrol);
hzminEdit = uicontrol('Style','edit',...
    'String','',...
    'Units','points',...
    'FontSize',fs,...
    'Callback',@ZLim1,...
    'Parent',hcontrol);
hzmax = uicontrol('Style','text',...
    'String','z max',...
    'Units','points',...
    'FontSize',fs,...
    'HorizontalAlignment','right',...
    'Parent',hcontrol);
hzmaxEdit = uicontrol('Style','edit',...
    'String','',...
    'Units','points',...
    'FontSize',fs,...
    'Callback',@ZLim2,...
    'Parent',hcontrol);

htmin = uicontrol('Style','text',...
    'String','t min',...
    'Units','points',...
    'FontSize',fs,...
    'HorizontalAlignment','right',...
    'Parent',hcontrol);
htminEdit = uicontrol('Style','edit',...
    'String','',...
    'Units','points',...
    'FontSize',fs,...
    'Callback',@axes1XLim1,...
    'Parent',hcontrol);
htmax = uicontrol('Style','text',...
    'String','t max',...
    'Units','points',...
    'FontSize',fs,...
    'HorizontalAlignment','right',...
    'Parent',hcontrol);
htmaxEdit = uicontrol('Style','edit',...
    'String','',...
    'Units','points',...
    'FontSize',fs,...
    'Callback',@axes1XLim2,...
    'Parent',hcontrol);

mog = obj.mogs(no);
av = [];
ap = [];
if ~isempty(obj.air)
    av = obj.air(mog.av);
    ap = obj.air(mog.ap);
end
if mog.data.rstepsz==0
    hoffsetEdit.String = '0.5';
else
    hoffsetEdit.String = num2str(mog.data.rstepsz*0.5);
end
indZop = [];
t0 = [];
updateZOP();
hzminEdit.String = num2str(haxes1.YLim(1));
hzmaxEdit.String = num2str(haxes1.YLim(2));
bhTx.scont = [];
bhRx.scont = [];
if nargin>=4
    bhTx = varargin{3}.boreholes(mog.Tx);
    bhRx = varargin{3}.boreholes(mog.Rx);
end

f.Visible = 'on';

    function resizeUI(varargin)
        f.Visible = 'off';
        
        width = f.Position(3);
        height = f.Position(4);
        
        axBorder = haxes1.Position(2);
        inWidth = width-2.5*axBorder;
        
        haxes1.Position = [axBorder axBorder 0.4*inWidth height-2*axBorder];
        haxes2.Position = [1.5*axBorder+0.4*inWidth axBorder 0.3*inWidth height-2*axBorder];
        hcontrol.Position = [2*axBorder+0.7*inWidth axBorder 0.3*inWidth height-1.75*axBorder];
        
        inWidth = inWidth-10;
        hBorder = 15;
        width = 0.3*inWidth-2*hBorder;
        height = height-2*axBorder;
        
        vSize = height/(17*0.3+16);
        vSpace = 0.3*vSize;
        
        hprint.Position = [hBorder 2*vSpace width vSize];
        hrays.Position = [hBorder 3*vSpace+vSize width vSize];
        
        hspreading.Position = [hBorder 5*vSpace+2*vSize width vSize];
        hamplitude.Position = [hBorder 5.5*vSpace+3*vSize width vSize];
        hcont.Position = [hBorder 6*vSpace+4*vSize width vSize];
        hvapp.Position = [hBorder 6.5*vSpace+5*vSize width vSize];
        
        hscale.Position = [hBorder 8.5*vSpace+6*vSize width 4*vSize];
        hscalePopup.Position = [hBorder vSpace width-2.3*hBorder vSize];
        hscaleEdit.Position = [width/4 2.5*vSpace+vSize width/2 vSize];
        
        
        hoffsetEdit.Position = [hBorder+width/3 10*vSpace+10*vSize width/3 vSize];
        hoffset.Position = [hBorder 10*vSpace+11*vSize width vSize];
        
        hzmax.Position = [hBorder 13*vSpace+12*vSize width/2 vSize];
        hzmaxEdit.Position = [2*hBorder+width/2 13*vSpace+12*vSize width/3 vSize];
        hzmin.Position = [hBorder 13.5*vSpace+13*vSize width/2 vSize];
        hzminEdit.Position = [2*hBorder+width/2 13.5*vSpace+13*vSize width/3 vSize];
        
        htmax.Position = [hBorder 14.5*vSpace+14*vSize width/2 vSize];
        htmaxEdit.Position = [2*hBorder+width/2 14.5*vSpace+14*vSize width/3 vSize];
        htmin.Position = [hBorder 15*vSpace+15*vSize width/2 vSize];
        htminEdit.Position = [2*hBorder+width/2 15*vSpace+15*vSize width/3 vSize];
        
        f.Visible = 'on';
    end
    function closeWindow(varargin)
        delete(f)
    end
    function updateZOP(varargin)
        
        [indZop, t0] = plotZOP(mog, av, ap, haxes1, true,...
            hspreading.Value, str2double(hoffsetEdit.String));
        clim = caxis(haxes1);
        Amax = clim(2);
        hscaleEdit.String = num2str(Amax);
        xlim = haxes1.XLim;
        if isempty(htminEdit.String)
            htminEdit.String = num2str(xlim(1));
        else
            axes1XLim1()
        end
        if isempty(htmaxEdit.String)
            htmaxEdit.String = num2str(xlim(2));
        else
            axes1XLim2()
        end
        
    end
    function setCaxis(varargin)
        if hscalePopup.Value==1
            caxis(haxes1,[-28000 28000])
            hscaleEdit.String = '28000';
        elseif hscalePopup.Value==2
            caxis(haxes1,[-7000 7000])
            hscaleEdit.String = '7000';
        elseif hscalePopup.Value==3
            caxis(haxes1,[-700 700])
            hscaleEdit.String = '700';
        end
    end
    function rays(varargin)
        xy = sqrt( (mog.data.Tx_x-mog.data.Rx_x).^2 + ...
            (mog.data.Tx_y-mog.data.Rx_y).^2 );
        figure
        plot([zeros(size(xy(indZop)))' xy(indZop)']', ...
            [mog.data.Tx_z(indZop)' mog.data.Rx_z(indZop)']','g')
        set(gca, 'DataAspectRatio',[1 1 1])
        axis tight
        xlabel('Tx-Rx Distance [m]')
        ylabel('Elevation [m]')
        
    end
    function print(varargin)
        hcontrol.Visible = 'off';
        printdlg(f)
        hcontrol.Visible = 'on';
    end
    function axes1XLim1(varargin)
        val = str2double(htminEdit.String);
        if val<haxes1.XLim(2)
            haxes1.XLim(1) = val;
        else
            warndlg('t min should be smaller than t max')
        end
    end
    function axes1XLim2(varargin)
        val = str2double(htmaxEdit.String);
        if val>haxes1.XLim(1)
            haxes1.XLim(2) = val;
        else
            warndlg('t max should be larger than t min')
        end
    end
    function ZLim1(varargin)
        val = str2double(hzminEdit.String);
        if val<haxes1.YLim(2)
            haxes1.YLim(1) = val;
            haxes2.YLim(1) = val;
        else
            warndlg('z min should be smaller than z max')
        end
    end
    function ZLim2(varargin)
        val = str2double(hzmaxEdit.String);
        if val>haxes1.YLim(1)
            haxes1.YLim(2) = val;
            haxes2.YLim(2) = val;
        else
            warndlg('z max should be larger than z min')
        end
    end
    function updateFig2(varargin)
        if hamplitude.Value==1 && hvapp.Value==1
            plotVelAmp()
        elseif hamplitude.Value==1
            plotAmp()
        elseif hvapp.Value==1
            plotAppVel()
        else
            cla(haxes2,'reset')
        end
    end

    function plot_scont(axes,varargin)
        hold(axes,'on')
        if ~isempty( bhTx.scont )
            plot(axes, 1./bhTx.scont.valeur, bhTx.scont.z,'b')
        end
        if ~isempty( bhRx.scont )
            plot(axes, 1./bhRx.scont.valeur, bhRx.scont.z,'g')
        end
        hold(axes,'off')
    end
    function plotAppVel(varargin)
        if isempty(mog)
            return
        end
        haxes2 = haxes2(1);
        cla(haxes2,'reset')
        lrai = sqrt( (mog.data.Tx_x(indZop)-mog.data.Rx_x(indZop)).^2 + ...
            (mog.data.Tx_y(indZop)-mog.data.Rx_y(indZop)).^2 + ...
            (mog.data.Tx_z(indZop)-mog.data.Rx_z(indZop)).^2 );
        tt = mog.tt(indZop)-t0(indZop);
        et = mog.et(indZop);
        vApp = lrai./tt;
        vMin = lrai./(tt+et);
        vMax = lrai./(tt-et);
        Tx_z = mog.data.Tx_z(indZop);
        ind = mog.tt(indZop)~=-1;
        plot(haxes2, vApp(ind), Tx_z(ind), 'o')
        z = Tx_z(ind);
        vMin = vMin(ind);
        vMax = vMax(ind);
        hold(haxes2, 'on')
        for n=1:length(vMin)
            plot(haxes2, [vMin(n) vMax(n)], [z(n) z(n)],'Color',[0.75 0.75 0.75])
        end
        hold(haxes2, 'off')
        if hcont.Value==1
            plot_scont(haxes2)
        end
        grid(haxes2, 'on')
        xlabel(haxes2, ['Velocity [',mog.data.cunits,'/',mog.data.tunits,']'])
        ylabel(haxes2, ['Elevation [',mog.data.cunits,']'])
        title(haxes2, 'Apparent velocity')
        haxes2.YLim = haxes1.YLim;
    end

    function plotAmp(varargin)
        if isempty(mog)
            return
        end
        haxes2 = haxes2(1);
        cla(haxes2,'reset')
        lrai = sqrt( (mog.data.Tx_x(indZop)-mog.data.Rx_x(indZop)).^2 + ...
            (mog.data.Tx_y(indZop)-mog.data.Rx_y(indZop)).^2 + ...
            (mog.data.Tx_z(indZop)-mog.data.Rx_z(indZop)).^2 );
        tmin = mog.amp_tmin(indZop);
        tmax = mog.amp_tmax(indZop);
        imin = zeros(size(tmin));
        imax = zeros(size(tmax));
        App = nan(size(tmin));
        for n=1:length(tmin), imin(n) = findnear(tmin(n), mog.data.timestp); end
        for n=1:length(tmax), imax(n) = findnear(tmax(n), mog.data.timestp); end
        for n=1:length(tmin)
            if imax(n)>imin(n)
                App(n) = lrai(n)*(max(mog.data.rdata(imin:imax,indZop(n))) - ...
                    min(mog.data.rdata(imin:imax,indZop(n))));
            end
        end
        Tx_z = mog.data.Tx_z(indZop);
        axes(haxes2)
        semilogx(App, Tx_z,'LineStyle','none','Marker','o','MarkerEdgeColor',[0 0.5 0])
        xlabel('Amplitude')
        ylabel('Elevation [m]')
        title('Peak-to-peak amplitude')
        if hcont.Value==1
            plot_scont(haxes2)
        end
        haxes2.YLim = haxes1.YLim;
        
    end
    function plotVelAmp(varargin)
        if isempty(mog)
            return
        end
        cla(haxes2,'reset')
        lrai = sqrt( (mog.data.Tx_x(indZop)-mog.data.Rx_x(indZop)).^2 + ...
            (mog.data.Tx_y(indZop)-mog.data.Rx_y(indZop)).^2 + ...
            (mog.data.Tx_z(indZop)-mog.data.Rx_z(indZop)).^2 );
        tt = mog.tt(indZop)-t0(indZop);
        et = mog.et(indZop);
        vApp = lrai./tt;
        vMin = lrai./(tt+et);
        vMax = lrai./(tt-et);
        Tx_z = mog.data.Tx_z(indZop);
        ind = mog.tt(indZop)~=-1;
        tmin = mog.amp_tmin(indZop);
        tmax = mog.amp_tmax(indZop);
        imin = zeros(size(tmin));
        imax = zeros(size(tmax));
        App = nan(size(tmin));
        for n=1:length(tmin), imin(n) = findnear(tmin(n), mog.data.timestp); end
        for n=1:length(tmax), imax(n) = findnear(tmax(n), mog.data.timestp); end
        for n=1:length(tmin)
            if imax(n)>imin(n)
                App(n) = lrai(n)*(max(mog.data.rdata(imin(n):imax(n),indZop(n))) - ...
                    min(mog.data.rdata(imin(n):imax(n),indZop(n))));
            end
        end
        
        axes(haxes2)
        [ax,h1,h2] = plotxx(vApp(ind), Tx_z(ind), App, Tx_z, 'plot', 'semilogx');
        set(get(ax(1),'Xlabel'),'String','Apparent velocity [m/ns]')
        set(get(ax(2),'Xlabel'),'String','Peak-to-peak amplitude')
        set(h1,'LineStyle','none','Marker','o','MarkerEdgeColor','b')
        set(h2,'LineStyle','none','Marker','o','MarkerEdgeColor',[0 0.5 0])
        set(ax(2),'YAxisLocation','right')
        ylabel('Elevation [m]')
        
        z = Tx_z(ind);
        vMin = vMin(ind);
        vMax = vMax(ind);
        hold(ax(1),'on')
        for n=1:length(vMin)
            plot([vMin(n) vMax(n)], [z(n) z(n)],'Color',[0.75 0.75 0.75])
        end
        hold(ax(1),'off')
        if hcont.Value==1
            plot_scont(ax(1))
        end
        ax(1).YLim = haxes1.YLim;
        ax(2).YLim = haxes1.YLim;
        haxes2 = ax;
    end
end