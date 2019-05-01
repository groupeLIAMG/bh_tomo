function createDt(obj,varargin)

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

width = 475;
height = 300*vScale;

% get screen size
su = get(groot,'Units');
set(groot,'Units','points')
scnsize = get(groot,'ScreenSize');
pos = [scnsize(3)/2-width/2 scnsize(4)/2-height/2 width height];
set(groot,'Units',su)       % Restore default root screen units

f = figure('Visible','off',...
    'Units','points',...
    'Position',pos,...
    'Tag','fig_bh_tomo2_createDt',...
    'Name',['Create ',char(916),' MOG'],...
    'NumberTitle','off',...
    'ToolBar','none',...,
    'MenuBar','None',...
    'SizeChangedFcn',@resizeUI,...
    'CloseRequestFcn',@closeWindow);

names_mog = cell(1,length(obj.mogs));
for n=1:length(obj.mogs)
    names_mog{n} = obj.mogs(n).name;
end
ids = [];

href = uicontrol('Style','text',...
    'Units','points',...
    'String','Minuend MOG',...
    'Fontsize',fs,...
    'HorizontalAlignment','center',...
    'Parent',f);
hrefPopup = uicontrol('Style','popupmenu',...
    'Units','points',...
    'String',names_mog,...
    'Fontsize',fs,...
    'Callback',@getCompat,...
    'Parent',f);

hcompat = uicontrol('Style','text',...
    'Units','points',...
    'String','Subtrahend MOG',...
    'Fontsize',fs,...
    'HorizontalAlignment','center',...
    'Parent',f);

hcompatList = uicontrol('Style','popupmenu',...
    'Units','points',...
    'String','-',...
    'Fontsize',fs,...
    'Parent',f);

hname = uicontrol('Style','text',...
    'Units','points',...
    'String','Name of Difference MOG',...
    'Fontsize',fs,...
    'HorizontalAlignment','right',...
    'Parent',f);
hnameEdit= uicontrol('Style','edit',...
    'Units','points',...
    'String','',...
    'Fontsize',fs,...
    'Parent',f);

hoff = uicontrol('Style','text',...
    'Units','points',...
    'String','Offset Tolerance',...
    'Fontsize',fs,...
    'HorizontalAlignment','right',...
    'Parent',f);
hoffEdit= uicontrol('Style','edit',...
    'Units','points',...
    'String','0.05',...
    'Fontsize',fs,...
    'Parent',f);

bgmethod = uibuttongroup('Visible', 'off', ...
    'Units','points',...
    'Title','Method',...
    'Fontsize',fs,...
    'Parent',f);
bxcorr = uicontrol(bgmethod,...
    'Style', 'radiobutton',...
    'String', 'Use X Correlation',...
    'Fontsize',fs,...
    'Units','normalized',...
    'HandleVisibility','off');
bpicks = uicontrol(bgmethod,...
    'Style', 'radiobutton',...
    'String', 'Use Picked Traveltimes',...
    'Fontsize',fs,...
    'Units','normalized',...
    'HandleVisibility','off');
hmaxlag = uicontrol('Style','text',...
    'Units','points',...
    'String',['Max lag (',obj.mogs(1).data.tunits,')'],...
    'Fontsize',fs,...
    'HorizontalAlignment','right',...
    'Units','normalized',...
    'Parent',bgmethod);
hmaxlagEdit= uicontrol('Style','edit',...
    'Units','points',...
    'String','10',...
    'Fontsize',fs,...
    'Units','normalized',...
    'Parent',bgmethod);
bgmethod.Visible = 'on';

hdone = uicontrol('Style','pushbutton',...
    'Units','points',...
    'String','Done',...
    'Fontsize',fs,...
    'Callback',@done,...
    'Parent',f);

hcancel = uicontrol('Style','pushbutton',...
    'Units','points',...
    'String','Cancel',...
    'Fontsize',fs,...
    'Callback',@closeWindow,...
    'Parent',f);


getCompat()

uiwait(f)

    function resizeUI(varargin)
        f.Visible = 'off';
        
        width = f.Position(3);
        height = f.Position(4);
        
        hSize = width/2.5;
        hSpace = width/25;
        
        vFac = 0.8*height/500;
        if vFac<1
            vFac = 1;
        end
        if ispc
            vFac = vFac*0.81;
        end
        vSize = 22*vFac;
        vBorderTop = 50*vFac;
        vBorder = 20*vFac;
        
        href.Position = [2*hSpace height-vBorderTop hSize vSize];
        hrefPopup.Position = [2*hSpace height-vBorderTop-vSize hSize vSize];
        
        hoff.Position = [2*hSpace height-vBorderTop-2.5*vSize hSize vSize];
        hoffEdit.Position = [3*hSpace+hSize height-vBorderTop-2.5*vSize hSize vSize];

        hname.Position = [2*hSpace height-vBorderTop-4*vSize hSize vSize];
        hnameEdit.Position = [3*hSpace+hSize height-vBorderTop-4*vSize hSize vSize];
        
        hcompat.Position = [3*hSpace+hSize height-vBorderTop hSize vSize];
        hcompatList.Position = [3*hSpace+hSize height-vBorderTop-vSize hSize vSize];
        
        bgmethod.Position = [2*hSpace vBorder+1.3*vSize 2*hSize+hSpace 4*vSize];
        bpicks.Position = [0.05 0.55 0.8 0.35];
        bxcorr.Position = [0.05 0.1 0.4 0.35];
        hmaxlag.Position = [0.52 0.05 0.2 0.35];
        hmaxlagEdit.Position = [0.75 0.1 0.2 0.35];

        hcancel.Position = [2*hSpace vBorder hSize vSize];
        
        hdone.Position = [3*hSpace+hSize vBorder hSize vSize];
        
        f.Visible = 'on';

    end

    function closeWindow(varargin)
        delete(f)
    end

    function done(varargin)

        if strcmp(hcompatList.String,'-')==1 
            warndlg('No MOG compatible')
            return
        end
        if isempty(hnameEdit.String)
            warndlg('Enter a name for new MOG')
            return
        end
        
        % Check if traveltimes were picked
        no = hrefPopup.Value;
        refMog = obj.mogs(no);
        ind = refMog.tt~=-1;
        if ~any(ind)
            warndlg(['Traveltimes were not picked for ',refMog.name])
            return
        end
        
        newMog = Mog(hnameEdit.String); % create MOG with new ID
        newMog.data = refMog.data.getSubset(ind);
        newMog.traces = refMog.traces(:,ind);
        newMog.data.comment = 'made by MogUI.createDt';
        newMog.av = [];
        newMog.ap = [];
        newMog.Tx = refMog.Tx;
        newMog.Rx = refMog.Rx;
        newMog.f_et = refMog.f_et;
        newMog.type = refMog.type;
        newMog.fac_dt = 1;
        newMog.user_fac_dt = 0;
        
        newMog.fw = refMog.fw;
        newMog.tt = refMog.tt(ind);
        newMog.et = refMog.et(ind);
        newMog.tt_done = refMog.tt_done(ind);
        newMog.ttTx = refMog.ttTx;
        newMog.ttTx_done = refMog.ttTx_done;
        newMog.amp_tmin = refMog.amp_tmin(ind);
        newMog.amp_tmax = refMog.amp_tmax(ind);
        newMog.amp_done = refMog.amp_done(ind);
        newMog.App = refMog.App(ind);
        newMog.fcentroid = refMog.fcentroid(ind);
        newMog.scentroid = refMog.scentroid(ind);
        newMog.tauApp = refMog.tauApp(ind);
        newMog.tauApp_et = refMog.tauApp_et(ind);
        newMog.tauFce = refMog.tauFce(ind);
        newMog.tauFce_et = refMog.tauFce_et(ind);
        newMog.tauHyb = refMog.tauHyb(ind);
        newMog.tauHyb_et = refMog.tauHyb_et(ind);
        newMog.Tx_z_orig = refMog.Tx_z_orig(ind);
        newMog.Rx_z_orig = refMog.Rx_z_orig(ind);
        newMog.TxCosDir = refMog.TxCosDir(ind,:);
        newMog.RxCosDir = refMog.RxCosDir(ind,:);
        newMog.in = refMog.in(ind);
        
        id = ids(hcompatList.Value);
        for mog=obj.mogs
            if id == mog.ID
                break
            end
        end
        
        maxOffset = str2double(hoffEdit.String);
        
        newMog.tt(:) = -1;
        newMog.tt_done(:) = false;
        nFound = 0;
        
        h = waitbar(0,'1','Name','Looking for Common Tx-Rx Pairs',...
            'CreateCancelBtn',...
            'setappdata(gcbf,''canceling'',1)');
        setappdata(h,'canceling',0)
        
        if bgmethod.SelectedObject == bpicks
            oldTT = refMog.getCorrectedTravelTimes(obj.air);
            oldTT = oldTT(ind);
            mogTT = mog.getCorrectedTravelTimes(obj.air);

            for n1=1:newMog.data.ntrace
                if getappdata(h,'canceling')
                    break
                end
                % Report current estimate in the waitbar's message field
                waitbar(n1/newMog.data.ntrace,h,sprintf('Trace no %d: Found %d Pairs',n1,nFound))
                for n2=1:mog.data.ntrace
                    if sqrt((mog.data.Tx_x(n2)-newMog.data.Tx_x(n1))^2 + ...
                            (mog.data.Tx_y(n2)-newMog.data.Tx_y(n1))^2 + ...
                            (mog.data.Tx_z(n2)-newMog.data.Tx_z(n1))^2) < maxOffset && ...
                            sqrt((mog.data.Rx_x(n2)-newMog.data.Rx_x(n1))^2 + ...
                            (mog.data.Rx_y(n2)-newMog.data.Rx_y(n1))^2 + ...
                            (mog.data.Rx_z(n2)-newMog.data.Rx_z(n1))^2) < maxOffset && ...
                            mog.tt_done(n2)
                        newMog.tt(n1) = oldTT(n1)-mogTT(n2);
                        newMog.tt_done(n1) = true;
                        nFound = nFound+1;
                        break
                    end
                end
            end

        elseif bgmethod.SelectedObject == bxcorr
            
            [t0new,~,~] = refMog.corr_t0(length(refMog.tt), obj.air(refMog.av),...
                                         obj.air(refMog.ap));
            t0new = t0new(ind);

            [t0mog,~,~] = mog.corr_t0(length(mog.tt), obj.air(mog.av),...
                                      obj.air(mog.ap));
            hf = figure();
            maxlag = round(str2double(hmaxlagEdit.String)/mog.data.timec);
            for n1=1:newMog.data.ntrace
                if getappdata(h,'canceling')
                    delete(h)
                    return
                end
                % Report current estimate in the waitbar's message field
                waitbar(n1/newMog.data.ntrace,h,sprintf('Trace no %d: Found %d Pairs',n1,nFound))
                for n2=1:mog.data.ntrace
                    if sqrt((mog.data.Tx_x(n2)-newMog.data.Tx_x(n1))^2 + ...
                            (mog.data.Tx_y(n2)-newMog.data.Tx_y(n1))^2 + ...
                            (mog.data.Tx_z(n2)-newMog.data.Tx_z(n1))^2) < maxOffset && ...
                       sqrt((mog.data.Rx_x(n2)-newMog.data.Rx_x(n1))^2 + ...
                            (mog.data.Rx_y(n2)-newMog.data.Rx_y(n1))^2 + ...
                            (mog.data.Rx_z(n2)-newMog.data.Rx_z(n1))^2) < maxOffset && ...
                            mog.tt_done(n2)
                        
                        trc1 = newMog.traces(:,n1);
                        trc2 = mog.traces(:,n2);
                        ns0 = round(t0new(n1)/mog.data.timec);
                        ns1 = round(t0mog(n2)/mog.data.timec);
                        dns = ns0-ns1;
                        trc2 = circshift(trc2, dns);
                        
                        figure(hf)
                        plot(newMog.data.timestp, trc1)
                        hold on
                        plot(mog.data.timestp, trc2)
                        hold off
                        legend(newMog.name,mog.name, 'Interpreter','none')
                        title(['t_0: tcr2 shifted by ',num2str(dns),' sample(s)'])
                        
                        [c, lags] = xcorr(trc1, trc2, maxlag);
                        [~, ixcm] = max(c);
                        lag = lags(ixcm);
                        newMog.tt(n1) = -lag * mog.data.timec;
                        
                        newMog.tt_done(n1) = true;
                        nFound = nFound+1;  
                        break
                    end
                end
            end

                        
        else
            warndlg('Select Method')
            delete(h)
            return
        end
        
        delete(h)       % DELETE the waitbar; don't try to CLOSE it.
        
        if nFound == 0
            warndlg('No Common Tx-Rx Found')
            return
        end
        
        debut = obj.mogs(end).no_traces(end);
        newMog.no_traces = debut + (1:newMog.data.ntrace);
        newMog.sorted = false;
        % append new MOG
        obj.mogs(end+1) = newMog;
            
        obj.notify('mogAdded')
        closeWindow()
    end
    function getCompat(varargin)
        no = hrefPopup.Value;
        refMog = obj.mogs(no);

        ids = [];
        nc = 1;
        compat{nc} = obj.mogs(no).name; %#ok<AGROW>
        ids(nc) = obj.mogs(no).ID;
        for nm=1:numel(obj.mogs)
            if nm==no
                continue
            end
            test1 = ( refMog.Tx==obj.mogs(nm).Tx && refMog.Rx==obj.mogs(nm).Rx );
            
            test2 = refMog.data.TxOffset==obj.mogs(nm).data.TxOffset && ...
                refMog.data.RxOffset==obj.mogs(nm).data.RxOffset;
            
            test3 = refMog.type == obj.mogs(nm).type;
            
            if test1 && test2 && test3
                nc = nc+1;
                compat{nc} = obj.mogs(nm).name; %#ok<AGROW>
                ids(nc) = obj.mogs(nm).ID;
            end
        end
        if nc==0
            warndlg('No compatible MOG found')
            hcompatList.String = '-';
        else
            hcompatList.String = compat;
        end

    end
end

