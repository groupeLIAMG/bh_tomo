function mergeMogs(obj,varargin)
% MERGEMOGS - Combine multiple MOGs into a single
%   MOGs should have the same Tx & Rx boreholes, ais shots data (if
%   present) and antenna offsets (if defined)
%
if isempty(obj.mogs)
    warndlg('No MOG in database')
    return
end
if numel(obj.mogs)==1
    warndlg('Only 1 MOG in database')
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
width = 450*vScale;
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
    'Tag','fig_bh_tomo2_prune',...
    'Name','Merge MOGs',...
    'NumberTitle','off',...
    'Resize','off',...
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
    'String','Reference MOG',...
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
    'String','Compatible MOGs',...
    'Fontsize',fs,...
    'HorizontalAlignment','center',...
    'Parent',f);

hcompatTable = uitable(f,...
    'ColumnName',{'','Name'},...
    'ColumnFormat',{'logical','char'},...
    'ColumnWidth',{25 'auto'},...
    'ColumnEditable',[true false],...
    'Units','points',...
    'Fontsize',fs,...
    'RowName',[]);

hkeep = uicontrol('Style','checkbox',...
    'Units','points',...
    'String','Erase MOGs after merge',...
    'Fontsize',fs,...
    'Parent',f);

hname = uicontrol('Style','text',...
    'Units','points',...
    'String','New MOG Name',...
    'Fontsize',fs,...
    'HorizontalAlignment','center',...
    'Parent',f);
hnameEdit= uicontrol('Style','edit',...
    'Units','points',...
    'String','',...
    'Fontsize',fs,...
    'Parent',f);

hmerge = uicontrol('Style','pushbutton',...
    'Units','points',...
    'String','Merge',...
    'Fontsize',fs,...
    'Callback',@doMerge,...
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
        
        hSize = width/3;
        hSpace = hSize/5;
        
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
        hkeep.Position = [2*hSpace height-vBorderTop-3*vSize hSize vSize];
        hname.Position = [2*hSpace height-vBorderTop-6*vSize hSize vSize];
        hnameEdit.Position = [2*hSpace height-vBorderTop-7*vSize hSize vSize];
        
        hcompat.Position = [3*hSpace+hSize height-vBorderTop hSize vSize];
        hcompatTable.Position = [3*hSpace+hSize height-vBorderTop-8*vSize hSize 8*vSize];
        
        hcancel.Position = [2*hSpace vBorder hSize vSize];
        
        hmerge.Position = [3*hSpace+hSize vBorder hSize vSize];
        
        f.Visible = 'on';
    end

    function closeWindow(varargin)
        delete(f)
    end

    function doMerge(varargin)
        
        % check if we actually have MOGs to merge
        tdata = hcompatTable.Data;
        ind = false(size(tdata,1),1);
        if isempty(tdata)
            warndlg('No MOGs compatible for merging')
            return
        else
            for nm=1:size(tdata,1)
                ind(nm) = tdata{nm,1};
            end
        end
        if ~any(ind)
            warndlg('No MOGs selected for merging')
            return
        end
        if isempty(hnameEdit.String)
            warndlg('Enter a name for new MOG')
            return
        end
        
        erase = false;
        if hkeep.Value==1
            choice = questdlg('MOGs will be erased after merge, proceed anyway?',...
                'bh_tomo_db',...
                'Don''t Erase','Cancel','Erase','Erase');
            switch choice
                case 'Don''t Erase'
                case 'Cancel'
                    return
                case 'Erase'
                    erase = true;
            end
        end
        
        % MERGE
        newMog = Mog(hnameEdit.String); % create MOG with new ID
        
        no = hrefPopup.Value;
        refMog = obj.mogs(no);
        newMog.data = refMog.data;
        newMog.av = refMog.av;
        newMog.ap = refMog.ap;
        newMog.Tx = refMog.Tx;
        newMog.Rx = refMog.Rx;
        newMog.f_et = refMog.f_et;
        newMog.type = refMog.type;
        newMog.fac_dt = refMog.fac_dt;
        newMog.user_fac_dt = refMog.user_fac_dt;
        
        newMog.fw = refMog.fw;
        newMog.tt = refMog.tt;
        newMog.et = refMog.et;
        newMog.tt_done = refMog.tt_done;
        newMog.ttTx = refMog.ttTx;
        newMog.ttTx_done = refMog.ttTx_done;
        newMog.amp_tmin = refMog.amp_tmin;
        newMog.amp_tmax = refMog.amp_tmax;
        newMog.amp_done = refMog.amp_done;
        newMog.App = refMog.App;
        newMog.fcentroid = refMog.fcentroid;
        newMog.scentroid = refMog.scentroid;
        newMog.tauApp = refMog.tauApp;
        newMog.tauApp_et = refMog.tauApp_et;
        newMog.tauFce = refMog.tauFce;
        newMog.tauFce_et = refMog.tauFce_et;
        newMog.tauHyb = refMog.tauHyb;
        newMog.tauHyb_et = refMog.tauHyb_et;
        newMog.Tx_z_orig = refMog.Tx_z_orig;
        newMog.Rx_z_orig = refMog.Rx_z_orig;
        newMog.in = refMog.in;
        
        try
            for nc=1:size(tdata,1)
                if tdata{nc,1}==true
                    % look for MOG with good ID
                    for no=1:length(obj.mogs)
                        if ids(nc) == obj.mogs(no).ID
                            
                            mog = obj.mogs(no);
                            
                            newMog.data.rdata = [newMog.data.rdata mog.data.rdata];
                            newMog.data.ntrace = newMog.data.ntrace + mog.data.ntrace;
                            newMog.data.Tx_x = [newMog.data.Tx_x mog.data.Tx_x];
                            newMog.data.Tx_y = [newMog.data.Tx_y mog.data.Tx_y];
                            newMog.data.Tx_z = [newMog.data.Tx_z mog.data.Tx_z];
                            newMog.data.Rx_x = [newMog.data.Rx_x mog.data.Rx_x];
                            newMog.data.Rx_y = [newMog.data.Rx_y mog.data.Rx_y];
                            newMog.data.Rx_z = [newMog.data.Rx_z mog.data.Rx_z];
                            newMog.data.tdata = [newMog.data.tdata mog.data.tdata];
                            
                            newMog.fw = [newMog.fw mog.fw];
                            newMog.tt = [newMog.tt mog.tt];
                            newMog.et = [newMog.et mog.et];
                            newMog.tt_done = [newMog.tt_done mog.tt_done];
                            newMog.ttTx = [newMog.ttTx mog.ttTx];
                            newMog.ttTx_done = [newMog.ttTx_done mog.ttTx_done];
                            newMog.amp_tmin = [newMog.amp_tmin mog.amp_tmin];
                            newMog.amp_tmax = [newMog.amp_tmax mog.amp_tmax];
                            newMog.amp_done = [newMog.amp_done mog.amp_done];
                            newMog.App = [newMog.App mog.App];
                            newMog.fcentroid = [newMog.fcentroid mog.fcentroid];
                            newMog.scentroid = [newMog.scentroid mog.scentroid];
                            newMog.tauApp = [newMog.tauApp mog.tauApp];
                            newMog.tauApp_et = [newMog.tauApp_et mog.tauApp_et];
                            newMog.tauFce = [newMog.tauFce mog.tauFce];
                            newMog.tauFce_et = [newMog.tauFce_et mog.tauFce_et];
                            newMog.tauHyb = [newMog.tauHyb mog.tauHyb];
                            newMog.tauHyb_et = [newMog.tauHyb_et mog.tauHyb_et];
                            newMog.Tx_z_orig = [newMog.Tx_z_orig mog.Tx_z_orig];
                            newMog.Rx_z_orig = [newMog.Rx_z_orig mog.Rx_z_orig];
                            newMog.in = [newMog.in mog.in];
                            
                        end
                    end
                end
            end
        catch ME
            warndlg({'Could not merge MOGs',ME.message})
            return
        end
        
        if erase==true
            keep = true(size(obj.mogs));
            for no=1:length(obj.mogs)
                if refMog.ID == obj.mogs(no).ID
                    keep(no) = false;
                end
            end
            for nc=1:size(tdata,1)
                if tdata{nc,1}==true
                    % look for MOG with good ID
                    for no=1:length(obj.mogs)
                        if ids(nc) == obj.mogs(no).ID
                            keep(no) = false;
                        end
                    end
                end
            end
            obj.mogs = obj.mogs(keep);
        end
        
        % append new MOG
        obj.mogs(end+1) = newMog;
            
        obj.notify('mogAdded')
        closeWindow()
        
    end
    function getCompat(varargin)
        no = hrefPopup.Value;
        refMog = obj.mogs(no);
        
        ids = [];
        nc = 0;
        for nm=1:numel(obj.mogs)
            if nm==no
                continue
            end
            test1 = ( refMog.Tx==obj.mogs(nm).Tx && refMog.Rx==obj.mogs(nm).Rx );
            test2 = false;
            if isempty(refMog.av) && isempty(obj.mogs(nm).av)
                test2 = true;
            elseif ~isempty(refMog.av) && isempty(obj.mogs(nm).av)
                test2 = false;
            elseif isempty(refMog.av) && ~isempty(obj.mogs(nm).av)
                test2 = false;
            else
                if refMog.av == obj.mogs(nm).av
                    test2 = true;
                end
            end
            test3 = false;
            if isempty(refMog.ap) && isempty(obj.mogs(nm).ap)
                test3 = true;
            elseif ~isempty(refMog.ap) && isempty(obj.mogs(nm).ap)
                test3 = false;
            elseif isempty(refMog.ap) && ~isempty(obj.mogs(nm).ap)
                test3 = false;
            else
                if refMog.ap == obj.mogs(nm).ap
                    test3 = true;
                end
            end
            
            test4 = refMog.data.TxOffset==obj.mogs(nm).data.TxOffset && ...
                refMog.data.RxOffset==obj.mogs(nm).data.RxOffset;
            
            test5 = refMog.type == obj.mogs(nm).type;
            
            if test1 && test2 && test3 && test4 && test5
                nc = nc+1;
                compat{nc,1} = true; %#ok<AGROW>
                compat{nc,2} = obj.mogs(nm).name; %#ok<AGROW>
                ids(nc) = obj.mogs(nm).ID;
            end
        end
        if nc==0
            warndlg('No compatible MOG found')
            hcompatTable.Data = {};
        else
            hcompatTable.Data = compat;
        end
    end

end