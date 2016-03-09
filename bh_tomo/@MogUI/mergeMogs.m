function mergeMogs(obj,varargin)
% MERGEMOGS - Combine multiple MOGs into a single
%   MOGs should have the same Tx & Rx boreholes, ais shots data (if
%   present) and antenna offsets (if defined)
%
if isempty(obj.mogs)
    warndlg('No MOG in database')
    return
end
% if numel(obj.mogs)==1
%     warndlg('Only 1 MOG in database')
%     return
% end
fs = 11;
if isa(varargin{1},'matlab.ui.control.UIControl')
    fs = varargin{1}.FontSize;
elseif ispc
    fs = 9;
end

width = 450;
height = 300;

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
    'Fontsize',fs,...
    'RowName',[]);

hkeep = uicontrol('Style','checkbox',...
    'Units','points',...
    'String','Erase MOGs after merge',...
    'Fontsize',fs,...
    'Parent',f);

hmerge = uicontrol('Style','pushbutton',...
    'Units','points',...
    'String','Merge',...
    'Fontsize',fs,...
    'Callback',@doMerge,...
    'Parent',f);

hdone = uicontrol('Style','pushbutton',...
    'Units','points',...
    'String','Done',...
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
        vSize = 22*vFac;
        vBorderTop = 50*vFac;
        vBorder = 20*vFac;

        href.Position = [2*hSpace height-vBorderTop hSize vSize];
        hrefPopup.Position = [2*hSpace height-vBorderTop-vSize hSize vSize];
        hkeep.Position = [2*hSpace height-vBorderTop-3*vSize hSize vSize];
        hcompat.Position = [3*hSpace+hSize height-vBorderTop hSize vSize];
        hcompatTable.Position = [3*hSpace+hSize height-vBorderTop-8*vSize hSize 8*vSize];
        
        hmerge.Position = [2*hSpace vBorder hSize vSize];
        
        hdone.Position = [3*hSpace+hSize vBorder hSize vSize];
        
        
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
        % TODO
        
        
    end
    function getCompat(varargin)
        no = hrefPopup.Value;
        refMog = obj.mogs(no);
        
        compat = cell(numel(obj.mogs)-1,2);
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
            elseif isempty(refMog.p) && ~isempty(obj.mogs(nm).ap)
                test3 = false;
            else
                if refMog.ap == obj.mogs(nm).ap
                    test3 = true;
                end
            end
            
            test4 = refMog.data.TxOffset==obj.mogs(nm).data.TxOffset && ...
                refMog.data.RxOffset==obj.mogs(nm).data.RxOffset;
            if test1 && test2 && test3 && test4
                nc = nc+1;
                compat{nc,1} = true;
                compat{nc,2} = obj.mogs(nm).name;
                ids(nc) = obj.mogs(nm).ID;
            end
        end
        if nc==0
            warndlg('No compatible MOG found')
            hcompatTable.Data = {};
        else
            hcompatTable.Data = compat{1:nc,:};
        end
    end

end