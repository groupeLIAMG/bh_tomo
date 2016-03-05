classdef MogUI < handle
    %MOGUI User interface to manage MOGs

    properties
        mogs
    end
    properties (Dependent)
        Position
        FontSize   % should be within [10 11 12]
    end
    properties (Access=private)
        data_rep
        handles
    end
    properties (SetAccess=immutable,Hidden=true)
        bhUI
    end
    events
        mogAdded, mogDeleted
    end
    
    methods
        function obj = MogUI(b,varargin)
            if isa(b,'BoreholeUI')
                obj.bhUI = b;
            else
                error('Give BoreholeUI object as input')
            end
            obj.mogs = Mog.empty;
            obj.handles.hp = uipanel(varargin{:},...
                'Title','MOGs',...
                'Visible','off',...
                'SizeChangedFcn', @obj.resizeUI);
            obj.data_rep = '';
            obj.addComponents();
            obj.resizeUI();
        end
        function set.mogs(obj, m)
            if isa(m, 'Mog')
                obj.mogs = m;
            else
                error('Mogs should be objects of class Mog')
            end
            obj.updateList(1)
            obj.updateUIfields()
        end
        function set.Position(obj, p)
            obj.handles.hp.Position = p;
        end
        function set.FontSize(obj, s)
            obj.handles.hp.FontSize = s+1;
            obj.handles.addMog.FontSize = s;
            obj.handles.removeMog.FontSize = s;
            obj.handles.listMogs.FontSize = s;
            obj.handles.textType.FontSize = s;
            obj.handles.textTx.FontSize = s;
            obj.handles.textRx.FontSize = s;
            obj.handles.popupType.FontSize = s;
            obj.handles.popupTx.FontSize = s;
            obj.handles.popupRx.FontSize = s;
            obj.handles.airBefore.FontSize = s;
            obj.handles.airAfter.FontSize = s;
            obj.handles.renameMog.FontSize = s;
            obj.handles.importMog.FontSize = s;
            obj.handles.spectra.FontSize = s;
            obj.handles.rawData.FontSize = s;
            obj.handles.zop.FontSize = s;
            obj.handles.statsTt.FontSize = s;
            obj.handles.statsAmp.FontSize = s;
            obj.handles.rays.FontSize = s;
            obj.handles.exportTt.FontSize = s;
            obj.handles.exportTau.FontSize = s;
            obj.handles.prune.FontSize = s;
            obj.handles.textDate.FontSize = s;
            obj.handles.editDate.FontSize = s;
            obj.handles.textFreq.FontSize = s;
            obj.handles.editFreq.FontSize = s;
            obj.handles.textFeedRx.FontSize = s;
            obj.handles.editFeedRx.FontSize = s;
            obj.handles.textFeedTx.FontSize = s;
            obj.handles.editFeedTx.FontSize = s;
            obj.handles.checkDtCorr.FontSize = s;
            obj.handles.editDtCorr.FontSize = s;
            obj.handles.textMultFac.FontSize = s;
            obj.handles.editMultFac.FontSize = s;
            obj.resizeUI();
        end
        function updateBHlist(obj,varargin)
            boreholes = obj.bhUI.boreholes;
            if ~isempty(boreholes)
                names = cell(1,length(boreholes));
                for n=1:length(boreholes)
                    names{n} = boreholes(n).name;
                end
                obj.handles.popupTx.String = names;
                obj.handles.popupRx.String = names;
            else
                obj.handles.popupTx.String = ' ';
                obj.handles.popupRx.String = ' ';
            end
            obj.handles.popupTx.Value = 1;
            obj.handles.popupRx.Value = 1;
        end
    end
    
    methods (Access=private)
        addComponents(obj)
        resizeUI(obj,varargin)
        addMog(obj,varargin)
        updateCoords(obj,varargin)
        [x, y, z, c] = projectBorehole(fdata, prof, nom)
        
        function removeMog(obj,varargin)
            no = obj.handles.listMogs.Value;
            if no>0 && no<=length(obj.mogs)
                ind=1:length(obj.mogs);
                ind = ind~=no;
                obj.mogs = obj.mogs(ind);
                
                obj.updateList()
                obj.updateUIfields()
                obj.notify('mogDeleted')
            end
        end
        function listMogs(obj,varargin)
            obj.updateUIfields()
        end
        function airBefore(obj,varargin)
        end
        function airAfter(obj,varargin)
        end
        function popupType(obj,varargin)
            no = obj.handles.listMogs.Value;
            if no>0 && no<=length(obj.mogs)
                obj.mogs(no).type = obj.handles.popupType.type;
                obj.updateCoords()
            end
        end
        function popupTx(obj,varargin)
            no = obj.handles.listMogs.Value;
            if no>0 && no<=length(obj.mogs)
                obj.mogs(no).Tx = obj.handles.popupTx.Value;
            end
        end
        function popupRx(obj,varargin)
            no = obj.handles.listMogs.Value;
            if no>0 && no<=length(obj.mogs)
                obj.mogs(no).Rx = obj.handles.popupRx.Value;
            end
        end
        function renameMog(obj,varargin)
            no = obj.handles.listMogs.Value;
            if no>0 && no<=length(obj.mogs)
                name = myinputdlg('Enter new name');
                if ~isempty(name)
                    obj.mogs(no).name = char(name);
                    obj.updateList()
                end
            end
        end
        function importMog(obj,varargin)
        end
        function spectra(obj,varargin)
        end
        function rawData(obj,varargin)
        end
        function zop(obj,varargin)
        end
        function statsTt(obj,varargin)
        end
        function statsAmp(obj,varargin)
        end
        function rays(obj,varargin)
        end
        function exportTt(obj,varargin)
        end
        function exportTau(obj,varargin)
        end
        function prune(obj,varargin)
        end
        function editDate(obj,varargin)
            no = obj.handles.listMogs.Value;
            if no>0 && no<=length(obj.mogs)
                obj.mogs(no).data = obj.handles.editDate.String;
            end
        end
        function editFreq(obj,varargin)
            no = obj.handles.listMogs.Value;
            if no>0 && no<=length(obj.mogs)
                obj.mogs(no).data.rnomfreq = str2double(obj.handles.editFreq.String);
            end
        end
        function editFeedRx(obj,varargin)
            no = obj.handles.listMogs.Value;
            if no>0 && no<=length(obj.mogs)
                obj.mogs(no).data.RxOffset = str2double(obj.handles.editFeedRx.String);
            end
        end
        function editFeedTx(obj,varargin)
            no = obj.handles.listMogs.Value;
            if no>0 && no<=length(obj.mogs)
                obj.mogs(no).data.TxOffset = str2double(obj.handles.editFeedTx.String);
            end
        end
        function checkDtCorr(obj,varargin)
            no = obj.handles.listMogs.Value;
            if no>0 && no<=length(obj.mogs)
                obj.mogs(no).user_fac_dt = obj.handles.checkDtCorr.Value;
            end
        end
        function editDtCorr(obj,varargin)
            no = obj.handles.listMogs.Value;
            if no>0 && no<=length(obj.mogs)
                obj.mogs(no).fac_dt = str2double(obj.handles.editDtCorr.String);
            end
        end
        function editMultFac(obj,varargin)
            no = obj.handles.listMogs.Value;
            if no>0 && no<=length(obj.mogs)
                obj.mogs(no).f_et = str2double(obj.handles.editMultFac.String);
            end
        end
        function updateList(obj,varargin)
            value = [];
            if nargin==2
                value = varargin{1};
            end
            if ~isempty(obj.mogs)
                names = cell(1,length(obj.mogs));
                for n=1:length(obj.mogs)
                    names{n} = obj.mogs(n).name;
                end
                obj.handles.listMogs.String = names;
                obj.handles.listMogs.Value = length(obj.mogs);
            else
                obj.handles.listMogs.String = '';
                obj.handles.listMogs.Value = 1;
            end
            if ~isempty(value)
                obj.handles.listMogs.Value = value;
            end
        end
        function updateUIfields(obj)
            no = obj.handles.listMogs.Value;
            if no>0 && no<=length(obj.mogs)
                obj.handles.editDate.String = obj.mogs(no).date;
                obj.handles.editMultFac.String = num2str( obj.mogs(no).f_et );
                obj.handles.editDtCorr.String = num2str( obj.mogs(no).fac_dt );
                obj.handles.editFeedTx.String = num2str( obj.mogs(no).data.TxOffset );
                obj.handles.editFeedRx.String = num2str( obj.mogs(no).data.RxOffset );
                obj.handles.editFreq.String = num2str( obj.mogs(no).data.rnomfreq );
                obj.handles.checkDtCorr.Value = obj.mogs(no).user_fac_dt;
                obj.handles.popupType.Value = obj.mogs(no).type;
                obj.handles.popupRx.Value = obj.mogs(no).Rx;
                obj.handles.popupTx.Value = obj.mogs(no).Tx;
            
            % TODO add air shot 
            end
        end
    end
end

