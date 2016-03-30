classdef MogUI < handle
    %MOGUI User interface to manage MOGs
    
    properties
        mogs       % MOGs data
        air        % air shots data
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
        mogAdded, mogDeleted, mogEdited
    end
    
    methods
        function obj = MogUI(b,varargin)
            if isa(b,'BoreholeUI')
                obj.bhUI = b;
            else
                error('Give BoreholeUI object as input')
            end
            obj.mogs = Mog.empty;
            obj.air = AirShots.empty;
            obj.handles.hp = uipanel(varargin{:},...
                'Title','MOGs',...
                'Visible','off',...
                'SizeChangedFcn', @obj.resizeUI);
            obj.data_rep = '';
            obj.addComponents();
            obj.handles.hp.Visible = 'on';
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
        function set.air(obj, a)
            if isa(a, 'AirShots')
                obj.air = a;
            else
                error('air should be objects of class AirShots')
            end
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
            obj.handles.mergeMogs.FontSize = s;
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
            no = obj.handles.listMogs.Value;
            if no>0 && no<=length(obj.mogs)
                obj.handles.popupTx.Value = obj.mogs(no).Tx;
                obj.handles.popupRx.Value = obj.mogs(no).Rx;     
            else
                obj.handles.popupTx.Value = 1;
                obj.handles.popupRx.Value = 1;
            end
        end
    end

    methods (Access=private)
        addComponents(obj)
        resizeUI(obj,varargin)
        addMog(obj,varargin)
        importMog(obj,varargin)
        updateCoords(obj,varargin)
        statsTt(obj,varargin)
        airBefore(obj,varargin)
        airAfter(obj,varargin)
        prune(obj,varargin)
        spectra(obj,varargin)
        zop(obj,varargin)
        mergeMogs(obj,varargin)
        
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
        function popupType(obj,varargin)
            no = obj.handles.listMogs.Value;
            if no>0 && no<=length(obj.mogs)
                obj.mogs(no).type = obj.handles.popupType.type;
                obj.updateCoords()
                obj.notify('mogEdited')
            end
        end
        function popupTx(obj,varargin)
            no = obj.handles.listMogs.Value;
            if no>0 && no<=length(obj.mogs)
                obj.mogs(no).Tx = obj.handles.popupTx.Value;
                obj.updateCoords(no);
                obj.notify('mogEdited')
            end
        end
        function popupRx(obj,varargin)
            no = obj.handles.listMogs.Value;
            if no>0 && no<=length(obj.mogs)
                obj.mogs(no).Rx = obj.handles.popupRx.Value;
                obj.updateCoords(no);
                obj.notify('mogEdited')
            end
        end
        function renameMog(obj,varargin)
            no = obj.handles.listMogs.Value;
            if no>0 && no<=length(obj.mogs)
                name = myinputdlg('Enter new name');
                if ~isempty(name)
                    obj.mogs(no).name = char(name);
                    obj.updateList()
                    obj.notify('mogEdited')
                end
            end
        end
        function rawData(obj,varargin)
            no = obj.handles.listMogs.Value;
            if no>0 && no<=length(obj.mogs)
                figure
                in = obj.mogs(no).in;
                imagesc(1:obj.mogs(no).no_traces(in),...
                    obj.mogs(no).data.timestp,...
                    obj.mogs(no).data.rdata(:,in))
                clim = caxis;
                cmax = max(abs(clim));
                caxis([-cmax cmax])
                colormap ramac_cmap(16)
                colorbar('peer',gca);
                xlabel('Trace no')
                ylabel('Time [ns]')
            end
        end
        function statsAmp(obj,varargin)
            no = obj.handles.listMogs.Value;
            if no>0 && no<=length(obj.mogs)
                mog = obj.mogs(no);
                hyp = sqrt( (mog.data.Tx_x-mog.data.Rx_x).^2 + ...
                    (mog.data.Tx_y-mog.data.Rx_y).^2 + ...
                    (mog.data.Tx_z-mog.data.Rx_z).^2 );
                dz = mog.data.Rx_z-mog.data.Tx_z;
                theta = 180/pi*asin(dz./hyp);
                
                h_stat = figure;
                set(h_stat, 'Position',[256 71 512 620]);
                
                ind = mog.amp_tmax ~= -1 & mog.tauApp ~= -1 & mog.amp_done;
                ind = ind & mog.in;
                
                subplot(321)
                plot(hyp(ind), mog.tauApp(ind),'o')
                xlabel('Straight ray length [m]')
                ylabel('\tau_a')
                title('amplitude - amp. ratio')
                
                subplot(322)
                plot(theta(ind), mog.tauApp(ind)./hyp(ind),'o')
                xlabel('Angle w/r to horizontal [deg]')
                ylabel('\alpha_a')
                title('amplitude - amp. ratio')
                
                ind = mog.amp_tmax ~= -1 & mog.tauFce ~= -1 & mog.amp_done;
                ind = ind & mog.in;
                
                subplot(323)
                plot(hyp(ind), mog.tauFce(ind),'o')
                xlabel('Straight ray length [m]')
                ylabel('\tau_a')
                title('amplitude - centroid freq.')
                
                subplot(324)
                plot(theta(ind), mog.tauFce(ind)./hyp(ind),'o')
                xlabel('Angle w/r to horizontal [deg]')
                ylabel('\alpha_a')
                title('amplitude - centroid freq.')
                
                ind = mog.amp_tmax ~= -1 & mog.tauHyb ~= -1 & mog.amp_done;
                ind = ind & mog.in;
                
                subplot(325)
                plot(hyp(ind), mog.tauHyb(ind),'o')
                xlabel('Straight ray length [m]')
                ylabel('\tau_a')
                title('amplitude - hybrid')
                
                subplot(326)
                plot(theta(ind), mog.tauHyb(ind)./hyp(ind),'o')
                xlabel('Angle w/r to horizontal [deg]')
                ylabel('\alpha_a')
                title('amplitude - hybrid')
            end
        end
        function rays(obj,varargin)
            no = obj.handles.listMogs.Value;
            if no>0 && no<=length(obj.mogs)
                mog = obj.mogs(no);
                boreholes = obj.bhUI.boreholes;
                ind1 = mog.tt==-1 & mog.in;
                ind2 = mog.tt~=-1 & mog.in;
                
                [~, a]=lsplane([boreholes(mog.Tx).fdata; boreholes(mog.Rx).fdata]);
                el = (pi-a(3))*180/pi - 15;
                az = atan( cos(a(2))/cos(a(1)) )*180/pi - 10;
                
                figure
                plot3([mog.data.Tx_x(ind1); mog.data.Rx_x(ind1)], ...
                    [mog.data.Tx_y(ind1); mog.data.Rx_y(ind1)], ...
                    [mog.data.Tx_z(ind1); mog.data.Rx_z(ind1)],'r')
                hold on
                plot3([mog.data.Tx_x(ind2); mog.data.Rx_x(ind2)], ...
                    [mog.data.Tx_y(ind2); mog.data.Rx_y(ind2)], ...
                    [mog.data.Tx_z(ind2); mog.data.Rx_z(ind2)],'g')
                hold off
                set(gca, 'DataAspectRatio',[1 1 1],'Units','normalized')
                title([num2str(100*sum(ind2)/length(ind2)), '%'])
                axis tight
                xlabel('Tx-Rx Distance [m]')
                zlabel('Elevation [m]')
                view(az,el)
            end
        end
        function exportTt(obj,varargin)
            no = obj.handles.listMogs.Value;
            if no>0 && no<=length(obj.mogs)
                [filename, pathname] = uiputfile({'*.dat';'*.*'}, 'Export traveltimes');
                if isequal(filename,0) || isequal(pathname,0)
                    return;
                end
                mog = obj.mogs(no);
                done = mog.tt_done & mog.in;
                if sum(done)==0
                    errordlg('Data not processed!!!')
                    return
                end
                ind = mog.tt~=-1 & mog.in;
                [tt,~] = mog.getCorrectedTravelTimes(obj.air);
                obj.handles.editDtCorr.String = num2str( mog.fac_dt );
                et = mog.et(ind);
                tt = tt(ind);
                
                data = [mog.data.Tx_x(ind)' mog.data.Tx_y(ind)' mog.data.Tx_z(ind)' ...
                    mog.data.Rx_x(ind)' mog.data.Rx_y(ind)' mog.data.Rx_z(ind)' ...
                    tt' et' mog.no_traces(ind)']; %#ok<NASGU>
                
                eval(['save ',pathname,filename,' data -ascii'])
            end
        end
        function exportTau(obj,varargin)
            no = obj.handles.listMogs.Value;
            if no>0 && no<=length(obj.mogs)
                mog = obj.mogs(no);
                tt = mog.tauApp;
                et = mog.tauApp_et;
                ind = tt~=-1 & mog.in;
                tt = tt(ind);
                et = et(ind);
                if isempty(tt)
                    errordlg('Data not processed!!!')
                    return
                end
                [filename, pathname] = uiputfile({'*.dat';'*.*'}, 'Export reduced amplitudes');
                if isequal(filename,0) || isequal(pathname,0)
                    return;
                end
                data = [mog.data.Tx_x(ind)' mog.data.Tx_y(ind)' mog.data.Tx_z(ind)' ...
                    mog.data.Rx_x(ind)' mog.data.Rx_y(ind)' mog.data.Rx_z(ind)' ...
                    tt' et' mog.no_traces(ind)']; %#ok<NASGU>
                
                eval(['save ',pathname,filename,' data -ascii'])
            end
        end
        function editDate(obj,varargin)
            no = obj.handles.listMogs.Value;
            if no>0 && no<=length(obj.mogs)
                obj.mogs(no).data = obj.handles.editDate.String;
                obj.notify('mogEdited')
            end
        end
        function editFreq(obj,varargin)
            no = obj.handles.listMogs.Value;
            if no>0 && no<=length(obj.mogs)
                obj.mogs(no).data.rnomfreq = str2double(obj.handles.editFreq.String);
                obj.notify('mogEdited')
            end
        end
        function editFeedRx(obj,varargin)
            no = obj.handles.listMogs.Value;
            if no>0 && no<=length(obj.mogs)
                obj.mogs(no).data.RxOffset = str2double(obj.handles.editFeedRx.String);
                obj.notify('mogEdited')
            end
        end
        function editFeedTx(obj,varargin)
            no = obj.handles.listMogs.Value;
            if no>0 && no<=length(obj.mogs)
                obj.mogs(no).data.TxOffset = str2double(obj.handles.editFeedTx.String);
                obj.notify('mogEdited')
            end
        end
        function checkDtCorr(obj,varargin)
            no = obj.handles.listMogs.Value;
            if no>0 && no<=length(obj.mogs)
                obj.mogs(no).user_fac_dt = obj.handles.checkDtCorr.Value;
                obj.notify('mogEdited')
            end
        end
        function editDtCorr(obj,varargin)
            no = obj.handles.listMogs.Value;
            if no>0 && no<=length(obj.mogs)
                obj.mogs(no).fac_dt = str2double(obj.handles.editDtCorr.String);
                obj.notify('mogEdited')
            end
        end
        function editMultFac(obj,varargin)
            no = obj.handles.listMogs.Value;
            if no>0 && no<=length(obj.mogs)
                obj.mogs(no).f_et = str2double(obj.handles.editMultFac.String);
                obj.notify('mogEdited')
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
                if ~isempty(obj.mogs(no).av)
                    obj.handles.textAirBefore.String = obj.air(obj.mogs(no).av).name;
                end
                if ~isempty(obj.mogs(no).ap)
                    obj.handles.textAirAfter.String = obj.air(obj.mogs(no).ap).name;
                end
            end
        end
    end
end

