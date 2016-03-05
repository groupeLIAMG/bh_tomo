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
    events
        mogAdded, mogDeleted
    end
    
    methods
        function obj = MogUI(varargin)
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
            boreholes = varargin{1}.boreholes;
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
        function addComponents(obj)
            
            obj.handles.addMog = uicontrol('Style','pushbutton',...
                'String','Add MOG',...
                'Units','points',...
                'Callback',@obj.addMog,...
                'Parent',obj.handles.hp);
            obj.handles.removeMog = uicontrol('Style','pushbutton',...
                'String','Remove MOG',...
                'Units','points',...
                'Callback',@obj.removeMog,...
                'Parent',obj.handles.hp);
            
            obj.handles.listMogs = uicontrol('Style','listbox',...
                'Max',1,'Min',0,...
                'Units','points',...
                'Callback',@obj.listMogs,...
                'Parent',obj.handles.hp);
            
            obj.handles.textType = uicontrol('Style','text',...
                'String','Type',...
                'HorizontalAlignment','right',...
                'Units','points',...
                'Parent',obj.handles.hp);
            obj.handles.textTx = uicontrol('Style','text',...
                'String','Tx',...
                'HorizontalAlignment','right',...
                'Units','points',...
                'Parent',obj.handles.hp);
            obj.handles.textRx = uicontrol('Style','text',...
                'String','Rx',...
                'HorizontalAlignment','right',...
                'Units','points',...
                'Parent',obj.handles.hp);

            obj.handles.popupType = uicontrol('Style','popupmenu',...
                'String',{'Crosshole','VSP/VRP'},...
                'Units','points',...
                'Callback',@obj.popupType,...
                'Parent',obj.handles.hp);
            obj.handles.popupTx = uicontrol('Style','popupmenu',...
                'String',' ',...
                'Units','points',...
                'Callback',@obj.popupTx,...
                'Parent',obj.handles.hp);
            obj.handles.popupRx = uicontrol('Style','popupmenu',...
                'String',' ',...
                'Units','points',...
                'Callback',@obj.popupRx,...
                'Parent',obj.handles.hp);

            obj.handles.airBefore = uicontrol('Style','pushbutton',...
                'String','Air Shot Before',...
                'Units','points',...
                'Callback',@obj.airBefore,...
                'Parent',obj.handles.hp);
            obj.handles.airAfter = uicontrol('Style','pushbutton',...
                'String','Air Shot After',...
                'Units','points',...
                'Callback',@obj.airAfter,...
                'Parent',obj.handles.hp);
            obj.handles.textAirBefore = uicontrol('Style','text',...
                'BackgroundColor',[1 1 1],...
                'Units','points',...
                'Parent',obj.handles.hp);
            obj.handles.textAirAfter = uicontrol('Style','text',...
                'BackgroundColor',[1 1 1],...
                'Units','points',...
                'Parent',obj.handles.hp);
            
            obj.handles.renameMog = uicontrol('Style','pushbutton',...
                'String','Rename',...
                'Units','points',...
                'Callback',@obj.renameMog,...
                'Parent',obj.handles.hp);
            obj.handles.importMog = uicontrol('Style','pushbutton',...
                'String','Import',...
                'Units','points',...
                'Callback',@obj.importMog,...
                'Parent',obj.handles.hp);
            obj.handles.spectra = uicontrol('Style','pushbutton',...
                'String','Spectra',...
                'Units','points',...
                'Callback',@obj.spectra,...
                'Parent',obj.handles.hp);
            obj.handles.rawData = uicontrol('Style','pushbutton',...
                'String','Raw Data',...
                'Units','points',...
                'Callback',@obj.rawData,...
                'Parent',obj.handles.hp);
            obj.handles.zop = uicontrol('Style','pushbutton',...
                'String','Trace ZOP',...
                'Units','points',...
                'Callback',@obj.zop,...
                'Parent',obj.handles.hp);

            obj.handles.statsTt = uicontrol('Style','pushbutton',...
                'String','Stats tt',...
                'Units','points',...
                'Callback',@obj.statsTt,...
                'Parent',obj.handles.hp);
            obj.handles.statsAmp = uicontrol('Style','pushbutton',...
                'String','Stats Ampl.',...
                'Units','points',...
                'Callback',@obj.statsAmp,...
                'Parent',obj.handles.hp);
            obj.handles.rays = uicontrol('Style','pushbutton',...
                'String','Ray Coverage',...
                'Units','points',...
                'Callback',@obj.rays,...
                'Parent',obj.handles.hp);
            obj.handles.exportTt = uicontrol('Style','pushbutton',...
                'String','Export tt',...
                'Units','points',...
                'Callback',@obj.exportTt,...
                'Parent',obj.handles.hp);
            obj.handles.exportTau = uicontrol('Style','pushbutton',...
                'String',['Export ' char(964)],...
                'Units','points',...
                'Callback',@obj.exportTau,...
                'Parent',obj.handles.hp);
            obj.handles.prune = uicontrol('Style','pushbutton',...
                'String','Prune',...
                'Units','points',...
                'Callback',@obj.prune,...
                'Parent',obj.handles.hp);

            obj.handles.textDate = uicontrol('Style','text',...
                'String','Date',...
                'HorizontalAlignment','right',...
                'Units','points',...
                'Parent',obj.handles.hp);
            obj.handles.editDate = uicontrol('Style','edit',...
                'Units','points',...
                'Callback',@obj.editDate,...
                'Parent',obj.handles.hp);
            
            obj.handles.textFreq = uicontrol('Style','text',...
                'String','Nominal Frequency of Antennas',...
                'HorizontalAlignment','right',...
                'Units','points',...
                'Parent',obj.handles.hp);
            obj.handles.editFreq = uicontrol('Style','edit',...
                'Units','points',...
                'Callback',@obj.editFreq,...
                'Parent',obj.handles.hp);
            obj.handles.textFeedRx = uicontrol('Style','text',...
                'String','Antenna Feedpoint Offset - Rx',...
                'HorizontalAlignment','right',...
                'Units','points',...
                'Parent',obj.handles.hp);
            obj.handles.editFeedRx = uicontrol('Style','edit',...
                'Units','points',...
                'Callback',@obj.editFeedRx,...
                'Parent',obj.handles.hp);
            obj.handles.textFeedTx = uicontrol('Style','text',...
                'String','Antenna Feedpoint Offset - Tx',...
                'HorizontalAlignment','right',...
                'Units','points',...
                'Parent',obj.handles.hp);
            obj.handles.editFeedTx = uicontrol('Style','edit',...
                'Units','points',...
                'Callback',@obj.editFeedTx,...
                'Parent',obj.handles.hp);
            obj.handles.checkDtCorr = uicontrol('Style','checkbox',...
                'String','Time Step Correction Factor',...
                'Units','points',...
                'Callback',@obj.checkDtCorr,...
                'Parent',obj.handles.hp);
            obj.handles.editDtCorr = uicontrol('Style','edit',...
                'Units','points',...
                'Callback',@obj.editDtCorr,...
                'Parent',obj.handles.hp);
            obj.handles.textMultFac = uicontrol('Style','text',...
                'String','Std dev. Multiplication Factor',...
                'HorizontalAlignment','right',...
                'Units','points',...
                'Parent',obj.handles.hp);
            obj.handles.editMultFac = uicontrol('Style','edit',...
                'Units','points',...
                'Callback',@obj.editMultFac,...
                'Parent',obj.handles.hp);

        end
        function resizeUI(obj,varargin)
            oldUnits = obj.handles.hp.Units;
            obj.handles.hp.Units = 'points';
            obj.handles.hp.Visible = 'off';

            p = obj.handles.hp.Position;
            width = p(3);  % prefered: 630
            height = p(4); % prefered: 380

            hSize = width/7;
            hSpace = width/120;
            hBorder = width/30;

            vFac = 0.8*height/360;
            if vFac<1
                vFac = 1;
            end
            vSize = 22*vFac;
            vSpace = 5*vFac;
            vBorderTop = 45*vFac;
            vBorder = 15*vFac;
            
            obj.handles.addMog.Position = [hBorder height-vBorderTop 1.5*hSize vSize];
            obj.handles.removeMog.Position = [hBorder+2*hSpace+1.5*hSize height-vBorderTop 1.5*hSize vSize];

            obj.handles.listMogs.Position = [1.5*hBorder vBorder+5*vSpace+5.5*vSize 2*hSpace+3*hSize-hBorder height-2*vBorder-vBorderTop-6.5*vSize];

            obj.handles.textType.Position = [hBorder vBorder+4*vSpace+4*vSize 1.25*hSize vSize];
            obj.handles.textTx.Position = [hBorder vBorder+3*vSpace+3*vSize 1.25*hSize vSize];
            obj.handles.textRx.Position = [hBorder vBorder+2*vSpace+2*vSize 1.25*hSize vSize];
            obj.handles.popupType.Position = [hBorder+2*hSpace+1.5*hSize vBorder+4*vSpace+4*vSize 1.5*hSize vSize];
            obj.handles.popupTx.Position = [hBorder+2*hSpace+1.5*hSize vBorder+3*vSpace+3*vSize 1.5*hSize vSize];
            obj.handles.popupRx.Position = [hBorder+2*hSpace+1.5*hSize vBorder+2*vSpace+2*vSize 1.5*hSize vSize];
            
            obj.handles.airBefore.Position = [hBorder vBorder+vSpace+vSize 1.5*hSize vSize];
            obj.handles.airAfter.Position = [hBorder vBorder 1.5*hSize vSize];
            obj.handles.textAirBefore.Position = [hBorder+2*hSpace+1.5*hSize vBorder+vSpace+vSize 1.5*hSize vSize];
            obj.handles.textAirAfter.Position = [hBorder+2*hSpace+1.5*hSize vBorder 1.5*hSize vSize];

            obj.handles.renameMog.Position = [width-hBorder-3*hSize-2*hSpace height-vBorderTop-vSize hSize vSize];
            obj.handles.importMog.Position = [width-hBorder-2*hSize-hSpace   height-vBorderTop-vSize hSize vSize];
            obj.handles.rawData.Position   = [width-hBorder-3*hSize-2*hSpace height-vBorderTop-vSpace-2*vSize hSize vSize];
            obj.handles.zop.Position       = [width-hBorder-2*hSize-hSpace   height-vBorderTop-vSpace-2*vSize hSize vSize];
            obj.handles.spectra.Position   = [width-hBorder-hSize            height-vBorderTop-vSpace-2*vSize hSize vSize];
            obj.handles.statsTt.Position   = [width-hBorder-3*hSize-2*hSpace height-vBorderTop-2*vSpace-3*vSize hSize vSize];
            obj.handles.statsAmp.Position  = [width-hBorder-2*hSize-hSpace   height-vBorderTop-2*vSpace-3*vSize hSize vSize];
            obj.handles.rays.Position      = [width-hBorder-hSize            height-vBorderTop-2*vSpace-3*vSize hSize vSize];
            obj.handles.exportTt.Position  = [width-hBorder-3*hSize-2*hSpace height-vBorderTop-3*vSpace-4*vSize hSize vSize];
            obj.handles.exportTau.Position = [width-hBorder-2*hSize-hSpace   height-vBorderTop-3*vSpace-4*vSize hSize vSize];
            obj.handles.prune.Position     = [width-hBorder-hSize            height-vBorderTop-3*vSpace-4*vSize hSize vSize];

            obj.handles.textDate.Position = [width-hBorder-2.5*hSize-hSpace vBorder hSize vSize];
            obj.handles.editDate.Position = [width-hBorder-1.5*hSize vBorder 1.5*hSize vSize];

            obj.handles.textFreq.Position = [width-hBorder-3*hSize-hSpace vBorderTop+5*(vSpace+vSize) 2*hSize vSize];
            obj.handles.editFreq.Position = [width-hBorder-hSize vBorderTop+5*(vSpace+vSize) hSize vSize];
            obj.handles.textFeedRx.Position = [width-hBorder-3*hSize-hSpace vBorderTop+4*(vSpace+vSize) 2*hSize vSize];
            obj.handles.editFeedRx.Position = [width-hBorder-hSize vBorderTop+4*(vSpace+vSize) hSize vSize];
            obj.handles.textFeedTx.Position = [width-hBorder-3*hSize-hSpace vBorderTop+3*(vSpace+vSize) 2*hSize vSize];
            obj.handles.editFeedTx.Position = [width-hBorder-hSize vBorderTop+3*(vSpace+vSize) hSize vSize];
            ext = obj.handles.checkDtCorr.Extent;
            obj.handles.checkDtCorr.Position = [width-hBorder-1.15*ext(3)-hSize-hSpace vBorderTop+2*(vSpace+vSize) 1.15*ext(3) vSize];
            obj.handles.editDtCorr.Position = [width-hBorder-hSize vBorderTop+2*(vSpace+vSize) hSize vSize];
            obj.handles.textMultFac.Position = [width-hBorder-3*hSize-hSpace vBorderTop+vSpace+vSize 2*hSize vSize];
            obj.handles.editMultFac.Position = [width-hBorder-hSize vBorderTop+vSpace+vSize hSize vSize];

            obj.handles.hp.Visible = 'on';
            obj.handles.hp.Units = oldUnits;
        end
        
        function addMog(obj,varargin)
            if strcmp(obj.handles.popupTx.String,' ')
                warndlg('Define boreholes first.')
                return
            end
            old_rep = pwd;
            if ~isempty(obj.data_rep)
                cd( obj.data_rep );
            end
            [file, rep, filterindex] = uigetfile({
                '*.rad;*.RAD', 'RAMAC file *.rad';...
                '*.mat', 'Magnetic NanoMappers data in Matlab format (*.mat)';...
                '*.hd;*.HD', 'EKKO file *.hd';...
                '*.mat', 'RAMAC file in Matlab format *.mat';...
                '*.dat', 'MSIS file *.dat';...
                '*.su;*.SU', 'Seismic Un*x file *.su';...
                '*.sgy;*.segy;*.SGY;*.SEGY', 'SEG Y file'},...
                'Open data file');
            cd( old_rep );
            if isequal(file,0) || isequal(rep,0)
                return
            end
            obj.data_rep = rep;

            no = length(obj.mogs)+1;
            switch filterindex
                case 1
                    name = file(1:end-4);
                    obj.mogs(no) = Mog(name);
                    obj.mogs(no).data = MogData();
                    obj.mogs(no).data.readRAMAC([rep,name]);
                case 2
                    name = file(1:end-4);
                    obj.mogs(no) = Mog(name);
                    obj.mogs(no).data = MogData();
                    obj.mogs(no).data.readNanoMapper([rep,name]);
                case 3
                    name = file(1:end-3);
                    obj.mogs(no) = Mog(name);
                    obj.mogs(no).data = MogData();
                    obj.mogs(no).data.readEKKO([rep,name]);
                case 4
                    name = file(1:end-4);
                    obj.mogs(no) = Mog(name);
                    obj.mogs(no).data = MogData();
                    obj.mogs(no).data.readMAT([rep,name]);
                case 5
                    name = file(1:end-4);
                    obj.mogs(no) = Mog(name);
                    obj.mogs(no).data = MogData();
                    obj.mogs(no).data.readMSIS([rep,name]);
                case 6
                    name = file(1:end-4);
                    obj.mogs(no) = Mog(name);
                    obj.mogs(no).data = MogData();
                    obj.mogs(no).data.readSU([rep,name]);
                case 7
                    name = file(1:end-4);
                    obj.mogs(no) = Mog(name);
                    obj.mogs(no).data = MogData();
                    obj.mogs(no).data.readSEGY([rep,name]);
            end
            obj.mogs(no).initialize();
            
            debut=0;
            if no>1
                debut = obj.mogs(no-1).no_traces( length(obj.mogs(no-1).no_traces) );
            end
            obj.mogs(no).no_traces = debut + (1:obj.mogs(no).data.ntrace);
            obj.mogs(no).sorted = false;

            
            obj.updateList()
        end
        function removeMog(obj,varargin)
        end
        function listMogs(obj,varargin)
        end
        function airBefore(obj,varargin)
        end
        function airAfter(obj,varargin)
        end
        function popupType(obj,varargin)
        end
        function popupTx(obj,varargin)
        end
        function popupRx(obj,varargin)
        end
        function renameMog(obj,varargin)
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
        end
        function editFreq(obj,varargin)
        end
        function editFeedRx(obj,varargin)
        end
        function editFeedTx(obj,varargin)
        end
        function checkDtCorr(obj,varargin)
        end
        function editDtCorr(obj,varargin)
        end
        function editMultFac(obj,varargin)
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
        end
    end
end

