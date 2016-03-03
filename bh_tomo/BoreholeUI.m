classdef BoreholeUI < handle
    %BOREHOLEUI User interface to manage Boreholes
    
    properties (SetAccess=private)
        boreholes
    end
    properties (Dependent)
        Position
        FontSize   % should remain be within [10 11 12]
    end
    properties (Access=private)
        handles
    end
    events
        coordinatesChanged, boreholeAdded, boreholeDeleted
    end
    methods
        function obj = BoreholeUI(varargin)
            obj.handles.hp = uipanel(varargin{:},...
            'Title','Boreholes',...
            'Visible','off',...
            'SizeChangedFcn', @obj.resizeUI);
            
            obj.addComponents();
            obj.resizeUI();
        end
        function set.Position(obj, p)
            obj.handles.hp.Position = p;
        end
        function set.FontSize(obj, s)
            obj.handles.hp.FontSize = s+1;
            obj.handles.addBH.FontSize = s;
            obj.handles.removeBH.FontSize = s;
            obj.handles.importBH.FontSize = s;
            obj.handles.plotBH.FontSize = s;
            obj.handles.listBH.FontSize = s;
            obj.handles.textCoord.FontSize = s;
            obj.handles.textX.FontSize = s;
            obj.handles.textY.FontSize = s;
            obj.handles.textZ.FontSize = s;
            obj.handles.textCollar.FontSize = s;
            obj.handles.textBottom.FontSize = s;
            obj.handles.topX.FontSize = s;
            obj.handles.topY.FontSize = s;
            obj.handles.topZ.FontSize = s;
            obj.handles.bottomX.FontSize = s;
            obj.handles.bottomY.FontSize = s;
            obj.handles.bottomZ.FontSize = s;
            obj.handles.textZwater.FontSize = s;
            obj.handles.Zwater.FontSize = s;
            obj.handles.textZsurf.FontSize = s;
            obj.handles.Zsurf.FontSize = s;
            obj.handles.textDiam.FontSize = s;
            obj.handles.diameter.FontSize = s;
            obj.handles.contS.FontSize = s;
            obj.handles.contA.FontSize = s;
            obj.resizeUI();
        end
    end
    
    methods (Access=private)
        function addComponents(obj)
            
            obj.handles.addBH = uicontrol('Style','pushbutton',...
                'String','Add',...
                'Units','points',...
                'Callback',@obj.addBH,...
                'Parent',obj.handles.hp);
            obj.handles.removeBH = uicontrol('Style','pushbutton',...
                'String','Remove',...
                'Units','points',...
                'Callback',@obj.removeBH,...
                'Parent',obj.handles.hp);
            obj.handles.importBH = uicontrol('Style','pushbutton',...
                'String','Import',...
                'Units','points',...
                'Callback',@obj.importBH,...
                'Parent',obj.handles.hp);
            obj.handles.plotBH = uicontrol('Style','pushbutton',...
                'String','Plot',...
                'Units','points',...
                'Callback',@obj.plotBH,...
                'Parent',obj.handles.hp);

            obj.handles.listBH = uicontrol('Style','listbox',...
                'Max',1,'Min',0,...
                'Units','points',...
                'Callback',@obj.listBH,...
                'Parent',obj.handles.hp);
            
            obj.handles.textCoord = uicontrol('Style','text',...
                'String','Coordinates',...
                'Units','points',...
                'HorizontalAlignment','right',...
                'Parent',obj.handles.hp);
            obj.handles.textX = uicontrol('Style','text',...
                'String','X:',...
                'Units','points',...
                'HorizontalAlignment','right',...
                'Parent',obj.handles.hp);
            obj.handles.textY = uicontrol('Style','text',...
                'String','Y:',...
                'Units','points',...
                'HorizontalAlignment','right',...
                'Parent',obj.handles.hp);
            obj.handles.textZ = uicontrol('Style','text',...
                'String','Elevation:',...
                'Units','points',...
                'HorizontalAlignment','right',...
                'Parent',obj.handles.hp);
            obj.handles.textCollar = uicontrol('Style','text',...
                'String','Collar',...
                'Units','points',...
                'HorizontalAlignment','center',...
                'Parent',obj.handles.hp);
            obj.handles.textBottom = uicontrol('Style','text',...
                'String','Bottom',...
                'Units','points',...
                'HorizontalAlignment','center',...
                'Parent',obj.handles.hp);
            
            obj.handles.topX = uicontrol('Style','edit',...
                'Units','points',...
                'Callback',@obj.topX,...
                'Parent',obj.handles.hp);
            obj.handles.topY = uicontrol('Style','edit',...
                'Units','points',...
                'Callback',@obj.topY,...
                'Parent',obj.handles.hp);
            obj.handles.topZ = uicontrol('Style','edit',...
                'Units','points',...
                'Callback',@obj.topZ,...
                'Parent',obj.handles.hp);
            obj.handles.bottomX = uicontrol('Style','edit',...
                'Units','points',...
                'Callback',@obj.bottomX,...
                'Parent',obj.handles.hp);
            obj.handles.bottomY = uicontrol('Style','edit',...
                'Units','points',...
                'Callback',@obj.bottomY,...
                'Parent',obj.handles.hp);
            obj.handles.bottomZ = uicontrol('Style','edit',...
                'Units','points',...
                'Callback',@obj.bottomZ,...
                'Parent',obj.handles.hp);

            obj.handles.textZsurf = uicontrol('Style','text',...
                'String','Elevation at surface:',...
                'Units','points',...
                'HorizontalAlignment','right',...
                'Parent',obj.handles.hp);
            obj.handles.Zsurf = uicontrol('Style','edit',...
                'Units','points',...
                'Callback',@obj.Zsurf,...
                'Parent',obj.handles.hp);

            obj.handles.textZwater = uicontrol('Style','text',...
                'String','Elev. water:',...
                'Units','points',...
                'HorizontalAlignment','right',...
                'Parent',obj.handles.hp);
            obj.handles.Zwater = uicontrol('Style','edit',...
                'Units','points',...
                'Callback',@obj.Zwater,...
                'Parent',obj.handles.hp);

            obj.handles.textDiam = uicontrol('Style','text',...
                'String','Diameter:',...
                'Units','points',...
                'HorizontalAlignment','right',...
                'Parent',obj.handles.hp);
            obj.handles.diameter = uicontrol('Style','edit',...
                'Units','points',...
                'Callback',@obj.diameter,...
                'Parent',obj.handles.hp);

            obj.handles.contS = uicontrol('Style','pushbutton',...
                'String','Constraints slown.',...
                'Units','points',...
                'Callback',@obj.contS,...
                'Parent',obj.handles.hp);
            obj.handles.contA = uicontrol('Style','pushbutton',...
                'String','Constraints atten.',...
                'Units','points',...
                'Callback',@obj.contA,...
                'Parent',obj.handles.hp);

            obj.handles.hp.Visible = 'on';
        end
        function resizeUI(obj,varargin)
            oldUnits = obj.handles.hp.Units;
            obj.handles.hp.Units = 'points';
            obj.handles.hp.Visible = 'off';

            p = obj.handles.hp.Position;
            width = p(3);  % prefered: 300
            height = p(4); % prefered: 360

            vSize = 20;
            hSize = width/5;
            hSpace = width/60;
            vSpace = 5;
            hBorder = width/15;
            vBorderTop = height/8;
            vBorder = 15;
            
            vSizeList = height-vBorderTop-2*vSpace-(vBorder+7*vSize+7*vSpace);
            
            obj.handles.addBH.Position = [hBorder height-vBorderTop hSize vSize];
            obj.handles.removeBH.Position = [hBorder+hSize+hSpace height-vBorderTop hSize vSize];
            obj.handles.importBH.Position = [hBorder+2*hSize+2*hSpace height-vBorderTop hSize vSize];
            obj.handles.plotBH.Position = [hBorder+3*hSize+3*hSpace height-vBorderTop hSize vSize];
            
            obj.handles.listBH.Position = [1.5*hBorder vBorder+7*vSize+7*vSpace width-3*hBorder vSizeList];
            obj.handles.textCoord.Position = [width/2-3/2*hSize-hSpace vBorder+6*vSize+6*vSpace hSize vSize];
            obj.handles.textX.Position = [width/2-3/2*hSize-hSpace vBorder+5*vSize+5*vSpace hSize vSize];
            obj.handles.textY.Position = [width/2-3/2*hSize-hSpace vBorder+4*vSize+4*vSpace hSize vSize];
            obj.handles.textZ.Position = [width/2-3/2*hSize-hSpace vBorder+3*vSize+3*vSpace hSize vSize];
            obj.handles.textCollar.Position = [width/2-hSize/2 vBorder+6*vSize+6*vSpace hSize vSize];
            obj.handles.textBottom.Position = [width/2+hSize/2+hSpace vBorder+6*vSize+6*vSpace hSize vSize];
            obj.handles.topX.Position = [width/2-hSize/2 vBorder+5*vSize+5*vSpace hSize vSize];
            obj.handles.topY.Position = [width/2-hSize/2 vBorder+4*vSize+4*vSpace hSize vSize];
            obj.handles.topZ.Position = [width/2-hSize/2 vBorder+3*vSize+3*vSpace hSize vSize];
            obj.handles.bottomX.Position = [width/2+hSize/2+hSpace vBorder+5*vSize+5*vSpace hSize vSize];
            obj.handles.bottomY.Position = [width/2+hSize/2+hSpace vBorder+4*vSize+4*vSpace hSize vSize];
            obj.handles.bottomZ.Position = [width/2+hSize/2+hSpace vBorder+3*vSize+3*vSpace hSize vSize];
            
            obj.handles.textZsurf.Position = [width/2-5/3*hSize-hSpace vBorder+2*vSize+2*vSpace 5/3*hSize vSize];
            obj.handles.Zsurf.Position = [width/2 vBorder+2*vSize+2*vSpace hSize vSize];
            
            obj.handles.textZwater.Position = [hBorder vBorder+vSize+vSpace hSize vSize];
            obj.handles.Zwater.Position = [hBorder+hSize+hSpace vBorder+vSize+vSpace hSize vSize];
            obj.handles.textDiam.Position = [hBorder+2*hSize+2*hSpace vBorder+vSize+vSpace 5/6*hSize vSize];
            obj.handles.diameter.Position = [hBorder+17/6*hSize+3*hSpace vBorder+vSize+vSpace hSize vSize];
            
            obj.handles.contS.Position = [width/2-2*hSize-hSpace-2 vBorder 2*hSize vSize];
            obj.handles.contA.Position = [width/2+hSpace vBorder 2*hSize vSize];
            
            obj.handles.hp.Visible = 'on';
            obj.handles.hp.Units = oldUnits;
        end
        
        function addBH(obj,varargin)
            name = inputdlg('Borehole name');
            if isempty(name)
                return;
            end
            if isempty(obj.boreholes)
                obj.boreholes = Borehole(char(name));
            else
                obj.boreholes(end+1) = Borehole(char(name));
            end
            names = cell(1,length(obj.boreholes));
            for n=1:length(obj.boreholes)
                names{n} = obj.boreholes(n).name;
            end
            obj.handles.listBH.String = names;
            obj.handles.listBH.Value = length(obj.boreholes);
            obj.updateEdits
            notify(obj,'boreholeAdded')
        end
        function removeBH(obj, varargin)
            no = obj.handles.listBH.Value;
            ind=1:length(obj.boreholes);
            ind = ind~=no;

            obj.boreholes = obj.boreholes(ind);
            names = cell(1,length(obj.boreholes));
            for n=1:length(obj.boreholes)
                names{n} = obj.boreholes(n).name;
            end
            obj.handles.listBH.String = names;
            obj.handles.listBH.Value = length(obj.boreholes);
            obj.updateEdits
            notify(obj,'boreholeDeleted')
        end
        function importBH(obj, varargin)
            [filename, pathname] = uigetfile('*.xyz', 'Import borehole');
            if isequal(filename,0) || isequal(pathname,0)
                return
            else
                name = filename(1:end-4);
                fdata = load([pathname,filename]);
            end
            if isempty(obj.boreholes)
                obj.boreholes = Borehole(char(name));
            else
                obj.boreholes(end+1) = Borehole(char(name));
            end
            obj.boreholes(end).X = fdata(1,1);
            obj.boreholes(end).Y = fdata(1,2);
            obj.boreholes(end).Z = fdata(1,3);
            obj.boreholes(end).Xmax = fdata(end,1);
            obj.boreholes(end).Ymax = fdata(end,2);
            obj.boreholes(end).Zmax = fdata(end,3);
            obj.boreholes(end).fdata = fdata;

            names = cell(1,length(obj.boreholes));
            for n=1:length(obj.boreholes)
                names{n} = obj.boreholes(n).name;
            end
            obj.handles.listBH.String = names;
            obj.handles.listBH.Value = length(obj.boreholes);
            obj.updateEdits
            notify(obj,'boreholeAdded')
        end
        function plotBH(obj, varargin)
            if isempty(obj.boreholes), return, end

            figure
            for n=1:length(obj.boreholes)
                plot3([obj.boreholes(n).X obj.boreholes(n).Xmax],...
                    [obj.boreholes(n).Y obj.boreholes(n).Ymax],...
                    [obj.boreholes(n).Z obj.boreholes(n).Zmax],'g','LineWidth',2)
                hold on
                plot3(obj.boreholes(n).X, obj.boreholes(n).Y, obj.boreholes(n).Z_surf,'ro')
                text(obj.boreholes(n).X, obj.boreholes(n).Y, ...
                    obj.boreholes(n).Z_surf, obj.boreholes(n).name)
            end
            grid on
            hold off
            xlabel('X [m]')
            ylabel('Y [m]')
            zlabel('Elevation [m]')
            set(gca,'DataAspectRatio',[1 1 1])

        end
        function listBH(obj,varargin)
            obj.updateEdits
        end
        function topX(obj, varargin)
            no = obj.handles.listBH.Value;
            if no>0 && no<=length(obj.boreholes)
                if obj.boreholes(no).X==obj.boreholes(no).fdata(1,1)
                    obj.boreholes(no).X = str2double( obj.handles.topX.String );
                    obj.boreholes(no).fdata(1,1) = obj.boreholes(no).X;
                else
                    % prepend fdata with new value
                    obj.boreholes(no).X = str2double( obj.handles.topX.String );
                    obj.boreholes(no).fdata = [obj.boreholes(no).X obj.boreholes(no).Y obj.boreholes(no).Z;
                        obj.boreholes(no).fdata];
                end
                notify(obj,'coordinatesChanged');
            end
        end
        function topY(obj, varargin)
            no = obj.handles.listBH.Value;
            if no>0 && no<=length(obj.boreholes)
                if obj.boreholes(no).Y==obj.boreholes(no).fdata(1,2)
                    obj.boreholes(no).Y = str2double( obj.handles.topY.String );
                    obj.boreholes(no).fdata(1,2) = obj.boreholes(no).Y;
                else
                    % prepend fdata with new value
                    obj.boreholes(no).Y = str2double( obj.handles.topY.String );
                    obj.boreholes(no).fdata = [obj.boreholes(no).X obj.boreholes(no).Y obj.boreholes(no).Z;
                        obj.boreholes(no).fdata];
                end
                notify(obj,'coordinatesChanged');
            end
        end
        function topZ(obj, varargin)
            no = obj.handles.listBH.Value;
            if no>0 && no<=length(obj.boreholes)
                if obj.boreholes(no).Z==obj.boreholes(no).fdata(1,3)
                    obj.boreholes(no).Z = str2double( obj.handles.topZ.String );
                    obj.boreholes(no).fdata(1,3) = obj.boreholes(no).Z;
                else
                    % prepend fdata with new value
                    obj.boreholes(no).Z = str2double( obj.handles.topZ.String );
                    obj.boreholes(no).fdata = [obj.boreholes(no).X obj.boreholes(no).Y obj.boreholes(no).Z;
                        obj.boreholes(no).fdata];                    
                end
                notify(obj,'coordinatesChanged');
            end

        end
        function bottomX(obj, varargin)
            no = obj.handles.listBH.Value;
            if no>0 && no<=length(obj.boreholes)
                if obj.boreholes(no).Xmax == obj.boreholes(no).fdata(end,1)
                    obj.boreholes(no).Xmax = str2double( obj.handles.bottomX.String );
                    obj.boreholes(no).fdata(end,1) = obj.boreholes(no).Xmax;
                else
                    % append new value to fdata
                    obj.boreholes(no).Xmax = str2double( obj.handles.bottomX.String );
                    obj.boreholes(no).fdata = [obj.boreholes(no).fdata;
                        obj.boreholes(no).X obj.boreholes(no).Y obj.boreholes(no).Z];
                end
                notify(obj,'coordinatesChanged');
            end
        end
        function bottomY(obj, varargin)
            no = obj.handles.listBH.Value;
            if no>0 && no<=length(obj.boreholes)
                if obj.boreholes(no).Ymax == obj.boreholes(no).fdata(end,2)
                    obj.boreholes(no).Ymax = str2double( obj.handles.bottomY.String );
                    obj.boreholes(no).fdata(end,2) = obj.boreholes(no).Ymax;
                else
                    % append new value to fdata
                    obj.boreholes(no).Ymax = str2double( obj.handles.bottomY.String );
                    obj.boreholes(no).fdata = [obj.boreholes(no).fdata;
                        obj.boreholes(no).X obj.boreholes(no).Y obj.boreholes(no).Z];
                end
                notify(obj,'coordinatesChanged');
            end
        end
        function bottomZ(obj, varargin)
            no = obj.handles.listBH.Value;
            if no>0 && no<=length(obj.boreholes)
                if obj.boreholes(no).Zmax == obj.boreholes(no).fdata(end,3)
                    obj.boreholes(no).Zmax = str2double( obj.handles.bottomZ.String );
                    obj.boreholes(no).fdata(end,3) = obj.boreholes(no).Zmax;
                else
                    % append new value to fdata
                    obj.boreholes(no).Zmax = str2double( obj.handles.bottomZ.String );
                    obj.boreholes(no).fdata = [obj.boreholes(no).fdata;
                        obj.boreholes(no).X obj.boreholes(no).Y obj.boreholes(no).Z];
                end
                notify(obj,'coordinatesChanged');
            end
        end
        function Zsurf(obj, varargin)
            no = obj.handles.listBH.Value;
            if no>0 && no<=length(obj.boreholes)
                obj.boreholes(no).Z_surf = str2double( obj.handles.Zsurf.String );
            end
        end
        function Zwater(obj, varargin)
            no = obj.handles.listBH.Value;
            if no>0 && no<=length(obj.boreholes)
                obj.boreholes(no).Z_water = str2double( obj.handles.Zwater.String );
            end
        end
        function diameter(obj, varargin)
            no = obj.handles.listBH.Value;
            if no>0 && no<=length(obj.boreholes)
                obj.boreholes(no).diam = str2double( obj.handles.diameter.String );
            end
        end
        function contS(obj, varargin)
            no = obj.handles.listBH.Value;
            if no>0 && no<=length(obj.boreholes)
                [file, rep] = uigetfile('*.con','Constraints file - Slowness');
                if file==0
                    return
                end
                cont = load([rep,file]);
                scont.x = obj.boreholes(no).X;    % vertical straight boreholes
                scont.y = obj.boreholes(no).Y;
                scont.z = obj.boreholes(no).Z - cont(:,1);
                scont.valeur = cont(:,2);
                if size(cont,2)==3
                    scont.variance = cont(:,3);
                else
                    scont.variance = zeros(size(cont(:,2)));
                end
                obj.boreholes(no).scont = scont;
            end
        end
        function contA(obj, varargin)
            no = obj.handles.listBH.Value;
            if no>0 && no<=length(obj.boreholes)
                [file, rep] = uigetfile('*.con','Constraints file - Attenuation');
                if file==0
                    return
                end
                cont = load([rep,file]);
                acont.x = obj.boreholes(no).X;    % vertical straight boreholes
                acont.y = obj.boreholes(no).Y;
                acont.z = obj.boreholes(no).Z - cont(:,1);
                acont.valeur = cont(:,2);
                if size(cont,2)==3
                    acont.variance = cont(:,3);
                else
                    acont.variance = zeros(size(cont(:,2)));
                end
                obj.boreholes(no).acont = acont;
            end
        end

        function updateEdits(obj)
            no = obj.handles.listBH.Value;
            if no>0 && no<=length(obj.boreholes)
                obj.handles.topX.String = num2str( obj.boreholes(no).X );
                obj.handles.topY.String = num2str( obj.boreholes(no).Y );
                obj.handles.topZ.String = num2str( obj.boreholes(no).Z );
                obj.handles.bottomX.String = num2str( obj.boreholes(no).Xmax );
                obj.handles.bottomY.String = num2str( obj.boreholes(no).Ymax );
                obj.handles.bottomZ.String = num2str( obj.boreholes(no).Zmax );
                obj.handles.Zsurf.String = num2str( obj.boreholes(no).Z_surf );
                obj.handles.Zwater.String = num2str( obj.boreholes(no).Z_water );
                obj.handles.diameter.String = num2str( obj.boreholes(no).diam );
            end
        end
    end
end
