classdef ModelUI < handle
    %MODELUI User interface to manage 2D and 3D models
    
    properties
        models
    end
    properties (Dependent)
        Position
        FontSize   % should be within [10 11 12]
    end
        properties (Access=private)
        handles
    end
    events
        modelAdded, modelDeleted
    end

    methods
        function obj = ModelUI(varargin)
            obj.handles.hp = uipanel(varargin{:},...
            'Title','Models',...
            'Visible','off',...
            'SizeChangedFcn', @obj.resizeUI);
            
            obj.addComponents();
            obj.resizeUI();
        end
        function set.models(obj,m)
            if isa(m, 'Model')
                obj.models = m;
            else
                error('Models should be objects of class Model')
            end
            obj.updateList(1)
        end
        function set.Position(obj, p)
            obj.handles.hp.Position = p;
        end
        function set.FontSize(obj, s)
            obj.handles.hp.FontSize = s+1;
            obj.handles.addModel.FontSize = s;
            obj.handles.removeModel.FontSize = s;
            obj.handles.listModels.FontSize = s;
            obj.handles.grid.FontSize = s;
            obj.handles.createGrid.FontSize = s;
            obj.handles.editGrid.FontSize = s;
            obj.handles.mogs.FontSize = s;
            obj.handles.addMOG.FontSize = s;
            obj.handles.removeMOG.FontSize = s;
            obj.handles.listMOGs.FontSize = s;
            
            obj.resizeUI();
        end
    end
   
    methods (Access=private)
        function addComponents(obj)
            
            obj.handles.addModel = uicontrol('Style','pushbutton',...
                'String','Add Model',...
                'Units','points',...
                'Callback',@obj.addModel,...
                'Parent',obj.handles.hp);
            obj.handles.removeModel = uicontrol('Style','pushbutton',...
                'String','Remove Model',...
                'Units','points',...
                'Callback',@obj.removeModel,...
                'Parent',obj.handles.hp);
            
            obj.handles.listModels = uicontrol('Style','listbox',...
                'Max',1,'Min',0,...
                'Units','points',...
                'Callback',@obj.listModels,...
                'Parent',obj.handles.hp);
            
            obj.handles.grid = uipanel(obj.handles.hp,'Title','Grid',...
                'Units','points',...
                'SizeChangedFcn', @obj.resizeUI);
            
            obj.handles.createGrid = uicontrol('Style','pushbutton',...
                'String','Create',...
                'Units','points',...
                'Callback',@obj.createGrid,...
                'Parent',obj.handles.grid);
            obj.handles.editGrid = uicontrol('Style','pushbutton',...
                'String','Edit',...
                'Units','points',...
                'Callback',@obj.createGrid,...
                'Parent',obj.handles.grid);
            
            obj.handles.mogs = uipanel(obj.handles.hp,'Title','MOGs',...
                'Units','points',...
                'SizeChangedFcn', @obj.resizeUI);
            obj.handles.addMOG = uicontrol('Style','pushbutton',...
                'String','Add MOG',...
                'Units','points',...
                'Callback',@obj.addMOG,...
                'Parent',obj.handles.mogs);
            obj.handles.removeMOG = uicontrol('Style','pushbutton',...
                'String','Remove MOG',...
                'Units','points',...
                'Callback',@obj.removeMOG,...
                'Parent',obj.handles.mogs);
            obj.handles.listMOGs = uicontrol('Style','listbox',...
                'Max',1,'Min',0,...
                'Units','points',...
                'Callback',@obj.listMOGs,...
                'Parent',obj.handles.mogs);

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

            obj.handles.addModel.Position = [hBorder height-vBorderTop 1.5*hSize vSize];
            obj.handles.removeModel.Position = [hBorder+2*hSpace+1.5*hSize height-vBorderTop 1.5*hSize vSize];
            
            obj.handles.listModels.Position = [1.5*hBorder 3*vSpace+vBorder+2.5*vSize 2*hSpace+3*hSize-hBorder height-2*vBorder-vBorderTop-3.5*vSize];
            
            hSize2 = 2*hSpace+3*hSize;
            obj.handles.grid.Position = [hBorder vBorder hSize2 2.5*vSize];
            obj.handles.createGrid.Position = [hBorder vSize/2 hSize2/2-3/2*hBorder vSize];
            obj.handles.editGrid.Position = [hSize2/2+hBorder/2 vSize/2 hSize2/2-3/2*hBorder vSize];
            
            hSize2 = width-hSize2-3*hBorder;
            vSize2 = height-3*vBorder;
            obj.handles.mogs.Position = [2*hBorder+2*hSpace+3*hSize vBorder hSize2 height-3*vBorder];

            obj.handles.addMOG.Position = [hBorder vSize2-vBorderTop hSize2/2-3/2*hBorder vSize];
            obj.handles.removeMOG.Position = [hSize2/2+hBorder/2 vSize2-vBorderTop hSize2/2-3/2*hBorder vSize];
            obj.handles.listMOGs.Position = [1.5*hBorder vBorder hSize2-3*hBorder vSize2-2*vBorder-vBorderTop];
            
            obj.handles.hp.Visible = 'on';
            obj.handles.hp.Units = oldUnits;
        end
        
        function addModel(obj,varargin)
            name = myinputdlg('Model name');
            if isempty(name)
                return;
            end
            if isempty(char(name))
                return;
            end
            if isempty(obj.models)
                obj.models = Model(char(name));
            else
                obj.models(end+1) = Model(char(name));
            end
            obj.updateList()
            obj.notify('modelAdded')
        end
        function removeModel(obj,varargin)
            no = obj.handles.listModels.Value;
            ind=1:length(obj.models);
            ind = ind~=no;

            obj.models = obj.models(ind);
            obj.updateList()
            obj.notify('modelDeleted')
        end
        function createGrid(obj,varargin)
        end
        function editGrid(obj,varargin)
        end
        function addMOG(obj,varargin)
        end
        function removeMOG(obj,varargin)
        end
        function listModels(obj,varargin)
        end
        function listMOGs(obj,varargin)
        end
        
        function updateList(obj,varargin)
            value = [];
            if nargin==2
                value = varargin{1};
            end
            if ~isempty(obj.models)
                names = cell(1,length(obj.models));
                for n=1:length(obj.models)
                    names{n} = obj.models(n).name;
                end
                obj.handles.listModels.String = names;
                obj.handles.listModels.Value = length(obj.models);
            else
                obj.handles.listModels.String = '';
                obj.handles.listModels.Value = 1;
            end
            if ~isempty(value)
                obj.handles.listModels.Value = value;
            end
        end
        
    end

end

