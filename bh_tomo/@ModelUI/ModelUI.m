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
    properties (SetAccess=immutable,Hidden=true)
        bhUI
        mogUI
    end
    events
        modelAdded, modelDeleted, modelEdited
    end

    methods
        function obj = ModelUI(b,m,varargin)
            if isa(b,'BoreholeUI')
                obj.bhUI = b;
            else
                error('Give BoreholeUI object as input')
            end
            if isa(m,'MogUI')
                obj.mogUI = m;
            else
                error('Give BoreholeUI object as input')
            end
            obj.models = Model.empty;
            obj.handles.hp = uipanel(varargin{:},...
            'Title','Models',...
            'Visible','off',...
            'SizeChangedFcn', @obj.resizeUI);
            
            obj.addComponents();
            obj.handles.hp.Visible = 'on';
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
        end
    end
   
    methods (Access=private)
        g = gridEditor(obj,no,varargin)
        c = constraintsEditor(obj,varargin)
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
                'Units','points');
            
            obj.handles.createGrid = uicontrol('Style','pushbutton',...
                'String','Create',...
                'Units','points',...
                'Callback',@obj.createGrid,...
                'Parent',obj.handles.grid);
            obj.handles.editGrid = uicontrol('Style','pushbutton',...
                'String','Edit',...
                'Units','points',...
                'Callback',@obj.editGrid,...
                'Parent',obj.handles.grid);
            
            obj.handles.mogs = uipanel(obj.handles.hp,'Title','MOGs',...
                'Units','points');
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
                'Parent',obj.handles.mogs);

        end
        function resizeUI(obj,varargin)
            if isempty(obj.handles.hp)
                return
            end
            obj.handles.hp.Visible = 'off';
            oldUnits = obj.handles.hp.Units;
            obj.handles.hp.Units = 'points';

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
            if ispc
                vFac = 0.81*vFac;
            end
            vSize = 22*vFac;
            vSpace = 5*vFac;
            vBorderTop = 45*vFac;
            vBorder = 15*vFac;

            obj.handles.addModel.Position = [hBorder height-1.2*vBorderTop 1.5*hSize vSize];
            obj.handles.removeModel.Position = [hBorder+2*hSpace+1.5*hSize height-1.2*vBorderTop 1.5*hSize vSize];
            
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
            no = obj.handles.listModels.Value;
            if no>0 && no<=numel(obj.models)
                g = obj.gridEditor(no);
                if ~isempty(g)
                    obj.models(no).grid = g;
                end
            end
        end
        function editGrid(obj,varargin)
            no = obj.handles.listModels.Value;
            if no>0 && no<=numel(obj.models)
                g = obj.gridEditor(no,obj.models(no).grid);
                if ~isempty(g)
                    obj.models(no).grid = g;
                end
            end
        end
        function addMOG(obj,varargin)
            no = obj.handles.listModels.Value;
            if no>0 && no<=numel(obj.models)
                no_mog = chooseMOG(obj.mogUI);
                if no_mog==0
                    return
                end
                for n=1:numel(obj.models(no).mogs)
                    if obj.mogUI.mogs(no_mog).ID == obj.mogUI.mogs(obj.models(no).mogs(n)).ID 
                        warndlg('MOG already added to model')
                        return
                    end
                end
                obj.models(no).mogs(end+1) = no_mog;
                obj.models(no).boreholes = unique([obj.models(no).boreholes ...
                    obj.mogUI.mogs(no_mog).Tx obj.mogUI.mogs(no_mog).Rx]);
                obj.updateListMog()
                obj.notify('modelEdited')
            end            
        end
        function removeMOG(obj,varargin)
            no = obj.handles.listModels.Value;
            if no>0 && no<=numel(obj.models)
                no_mog = obj.handles.listMOGs.Value;
                ind=1:length(obj.models(no).mogs);
                ind = ind~=no_mog;
                obj.models(no).mogs = obj.models(no).mogs(ind);
                obj.handles.listMOGs.Value = 1;
                obj.updateListMog()
                obj.notify('modelEdited')
            end
        end
        function listModels(obj,varargin)
            obj.handles.listMOGs.Value = 1;
            obj.updateListMog();
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
            obj.updateListMog()
        end
        function updateListMog(obj,varargin)
            no = obj.handles.listModels.Value;
            if no>0 && no<=numel(obj.models)
                names_mog = cell(1,length(obj.models(no).mogs));
                for n=1:length(obj.models(no).mogs)
                    if isempty(obj.mogUI.mogs( obj.models(no).mogs(n) ).date)
                        names_mog{n} = obj.mogUI.mogs( obj.models(no).mogs(n) ).name;
                    else
                        names_mog{n} = [obj.mogUI.mogs( obj.models(no).mogs(n) ).name,' - ',...
                            obj.mogUI.mogs( obj.models(no).mogs(n) ).date];
                    end
                end
                obj.handles.listMOGs.String = names_mog;
            end
        end
        
    end

end

