classdef CovarianceUI < handle
    %COVARIANCEUI User interface for covariance models
    
    properties
        range
        angle   % angle in degrees
        sill
        fixRange
        fixAngle
        fixSill
    end
    properties (Access=private)
        handles
    end
    events
        covarianceEdited, typeChanged
    end
    
    methods
        function obj = CovarianceUI(r,a,s,varargin)
            if numel(r)==2 && numel(a)~=1
                error('angle must be a scalar for 2D media')
            elseif numel(r)==3 && numel(a)~=3
                error('3 angle values needed for 3D media')
            end
            
            obj.range = r;
            obj.angle = a;
            obj.sill = s;
            obj.fixRange = false(size(r));
            obj.fixAngle = false(size(a));
            obj.fixSill = false(size(s));
            
            obj.handles.hp = uipanel(varargin{:},...
                'Visible','off',...
                'BorderType','line');
            
            obj.addComponents();
            obj.handles.hp.Visible = 'on';
            obj.handles.hlist = [];
        end
        function delete(obj)
            delete(obj.handles.hp)
            delete(obj.handles.hlist)
        end
        function createList(obj,varargin)
            [~,ll]=enumeration('CovarianceModels');
            l=cell(numel(ll)-1,1); % we don't want nugget (last in enum list)
            for n=1:numel(ll)-1;
                l{n}=ll{n};
            end
            obj.handles.hlist = uicontrol('Style','popupmenu',...
                'String',l,...
                'callback',@obj.changeType,...
                varargin{:});
        end
        function t = getType(obj)
            if ~isempty(obj.handles.hlist)
                types = enumeration('CovarianceModels');
                t = types(obj.handles.hlist.Value);
            else
                t = CovarianceModels.empty;
            end
        end
        function setType(obj,t)
            obj.handles.hlist.Value = int8(t);
        end
        function setVisible(obj,v)
            obj.handles.hp.Visible = v;
            obj.handles.hlist.Visible = v;
        end
        function refresh(obj)
            obj.handles.editRangeX.String = num2str(obj.range(1));
            if length(obj.range)==3
                obj.handles.editRangeY.String = num2str(obj.range(2));
            end
            obj.handles.editRangeZ.String = num2str(obj.range(end));
            obj.handles.editAngleX.String = num2str(obj.angle(1));
            if length(obj.angle)==3
                obj.handles.editAngleY.String = num2str(obj.angle(2));
                obj.handles.editAngleZ.String = num2str(obj.angle(3));
            end
            obj.handles.editSill.String = num2str(obj.sill);
        end
    end
    
    methods (Access=private)
        function editRangeX(obj,varargin)
            obj.range(1) = str2double(obj.handles.editRangeX.String);
            obj.notify('covarianceEdited')
        end
        function editRangeY(obj,varargin)
            obj.range(2) = str2double(obj.handles.editRangeY.String);
            obj.notify('covarianceEdited')
        end
        function editRangeZ(obj,varargin)
            obj.range(end) = str2double(obj.handles.editRangeZ.String);
            obj.notify('covarianceEdited')
        end
        function editAngleX(obj,varargin)
            obj.angle(1) = str2double(obj.handles.editAngleX.String);
            obj.notify('covarianceEdited')
        end
        function editAngleY(obj,varargin)
            obj.angle(2) = str2double(obj.handles.editAngleY.String);
            obj.notify('covarianceEdited')
        end
        function editAngleZ(obj,varargin)
            obj.angle(3) = str2double(obj.handles.editAngleZ.String);
            obj.notify('covarianceEdited')
        end
        function editSill(obj,varargin)
            obj.sill = str2double(obj.handles.editSill.String);
            obj.notify('covarianceEdited')
        end
        
        function fixRangeX(obj,varargin)
            obj.fixRange(1) = obj.handles.fixRangeX.Value;
        end
        function fixRangeY(obj,varargin)
            obj.fixRange(2) = obj.handles.fixRangeY.Value;
        end
        function fixRangeZ(obj,varargin)
            obj.fixRange(end) = obj.handles.fixRangeZ.Value;
        end
        function fixAngleX(obj,varargin)
            obj.fixAngle(1) = obj.handles.fixAngleX.Value;
        end
        function fixAngleY(obj,varargin)
            obj.fixAngle(2) = obj.handles.fixAngleY.Value;
        end
        function fixAngleZ(obj,varargin)
            obj.fixAngle(3) = obj.handles.fixAngleZ.Value;
        end
        function fixSillVal(obj,varargin)
            obj.fixSill = obj.handles.fixSill.Value;
        end
        
        function addComponents(obj)
            d = numel(obj.range);
            
            nLines=4;
            if d==3
                nLines=7;
            end
            vSizeTot = nLines*22 + (nLines+1)*5;
            vSize = 22/vSizeTot;
            vSpace = 5/vSizeTot;
            hSize = 0.7;
            hSize2 = 0.21;
            hSpace = 0.03;
            
            obj.handles.editRangeX = uicontrol('Style','edit',...
                'String',num2str(obj.range(1)),...
                'Units','normalized',...
                'Callback',@obj.editRangeX,...
                'Parent',obj.handles.hp);
            obj.handles.editRangeZ = uicontrol('Style','edit',...
                'String',num2str(obj.range(end)),...
                'Units','normalized',...
                'Callback',@obj.editRangeZ,...
                'Parent',obj.handles.hp);
            obj.handles.editAngleX = uicontrol('Style','edit',...
                'String',num2str(obj.angle(1)),...
                'Units','normalized',...
                'Callback',@obj.editAngleX,...
                'Parent',obj.handles.hp);
            obj.handles.editSill = uicontrol('Style','edit',...
                'String',num2str(obj.sill),...
                'Units','normalized',...
                'Callback',@obj.editSill,...
                'Parent',obj.handles.hp);
            
            obj.handles.fixRangeX = uicontrol('Style','checkbox',...
                'Units','normalized',...
                'Callback',@obj.fixRangeX,...
                'Parent',obj.handles.hp);
            obj.handles.fixRangeZ = uicontrol('Style','checkbox',...
                'Units','normalized',...
                'Callback',@obj.fixRangeZ,...
                'Parent',obj.handles.hp);
            obj.handles.fixAngleX = uicontrol('Style','checkbox',...
                'Units','normalized',...
                'Callback',@obj.fixAngleX,...
                'Parent',obj.handles.hp);
            obj.handles.fixSill = uicontrol('Style','checkbox',...
                'Units','normalized',...
                'Callback',@obj.fixSillVal,...
                'Parent',obj.handles.hp);

            
            if d == 3
                obj.handles.editRangeY = uicontrol('Style','edit',...
                    'String',num2str(obj.range(2)),...
                    'Units','normalized',...
                    'Position',[hSpace 6*vSpace+5*vSize hSize vSize],...
                    'Callback',@obj.editRangeY,...
                    'Parent',obj.handles.hp);
                obj.handles.editAngleY = uicontrol('Style','edit',...
                    'String',num2str(obj.angle(2)),...
                    'Units','normalized',...
                    'Position',[hSpace 3*vSpace+2*vSize hSize vSize],...
                    'Callback',@obj.editAngleY,...
                    'Parent',obj.handles.hp);
                obj.handles.editAngleZ = uicontrol('Style','edit',...
                    'String',num2str(obj.angle(3)),...
                    'Units','normalized',...
                    'Position',[hSpace 2*vSpace+  vSize hSize vSize],...
                    'Callback',@obj.editAngleZ,...
                    'Parent',obj.handles.hp);
                
                obj.handles.fixRangeY = uicontrol('Style','checkbox',...
                    'Units','normalized',...
                    'Position',[2*hSpace+hSize 6*vSpace+5*vSize hSize2 vSize],...
                    'Callback',@obj.fixRangeY,...
                    'Parent',obj.handles.hp);
                obj.handles.fixAngleY = uicontrol('Style','checkbox',...
                    'Units','normalized',...
                    'Position',[2*hSpace+hSize 3*vSpace+2*vSize hSize2 vSize],...
                    'Callback',@obj.fixAngleY,...
                    'Parent',obj.handles.hp);
                obj.handles.fixAngleZ = uicontrol('Style','checkbox',...
                    'Units','normalized',...
                    'Position',[2*hSpace+hSize 2*vSpace+  vSize hSize2 vSize],...
                    'Callback',@obj.fixAngleZ,...
                    'Parent',obj.handles.hp);

                obj.handles.editRangeX.Position = [hSpace 7*vSpace+6*vSize hSize vSize];
                obj.handles.editRangeZ.Position = [hSpace 5*vSpace+4*vSize hSize vSize];
                obj.handles.editAngleX.Position = [hSpace 4*vSpace+3*vSize hSize vSize];
                obj.handles.editSill.Position =   [hSpace   vSpace         hSize vSize];

                obj.handles.fixRangeX.Position = [2*hSpace+hSize 7*vSpace+6*vSize hSize2 vSize];
                obj.handles.fixRangeZ.Position = [2*hSpace+hSize 5*vSpace+4*vSize hSize2 vSize];
                obj.handles.fixAngleX.Position = [2*hSpace+hSize 4*vSpace+3*vSize hSize2 vSize];
                obj.handles.fixSill.Position =   [2*hSpace+hSize   vSpace         hSize2 vSize];
            else
                obj.handles.editRangeX.Position = [hSpace 4*vSpace+3*vSize hSize vSize];
                obj.handles.editRangeZ.Position = [hSpace 3*vSpace+2*vSize hSize vSize];
                obj.handles.editAngleX.Position = [hSpace 2*vSpace+  vSize hSize vSize];
                obj.handles.editSill.Position =   [hSpace   vSpace         hSize vSize];

                obj.handles.fixRangeX.Position = [2*hSpace+hSize 4*vSpace+3*vSize hSize2 vSize];
                obj.handles.fixRangeZ.Position = [2*hSpace+hSize 3*vSpace+2*vSize hSize2 vSize];
                obj.handles.fixAngleX.Position = [2*hSpace+hSize 2*vSpace+  vSize hSize2 vSize];
                obj.handles.fixSill.Position =   [2*hSpace+hSize   vSpace         hSize2 vSize];
            end
        end
        function changeType(obj,varargin)
            obj.notify('typeChanged')
        end
    end
end

