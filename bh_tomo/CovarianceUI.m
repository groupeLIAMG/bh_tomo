classdef CovarianceUI < handle
    %COVARIANCEUI User interface for covariance models
    
    properties (Access=private)
        cov
        guiHandles
    end
    
    methods
        function obj = CovarianceUI(c)
            obj.cov = c;
        end
        
        function createUI(obj, parent, position)
            obj.guiHandles.hp = uipanel('Parent',parent,...
                'Position',position,...
                'Title','Covariance parameters',...
                'Resize','on',...
                'SizeChangedFcn',@obj.resizeUI);
            
            obj.guiHandles.textRangeX = uicontrol('Style','text',...
                'String','Range X',...
                'Units','points',...
                'HorizontalAlignment','right',...
                'Parent',obj.guiHandles.hp);
            obj.guiHandles.editRangeX = uicontrol('Style','edit',...
                'String',num2str(obj.cov.range(1)),...
                'Units','points',...
                'Callback',@obj.editRangeX,...
                'Parent',obj.guiHandles.hp);
            
            obj.guiHandles.textRangeZ = uicontrol('Style','text',...
                'String','Range Z',...
                'Units','points',...
                'HorizontalAlignment','right',...
                'Parent',obj.guiHandles.hp);
            obj.guiHandles.editRangeZ = uicontrol('Style','edit',...
                'String',num2str(obj.cov.range(end)),...
                'Units','points',...
                'Callback',@obj.editRangeZ,...
                'Parent',obj.guiHandles.hp);
            
            obj.guiHandles.textAngleX = uicontrol('Style','text',...
                'String','Angle',...
                'Units','points',...
                'HorizontalAlignment','right',...
                'Parent',obj.guiHandles.hp);
            obj.guiHandles.editAngleX = uicontrol('Style','edit',...
                'String',num2str(obj.cov.angle(1)),...
                'Units','points',...
                'Callback',@obj.editAngleX,...
                'Parent',obj.guiHandles.hp);

            obj.guiHandles.textSill = uicontrol('Style','text',...
                'String','Sill',...
                'Units','points',...
                'HorizontalAlignment','right',...
                'Parent',obj.guiHandles.hp);
            obj.guiHandles.editSill = uicontrol('Style','edit',...
                'String',num2str(obj.cov.sill),...
                'Units','points',...
                'Callback',@obj.editSill,...
                'Parent',obj.guiHandles.hp);

            d = numel(obj.cov.range);
            if d == 3
                obj.guiHandles.textRangeY = uicontrol('Style','text',...
                    'String','Range Y',...
                    'Units','points',...
                    'HorizontalAlignment','right',...
                    'Parent',obj.guiHandles.hp);
                obj.guiHandles.editRangeY = uicontrol('Style','edit',...
                    'String',num2str(obj.cov.range(2)),...
                    'Units','points',...
                    'Callback',@obj.editRangeY,...
                    'Parent',obj.guiHandles.hp);
                
                obj.guiHandles.textAngleX.String = 'Angle X';
                
                obj.guiHandles.textAngleY = uicontrol('Style','text',...
                    'String','Angle Y',...
                    'Units','points',...
                    'HorizontalAlignment','right',...
                    'Parent',obj.guiHandles.hp);
                obj.guiHandles.editAngleY = uicontrol('Style','edit',...
                    'String',num2str(obj.cov.angle(2)),...
                    'Units','points',...
                    'Callback',@obj.editAngleY,...
                    'Parent',obj.guiHandles.hp);

                obj.guiHandles.textAngleZ = uicontrol('Style','text',...
                    'String','Angle Z',...
                    'Units','points',...
                    'HorizontalAlignment','right',...
                    'Parent',obj.guiHandles.hp);
                obj.guiHandles.editAngleZ = uicontrol('Style','edit',...
                    'String',num2str(obj.cov.angle(3)),...
                    'Units','points',...
                    'Callback',@obj.editAngleZ,...
                    'Parent',obj.guiHandles.hp);
            end
        end
        function hp = getUI(obj, parent, position)
            obj.createUI(parent, position)
            hp = obj.guiHandles.hp;
        end
    end
    
    methods (Access=private)
        function editRangeX(obj,varargin)
            obj.cov.range(1) = str2double(obj.guiHandles.editRangeX.String);
        end
        function editRangeY(obj,varargin)
            obj.cov.range(2) = str2double(obj.guiHandles.editRangeY.String);
        end
        function editRangeZ(obj,varargin)
            obj.cov.range(end) = str2double(obj.guiHandles.editRangeZ.String);
        end
        function editAngleX(obj,varargin)
            obj.cov.angle(1) = str2double(obj.guiHandles.editAngleX.String);
        end
        function editAngleY(obj,varargin)
            obj.cov.angle(2) = str2double(obj.guiHandles.editAngleY.String);
        end
        function editAngleZ(obj,varargin)
            obj.cov.angle(3) = str2double(obj.guiHandles.editAngleZ.String);
        end
        function editSill(obj,varargin)
            obj.cov.sill = str2double(obj.guiHandles.editSill.String);
        end
        function resizeUI(obj,varargin)
            
            if ~isempty(obj.guiHandles)

                obj.guiHandles.hp.Units = 'points';
            
                % Get figure width and height
                pwidth = obj.guiHandles.hp.Position(3);
                pheight = obj.guiHandles.hp.Position(4);
                
                d = numel(obj.cov.range);
                if d==2
                    numLines = 4;
                else % d == 3
                    numLines = 7;
                end
                
                fac = 3;
                vSizeAll = (numLines-1)*obj.guiHandles.editRangeX.Position(4)/fac + ...
                    (numLines+1.5)*obj.guiHandles.editRangeX.Position(4);
                vSpace = obj.guiHandles.editRangeX.Position(4)/fac;
                vSize = obj.guiHandles.editRangeX.Position(4);
                if vSizeAll>pheight
                    vSize = fac*pheight / (numLines-1+fac*(numLines+1.5));
                    obj.guiHandles.textRangeX.Position(4) = vSize;
                    obj.guiHandles.textRangeZ.Position(4) = vSize;
                    obj.guiHandles.textAngleX.Position(4) = vSize;
                    obj.guiHandles.textSill.Position(4) = vSize;
                    obj.guiHandles.editRangeX.Position(4) = vSize;
                    obj.guiHandles.editRangeZ.Position(4) = vSize;
                    obj.guiHandles.editAngleX.Position(4) = vSize;
                    obj.guiHandles.editSill.Position(4) = vSize;
                    if d==3
                        obj.guiHandles.textRangeY.Position(4) = vSize;
                        obj.guiHandles.textAngleY.Position(4) = vSize;
                        obj.guiHandles.textAngleZ.Position(4) = vSize;
                        obj.guiHandles.editRangeY.Position(4) = vSize;
                        obj.guiHandles.editAngleY.Position(4) = vSize;
                        obj.guiHandles.editAngleZ.Position(4) = vSize;
                    end
                    vSpace = vSize/fac;
                end
                n=2.5;
                obj.guiHandles.textRangeX.Position(2) = pheight-n*vSize;
                obj.guiHandles.editRangeX.Position(2) = pheight-n*vSize;
                if d==3
                    n = n+1;
                    obj.guiHandles.textRangeY.Position(2) = pheight-n*vSize-(n-2.5)*vSpace;
                    obj.guiHandles.editRangeY.Position(2) = pheight-n*vSize-(n-2.5)*vSpace;
                end
                n = n+1;
                obj.guiHandles.textRangeZ.Position(2) = pheight-n*vSize-(n-2.5)*vSpace;
                obj.guiHandles.editRangeZ.Position(2) = pheight-n*vSize-(n-2.5)*vSpace;
                n = n+1;
                obj.guiHandles.textAngleX.Position(2) = pheight-n*vSize-(n-2.5)*vSpace;
                obj.guiHandles.editAngleX.Position(2) = pheight-n*vSize-(n-2.5)*vSpace;
                if d==3
                    n = n+1;
                    obj.guiHandles.textAngleY.Position(2) = pheight-n*vSize-(n-2.5)*vSpace;
                    obj.guiHandles.editAngleY.Position(2) = pheight-n*vSize-(n-2.5)*vSpace;
                    n = n+1;
                    obj.guiHandles.textAngleZ.Position(2) = pheight-n*vSize-(n-2.5)*vSpace;
                    obj.guiHandles.editAngleZ.Position(2) = pheight-n*vSize-(n-2.5)*vSpace;
                end                    
                n = n+1;
                obj.guiHandles.textSill.Position(2) = pheight-n*vSize-(n-2.5)*vSpace;
                obj.guiHandles.editSill.Position(2) = pheight-n*vSize-(n-2.5)*vSpace;

                hSpace = obj.guiHandles.textRangeX.Position(1)/2;
                hSize = obj.guiHandles.textRangeX.Position(3);
                hSizeAll = 5*hSpace + 2*hSize;
                if hSizeAll>pwidth
                    ratio = pwidth/hSizeAll;
                    hSpace = hSpace*ratio;
                    hSize = hSize*ratio;

                    obj.guiHandles.textRangeX.Position(3) = hSize;
                    obj.guiHandles.textRangeZ.Position(3) = hSize;
                    obj.guiHandles.textAngleX.Position(3) = hSize;
                    obj.guiHandles.textSill.Position(3) = hSize;
                    obj.guiHandles.editRangeX.Position(3) = hSize;
                    obj.guiHandles.editRangeZ.Position(3) = hSize;
                    obj.guiHandles.editAngleX.Position(3) = hSize;
                    obj.guiHandles.editSill.Position(3) = hSize;

                    obj.guiHandles.textRangeX.Position(1) = 2*hSpace;
                    obj.guiHandles.textRangeZ.Position(1) = 2*hSpace;
                    obj.guiHandles.textAngleX.Position(1) = 2*hSpace;
                    obj.guiHandles.textSill.Position(1) = 2*hSpace;
                    if d==3
                        obj.guiHandles.textRangeY.Position(3) = hSize;
                        obj.guiHandles.textAngleY.Position(3) = hSize;
                        obj.guiHandles.textAngleZ.Position(3) = hSize;
                        obj.guiHandles.editRangeY.Position(3) = hSize;
                        obj.guiHandles.editAngleY.Position(3) = hSize;
                        obj.guiHandles.editAngleZ.Position(3) = hSize;

                        obj.guiHandles.textRangeY.Position(1) = 2*hSpace;
                        obj.guiHandles.textAngleY.Position(1) = 2*hSpace;
                        obj.guiHandles.textAngleZ.Position(1) = 2*hSpace;
                        obj.guiHandles.editRangeY.Position(1) = 2*hSpace;
                        obj.guiHandles.editAngleY.Position(1) = 2*hSpace;
                        obj.guiHandles.editAngleZ.Position(1) = 2*hSpace;
                    end
                end
                obj.guiHandles.editRangeX.Position(1) = 3*hSpace+hSize;
                obj.guiHandles.editRangeZ.Position(1) = 3*hSpace+hSize;
                obj.guiHandles.editAngleX.Position(1) = 3*hSpace+hSize;
                obj.guiHandles.editSill.Position(1) = 3*hSpace+hSize;
                if d==3
                    obj.guiHandles.editRangeY.Position(1) = 3*hSpace+hSize;
                    obj.guiHandles.editAngleY.Position(1) = 3*hSpace+hSize;
                    obj.guiHandles.editAngleZ.Position(1) = 3*hSpace+hSize;
                end                    
            end
            
        end
    end
end

