classdef GridViewer < handle
    %UNTITLED Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        slider1
        slider2
    end
    properties (Access = private, Hidden = true)
        x
        y
        z
        data
        grid
        himsc1
        himsc2
        htext1
        htext2
    end
    
    methods
        function obj = GridViewer(g)
            obj.grid = g;
            obj.x=(g.grx(1)+g.dx/2):g.dx:(g.grx(end)-g.dx/3);
            obj.z=(g.grz(1)+g.dz/2):g.dz:(g.grz(end)-g.dz/3);
            if strcmp(g.type, '3D')
                obj.y=(g.gry(1)+g.dy/2):g.dy:(g.gry(end)-g.dy/3);
            end
        end
        function createSliders(obj,varargin)
            ny = length(obj.y);
            obj.slider1 = uicontrol('Style','slider',...
                'Min',1,...
                'Max',ny,...
                'SliderStep',[1/ny 1/ny],...
                'Value',round(ny/2),...
                'Callback',@obj.update3D,...
                varargin{:});
            obj.slider2 = uicontrol('Style','slider',...
                'Min',1,...
                'Max',ny,...
                'SliderStep',[1/ny 1/ny],...
                'Value',round(ny/2),...
                'Callback',@obj.update3D,...
                varargin{:});
        end
        
        function plotTomo(obj,tomo,ti,xl,yl,axes,varargin)
            if strcmp(obj.grid.type,'3D')
                obj.plotTomo3D(tomo,ti,xl,yl,axes,varargin{:})
            else
                obj.plotTomo2D(tomo,ti,xl,yl,axes)
            end
        end
        function plotTomo2D(obj,tomo,ti,xl,yl,axes)
            imagesc(obj.x,obj.z,reshape(tomo,length(obj.z),length(obj.x)),'Parent',axes)
            set(axes,'DataAspectRatio',[1 1 1],'YDir','normal')
            title(axes,ti,'FontSize',14)
            xlabel(axes,xl)
            if ~isempty(yl)
                ylabel(axes,yl)
            end
        end
        function plotTomo3D(obj,tomo,ti,xl,yl,axes,varargin)
            no=1;
            if nargin>=7
                no = varargin{1};
            end
            
            obj.data(no,:,:,:) = reshape(tomo,length(obj.x),length(obj.y),length(obj.z));
            if no==1
                iy = round(obj.slider1.Value);
            else
                iy = round(obj.slider2.Value);
            end
            sl = squeeze(obj.data(no,:,iy,:))';
            
            if no==1
                obj.himsc1 = imagesc(obj.x,obj.z,sl,'Parent',axes);
                obj.htext1 = text(obj.x(1),obj.z(end)+obj.grid.dz,...
                    ['Y = ',num2str(obj.y(iy))],...
                    'FontSize',14,...
                    'HorizontalAlignment','center');
            else
                obj.himsc2 = imagesc(obj.x,obj.z,sl,'Parent',axes);
                obj.htext1 = text(obj.x(1),obj.z(end)+obj.grid.dz,...
                    ['Y = ',num2str(obj.y(iy))],...
                    'FontSize',14,...
                    'HorizontalAlignment','center');
            end
            set(axes,'DataAspectRatio',[1 1 1],'YDir','normal')
            title(axes,ti,'FontSize',14)
            xlabel(axes,xl)
            if ~isempty(yl)
                ylabel(axes,yl)
            end
        end
        function update3D(obj,src,varargin)
            if isempty(obj.data)
                return
            end
            if src==obj.slider1
                no=1;
                iy = round(obj.slider1.Value);
                sl = squeeze(obj.data(no,:,iy,:))';
                obj.himsc1.CData = sl;
                obj.htext1.String = ['Y = ',num2str(obj.y(iy))];
            elseif src==obj.slider2
                no=2;
                iy = round(obj.slider2.Value);
                sl = squeeze(obj.data(no,:,iy,:))';
                obj.himsc2.CData = sl;
                obj.htext2.String = ['Y = ',num2str(obj.y(iy))];
            end
        end
    end
    
end

