classdef GridViewer < handle
    %UNTITLED Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        slider1
        slider2
        simuSlider
    end
    properties (Access = private, Hidden = true)
        x
        y
        z
        data
        simu
        varSimu
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
            title(axes,ti,'FontSize',14,'Interpreter','none')
            xlabel(axes,xl,'FontSize',12)
            if ~isempty(yl)
                ylabel(axes,yl,'FontSize',12)
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
                    'HorizontalAlignment','center',...
                    'Parent',axes);
            else
                obj.himsc2 = imagesc(obj.x,obj.z,sl,'Parent',axes);
                obj.htext2 = text(obj.x(1),obj.z(end)+obj.grid.dz,...
                    ['Y = ',num2str(obj.y(iy))],...
                    'FontSize',14,...
                    'HorizontalAlignment','center',...
                    'Parent',axes);
            end
            set(axes,'DataAspectRatio',[1 1 1],'YDir','normal')
            title(axes,ti,'FontSize',14,'Interpreter','none')
            xlabel(axes,xl,'FontSize',12)
            if ~isempty(yl)
                zlabel(axes,yl,'FontSize',12)
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
        function plotSimu(obj,simu,nf)
            nsimu = size(simu,2);
            v = var(simu,0,2);
            
            obj.simuSlider = uicontrol('Style','slider',...
                'Min',1,...
                'Max',nsimu,...
                'SliderStep',[1/nsimu 1/nsimu],...
                'Value',round(nsimu/2),...
                'Callback',@obj.updateSimu,...
                'Parent',nf);

            if strcmp(obj.grid.type,'3D')
                obj.varSimu = reshape(v,length(obj.x),length(obj.y),length(obj.z));
                obj.simu = reshape(simu,length(obj.x),length(obj.y),length(obj.z),nsimu);
                obj.plotSimu3D(nf)
                obj.slider1.Callback = @obj.updateSimu;
            else
                v = reshape(v,length(obj.z),length(obj.x));
                obj.simu = reshape(simu,length(obj.z),length(obj.x),nsimu);
                obj.plotSimu2D(nf,v)
            end
        end
        function plotSimu2D(obj,nf,v)
            ax = subplot(1,2,1,'Parent',nf);
            
            is = round(obj.simuSlider.Value);
            s = squeeze(obj.simu(:,:,is));
            obj.himsc1 = imagesc(obj.x,obj.z,s,'Parent',ax);
            set(ax,'DataAspectRatio',[1 1 1],'YDir','normal')
            title(ax,['Simulation no ',num2str(is)],'FontSize',14);
            colorbar('peer',ax)
            
            ax = subplot(1,2,2,'Parent',nf);
            
            imagesc(obj.x,obj.z,v,'Parent',ax);
            set(ax,'DataAspectRatio',[1 1 1],'YDir','normal')
            title(ax,'Variance of Simulations','FontSize',14);
            colorbar('peer',ax)
        end
        function plotSimu3D(obj,nf)
            ax = subplot(1,2,1,'Parent',nf);
            
            is = round(obj.simuSlider.Value);
            iy = round(obj.slider1.Value);
            s = squeeze(obj.simu(:,iy,:,is))';
            obj.himsc1 = imagesc(obj.x,obj.z,s,'Parent',ax);
            obj.htext1 = text(obj.x(1),obj.z(end)+obj.grid.dz,...
                ['Y = ',num2str(obj.y(iy))],...
                'FontSize',14,...
                'HorizontalAlignment','center');
            set(ax,'DataAspectRatio',[1 1 1],'YDir','normal')
            title(ax,['Simulation no ',num2str(is)],'FontSize',14);
            colorbar('peer',ax)
            
            ax = subplot(1,2,2,'Parent',nf);
            v = squeeze(obj.varSimu(:,iy,:))';
            obj.himsc2 = imagesc(obj.x,obj.z,v,'Parent',ax);
            obj.htext2 = text(obj.x(1),obj.z(end)+obj.grid.dz,...
                ['Y = ',num2str(obj.y(iy))],...
                'FontSize',14,...
                'HorizontalAlignment','center');
            set(ax,'DataAspectRatio',[1 1 1],'YDir','normal')
            title(ax,'Variance of Simulations','FontSize',14);
            colorbar('peer',ax)
        end
        function updateSimu(obj,varargin)
            is = round(obj.simuSlider.Value);
            if strcmp(obj.grid.type,'3D')
                iy = round(obj.slider1.Value);
                s = squeeze(obj.simu(:,iy,:,is))';
                obj.himsc1.CData = s;
                obj.htext1.String = ['Y = ',num2str(obj.y(iy))];
                v = squeeze(obj.varSimu(:,iy,:))';
                obj.himsc2.CData = v;
                obj.htext2.String = ['Y = ',num2str(obj.y(iy))];
            else
                s = squeeze(obj.simu(:,:,is));
                obj.himsc1.CData = s;
            end
            title(obj.himsc1.Parent,['Simulation no ',num2str(is)],'FontSize',14);
        end
    end
    
end

