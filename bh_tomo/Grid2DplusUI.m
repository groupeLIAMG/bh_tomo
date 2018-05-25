classdef Grid2DplusUI <handle
    %UNTITLED2 Summary of this class goes here
    %   Detailed explanation goes here
    
    properties (Dependent)
        Position
    end
    properties (SetAccess=private)
        type
        grid
        dx
        dz
        axPanel
    end
    properties (Access=private)
        handles
        data
        haxes
        
        order_boreholes
        n_planes
        order_planes
        xshift
    end
    events
        gridEdited
    end
    
    methods
        function obj = Grid2DplusUI(f,type,fs,g,data)
            
            obj.type = type;
            if isa(g,'Grid2D')
                obj.grid = g;
            else
                error('grid not an instance of Grid2D')
            end
            obj.data = data;
            
            if isempty(g.grx) % grid not initialized yet
                obj.dx = 1;
                obj.dz = 1;
            else
                obj.dx = g.grx(2)-g.grx(1);
                obj.dz = g.grz(2)-g.grz(1);
            end
            
            obj.addComponents(f,fs)
            obj.handles.params.Visible = 'on';
            obj.axPanel.Visible = 'on';
            
            obj.updateProj()
            
        end
        function set.Position(obj,p)
            obj.handles.params.Position = p;
        end
        function set.dx(obj,dx)
            if isnumeric(dx)
                obj.dx = dx;
            else
                error('Cell size should be numeric')
            end
        end
        function set.dz(obj,dz)
            if isnumeric(dz)
                obj.dz = dz;
            else
                error('Cell size should be numeric')
            end
        end
        function [xyz_p,dist] = project(obj,xyz)
            [xyz_p, xyz_no_plan] = Grid.proj_planes(xyz, obj.data.planes);
            dist = sqrt(sum((xyz-xyz_p).^2,2));
            [az,dip] = obj.get_azimuth_dip();
            for n=1:obj.data.n_planes
                no_f = obj.data.order_boreholes(n);
                no_plan = obj.data.order_planes(n);
                ind = xyz_no_plan==no_plan;
                oo = [obj.data.boreholes(no_f).X obj.data.boreholes(no_f).Y 0];
                xyz_p(ind,:) = Grid.transl_rotat(xyz_p(ind,:), oo, az(no_plan), dip(no_plan));
                if n>1
                    xyz_p(ind,1) = xyz_p(ind,1)+obj.data.planes(n-1).l;
                end
            end
            xyz_p(:,1) = xyz_p(:,1)-obj.xshift;
        end
    end
    
    methods (Static)
        function [grid,data] = buildGrid(data)
            grid = Grid2D();
            grid.bord = [1 1 1 1];
            grid.flip = 0;
            grid.type = '2D+';
            grid.nthreads = data.nthreads;
            data.order_boreholes = Grid.boreholes_order( data.boreholes );
            data.n_planes = length(data.boreholes)-1;
            data.order_planes = 1:data.n_planes;
            
            n1 = data.order_boreholes(1);
            grid.borehole_x0 = 1;
            grid.x0 = [data.boreholes(n1).X data.boreholes(n1).Y data.boreholes(n1).Z];
            for n=1:data.n_planes
                n1 = data.order_boreholes(n);
                n2 = data.order_boreholes(n+1);
                np = data.order_planes(n);
                X = [data.boreholes(n1).X   data.boreholes(n1).Y    data.boreholes(n1).Z;
                    data.boreholes(n1).Xmax data.boreholes(n1).Ymax data.boreholes(n1).Zmax;
                    data.boreholes(n2).X    data.boreholes(n2).Y    data.boreholes(n2).Z;
                    data.boreholes(n2).Xmax data.boreholes(n2).Ymax data.boreholes(n2).Zmax];
                % trouve le plan qui passe par les 4 points
                [data.planes(np).x0, data.planes(np).a]=Grid.lsplane(X);
                % data.planes(np).x0 : Centroid of the data = point on the best-fit plane
                % data.planes(np).a  : Direction cosines of the normal to the best-fit plane
                if data.planes(np).a(3)<0, data.planes(np).a=-data.planes(np).a; end
                % distance horizontale entre les trous
                data.planes(np).l = sqrt( (data.boreholes(n2).X-data.boreholes(n1).X)^2 +...
                    (data.boreholes(n2).Y-data.boreholes(n1).Y)^2 );
            end
            [~, data.Tx_no_plan] = Grid.proj_planes(data.Tx, data.planes);
            [~, data.Rx_no_plan] = Grid.proj_planes(data.Rx, data.planes);
        end
    end
    
    methods (Access=private)
        function addComponents(obj,f,fs)
            obj.handles.params = uipanel(f,'Title','Grid Parameters',...
                'FontSize',fs+1,...
                'Units','points',...
                'Visible','off',...
                'SizeChangedFcn',@obj.resizeUI);
            obj.handles.x = uicontrol('Style','text',...
                'String','X',...
                'FontSize',fs,...
                'Units','points',...
                'HorizontalAlignment','right',...
                'Parent',obj.handles.params);
            obj.handles.z = uicontrol('Style','text',...
                'String','Z',...
                'FontSize',fs,...
                'Units','points',...
                'HorizontalAlignment','right',...
                'Parent',obj.handles.params);
            obj.handles.padm = uicontrol('Style','text',...
                'String','Padding -',...
                'FontSize',fs,...
                'Units','points',...
                'HorizontalAlignment','center',...
                'Parent',obj.handles.params);
            obj.handles.padp = uicontrol('Style','text',...
                'String','Paddin +',...
                'FontSize',fs,...
                'Units','points',...
                'HorizontalAlignment','center',...
                'Parent',obj.handles.params);
            obj.handles.step = uicontrol('Style','text',...
                'String','Cell Size',...
                'FontSize',fs,...
                'Units','points',...
                'HorizontalAlignment','center',...
                'Parent',obj.handles.params);
            obj.handles.nxm = uicontrol('Style','edit',...
                'String',num2str(obj.grid.bord(1)),...
                'FontSize',fs,...
                'Units','points',...
                'Callback',@obj.inputEdit,...
                'Parent',obj.handles.params);
            obj.handles.nxp = uicontrol('Style','edit',...
                'String',num2str(obj.grid.bord(2)),...
                'FontSize',fs,...
                'Units','points',...
                'Callback',@obj.inputEdit,...
                'Parent',obj.handles.params);
            obj.handles.nzm = uicontrol('Style','edit',...
                'String',num2str(obj.grid.bord(3)),...
                'FontSize',fs,...
                'Units','points',...
                'Callback',@obj.inputEdit,...
                'Parent',obj.handles.params);
            obj.handles.nzp = uicontrol('Style','edit',...
                'String',num2str(obj.grid.bord(4)),...
                'FontSize',fs,...
                'Units','points',...
                'Callback',@obj.inputEdit,...
                'Parent',obj.handles.params);
            obj.handles.dx = uicontrol('Style','edit',...
                'String',num2str(obj.dx),...
                'FontSize',fs,...
                'Units','points',...
                'Callback',@obj.inputEdit,...
                'Parent',obj.handles.params);
            obj.handles.dz = uicontrol('Style','edit',...
                'String',num2str(obj.dz),...
                'FontSize',fs,...
                'Units','points',...
                'Callback',@obj.inputEdit,...
                'Parent',obj.handles.params);
            obj.handles.bh_x0 = uicontrol('Style','text',...
                'String','Borehole at Origin',...
                'FontSize',fs,...
                'Units','points',...
                'HorizontalAlignment','right',...
                'Parent',obj.handles.params);
            
            bhNames = cell(2,1);
            bhNames{1} = obj.data.boreholes(obj.data.order_boreholes(1)).name;
            bhNames{2} = obj.data.boreholes(obj.data.order_boreholes(end)).name;
            obj.handles.bhList = uicontrol('Style','popupmenu',...
                'String',bhNames,...
                'Value',1,...
                'FontSize',fs,...
                'Units','points',...
                'Callback',@obj.inputEdit,...
                'Parent',obj.handles.params);
            obj.handles.origin = uicontrol('Style','text',...
                'String','Origin',...
                'FontSize',fs,...
                'Units','points',...
                'HorizontalAlignment','right',...
                'Parent',obj.handles.params);
            obj.handles.x0x = uicontrol('Style','text',...
                'String',num2str(obj.grid.x0(1)),...
                'FontSize',fs,...
                'Units','points',...
                'HorizontalAlignment','center',...
                'BackgroundColor',[1 1 1],...
                'Parent',obj.handles.params);
            obj.handles.x0y = uicontrol('Style','text',...
                'String',num2str(obj.grid.x0(2)),...
                'FontSize',fs,...
                'Units','points',...
                'HorizontalAlignment','center',...
                'BackgroundColor',[1 1 1],...
                'Parent',obj.handles.params);
            obj.handles.x0z = uicontrol('Style','text',...
                'String',num2str(obj.grid.x0(3)),...
                'FontSize',fs,...
                'Units','points',...
                'HorizontalAlignment','center',...
                'BackgroundColor',[1 1 1],...
                'Parent',obj.handles.params);
            
            obj.handles.flip = uicontrol('Style','checkbox',...
                'String','Flip Horizontally',...
                'Value',obj.grid.flip,...
                'FontSize',fs,...
                'Units','points',...
                'Callback',@obj.inputEdit,...
                'Parent',obj.handles.params);
            
            obj.handles.numCores = uicontrol('Style','text',...
                'String','Nb of Cores',...
                'FontSize',fs,...
                'Units','points',...
                'HorizontalAlignment','right',...
                'Parent',obj.handles.params);
            
            ncores = feature('NumCores');
            cnc = cell(1, ncores);
            for n=1:ncores
                cnc{n} = sprintf('%d', n);
            end
            obj.handles.ncoreList = uicontrol('Style','popupmenu',...
                'String',cnc,...
                'Value',obj.grid.nthreads,...
                'FontSize',fs,...
                'Units','points',...
                'Callback',@obj.inputEdit,...
                'Parent',obj.handles.params);
            
            obj.handles.showFit = uicontrol('Style','pushbutton',...
                'String','Adjustment of Best-Fit Plane',...
                'FontSize',fs,...
                'Units','points',...
                'Callback',@obj.showFit,...
                'Parent',obj.handles.params);
            
            obj.axPanel = uipanel(f,'Title','Grid View',...
                'FontSize',fs+1,...
                'Units','points',...
                'Visible','off',...
                'SizeChangedFcn',@obj.resizeAx);
            obj.haxes = axes('Units','points','Parent',obj.axPanel);

            
        end
        function resizeUI(obj,varargin)
            obj.handles.params.Visible = 'off';
            width = obj.handles.params.Position(3);
            height= obj.handles.params.Position(4);
            
            hSize = width/5;
            hSpace = width/60;
            hBorder = width/15;
            
            vFac = 1;
            if ispc
                vFac = 0.81*vFac;
            end
            vSize = 22*vFac;
            vSpace = 5*vFac;
            vBorderTop = 45*vFac;
            vBorder = 15*vFac;
            
            obj.handles.padm.Position = [hBorder+hSize+hSpace height-vBorderTop hSize vSize];
            obj.handles.padp.Position = [hBorder+2*hSize+2*hSpace height-vBorderTop hSize vSize];
            obj.handles.step.Position = [hBorder+3*hSize+3*hSpace height-vBorderTop hSize vSize];
            
            obj.handles.x.Position = [hBorder height-vBorderTop-vSize hSize vSize];
            obj.handles.nxm.Position = [hBorder+hSize+hSpace height-vBorderTop-vSize hSize vSize];
            obj.handles.nxp.Position = [hBorder+2*hSize+2*hSpace height-vBorderTop-vSize hSize vSize];
            obj.handles.dx.Position = [hBorder+3*hSize+3*hSpace height-vBorderTop-vSize hSize vSize];
            
            obj.handles.z.Position = [hBorder height-vBorderTop-2*vSize-vSpace hSize vSize];
            obj.handles.nzm.Position = [hBorder+hSize+hSpace height-vBorderTop-2*vSize-vSpace hSize vSize];
            obj.handles.nzp.Position = [hBorder+2*hSize+2*hSpace height-vBorderTop-2*vSize-vSpace hSize vSize];
            obj.handles.dz.Position = [hBorder+3*hSize+3*hSpace height-vBorderTop-2*vSize-vSpace hSize vSize];
            
            obj.handles.bh_x0.Position = [hBorder height-vBorderTop-3*vSize-vSpace-2*vBorder 2*hSize vSize];
            obj.handles.bhList.Position = [hBorder+hSpace+2*hSize height-vBorderTop-3*vSize-vSpace-2*vBorder 2*hSize vSize];
            
            obj.handles.origin.Position = [hBorder height-vBorderTop-4*vSize-2*vSpace-2*vBorder hSize vSize];
            obj.handles.x0x.Position = [hBorder+hSize+hSpace height-vBorderTop-4*vSize-2*vSpace-2*vBorder hSize vSize];
            obj.handles.x0y.Position = [hBorder+2*hSize+2*hSpace height-vBorderTop-4*vSize-2*vSpace-2*vBorder hSize vSize];
            obj.handles.x0z.Position = [hBorder+3*hSize+3*hSpace height-vBorderTop-4*vSize-2*vSpace-2*vBorder hSize vSize];
            
            obj.handles.flip.Position = [hBorder 2*vBorder+vSize 1.75*hSize vSize];
            obj.handles.numCores.Position = [2*hBorder+1.75*hSize 2*vBorder+vSize 1.25*hSize vSize];
            obj.handles.ncoreList.Position = [2.5*hBorder+3*hSize 2*vBorder+vSize 0.9*hSize vSize];
            
            obj.handles.showFit.Position = [hSize vBorder 3*hSize vSize];
            
            obj.handles.params.Visible = 'on';
        end
        function resizeAx(obj,varargin)
            obj.axPanel.Visible = 'off';
            width = obj.axPanel.Position(3);
            height= obj.axPanel.Position(4);

            axBorder = obj.haxes.Position(2);
            
            obj.haxes.Position = [axBorder axBorder width-1.5*axBorder height-1.5*axBorder];
            obj.axPanel.Visible = 'on';
        end
        function inputEdit(obj,varargin)
            obj.dx = str2double( obj.handles.dx.String );
            obj.dz = str2double( obj.handles.dz.String );
            obj.grid.bord(1) = str2double( obj.handles.nxm.String );
            obj.grid.bord(2) = str2double( obj.handles.nxp.String );
            obj.grid.bord(3) = str2double( obj.handles.nzm.String );
            obj.grid.bord(4) = str2double( obj.handles.nzp.String );
            obj.grid.flip = obj.handles.flip.Value;
            val = obj.handles.bhList.Value;
            if val~=obj.grid.borehole_x0
                obj.grid.borehole_x0 = val;
                obj.data.order_boreholes = fliplr(obj.data.order_boreholes);
                obj.data.order_planes = fliplr(obj.data.order_planes);
                obj.data.planes = fliplr(obj.data.planes);
            end
            bhNo = obj.data.order_boreholes(1);
            obj.grid.x0(1) = obj.data.boreholes(bhNo).X;
            obj.grid.x0(2) = obj.data.boreholes(bhNo).Y;
            obj.grid.x0(3) = obj.data.boreholes(bhNo).Z;
            obj.handles.x0x.String = num2str(obj.grid.x0(1));
            obj.handles.x0y.String = num2str(obj.grid.x0(2));
            obj.handles.x0z.String = num2str(obj.grid.x0(3));
            obj.grid.nthreads = obj.handles.ncoreList.Value;
            
            obj.updateProj();
            
            obj.notify('gridEdited')
        end
        function updateProj(obj,varargin)

            z0 = obj.grid.x0(3);
            [az,dip] = obj.get_azimuth_dip();
            l = 0;
            for n=1:obj.data.n_planes
                no_f = obj.data.order_boreholes(n);
                no_plan = obj.data.order_planes(n);
                ind = obj.data.Tx_no_plan==no_plan;
                oo = [obj.data.boreholes(no_f).X obj.data.boreholes(no_f).Y 0];

                obj.grid.Tx(ind,:) = Grid.transl_rotat(obj.data.Tx(ind,:), oo,...
                    az(no_plan), dip(no_plan));
                obj.grid.TxCosDir(ind,:) = Grid.transl_rotat(obj.data.TxCosDir(ind,:),...
                    [0 0 0], az(no_plan), dip(no_plan));
                if n>1
                    l = l + obj.data.planes(n-1).l;
                    obj.grid.Tx(ind,1) = obj.grid.Tx(ind,1)+l;
                end
                
%                 figure
%                 subplot(121)
%                 plot(obj.data.Tx(ind,1),obj.data.Tx(ind,3),'ro')
%                 hold on
%                 
%                 subplot(122)
%                 plot(obj.grid.Tx(ind,1),obj.grid.Tx(ind,3),'ro')
%                 hold on
                
                ind = obj.data.Rx_no_plan==no_plan;
                obj.grid.Rx(ind,:) = Grid.transl_rotat(obj.data.Rx(ind,:), oo,...
                    az(no_plan), dip(no_plan));
                obj.grid.RxCosDir(ind,:) = Grid.transl_rotat(obj.data.RxCosDir(ind,:),...
                    [0 0 0], az(no_plan), dip(no_plan));
                if n>1
                    obj.grid.Rx(ind,1) = obj.grid.Rx(ind,1)+l;
                end
                
%                 subplot(121)
%                 plot(obj.data.Rx(ind,1),obj.data.Rx(ind,3),'x')
%                 hold off
% 
%                 subplot(122)
%                 plot(obj.grid.Rx(ind,1),obj.grid.Rx(ind,3),'x')
%                 hold off
            end
            obj.grid.Tx(:,3) = obj.grid.Tx(:,3);
            obj.grid.Rx(:,3) = obj.grid.Rx(:,3);
            
            % y = 0 in 2D
            obj.grid.Tx(:,2) = 0;
            obj.grid.Rx(:,2) = 0;
            
            obj.grid.in = obj.data.in;
            
            obj.xshift = min([obj.grid.Tx(:,1);obj.grid.Rx(:,1)]);
            obj.grid.Tx(:,1) = obj.grid.Tx(:,1)-obj.xshift;
            obj.grid.Rx(:,1) = obj.grid.Rx(:,1)-obj.xshift;
            xmin = 0;
            xmax = max([obj.grid.Tx(:,1);obj.grid.Rx(:,1)]);
            
            zmin = min([obj.grid.Tx(:,3);obj.grid.Rx(:,3)]);
            zmax = max([obj.grid.Tx(:,3);obj.grid.Rx(:,3)]);
            
            zmin = z0 - ceil( (z0-zmin)/obj.dz )*obj.dz;
            zmax = z0 + ceil( (zmax-z0)/obj.dz )*obj.dz;
            
            nx = ceil((xmax-xmin)/obj.dx);
            nz = ceil((zmax-zmin)/obj.dz);
            nxm = obj.grid.bord(1);
            nxp = obj.grid.bord(2);
            nzm = obj.grid.bord(3);
            nzp = obj.grid.bord(4);
            obj.grid.grx = xmin+obj.dx*((0-nxm):(nx+nxp));
            obj.grid.grz = zmin+obj.dz*((0-nzm):(nz+nzp));
            obj.grid.grx = obj.grid.grx';
            obj.grid.grz = obj.grid.grz';

            if obj.grid.flip==1
                obj.grid.grx = flipud(-obj.grid.grx);
                obj.grid.Tx(obj.grid.in,1) = -obj.grid.Tx(obj.grid.in,1);
                obj.grid.Rx(obj.grid.in,1) = -obj.grid.Rx(obj.grid.in,1);
                obj.grid.TxCosDir(obj.grid.in,1) = -obj.grid.TxCosDir(obj.grid.in,1);
                obj.grid.RxCosDir(obj.grid.in,1) = -obj.grid.RxCosDir(obj.grid.in,1);
            end
            
            % show grid
            
            xlim = [min(obj.grid.grx) max(obj.grid.grx)];
            zlim = [min(obj.grid.grz) max(obj.grid.grz)];
            xx1 = repmat(obj.grid.grx(:).',3,1) ;
            zz1 = repmat([zlim(:) ; nan],1,numel(obj.grid.grx)) ;
            xx2 = repmat([xlim(:) ; nan],1,numel(obj.grid.grz)) ;
            zz2 = repmat(obj.grid.grz(:).',3,1) ;
            xx1 = [xx1 xx2];
            zz1 = [zz1 zz2];
            
            plot(obj.haxes,xx1, zz1,'Color',[0.5 0.5 0.5])
            hold(obj.haxes,'on')
            plot(obj.haxes,obj.grid.Tx(obj.grid.in,1),obj.grid.Tx(obj.grid.in,3),'b*')
            plot(obj.haxes,obj.grid.Rx(obj.grid.in,1),obj.grid.Rx(obj.grid.in,3),'gx')
            hold(obj.haxes,'off')
            axis(obj.haxes,'equal')
            axis(obj.haxes,'tight')
            xlabel(obj.haxes,'X')
            ylabel(obj.haxes,'Z')
        end
        
        function [az,dip] = get_azimuth_dip(obj)
            az = zeros(1,obj.data.n_planes);
            dip = zeros(1,obj.data.n_planes);
            for nn=1:obj.data.n_planes
                no1 = obj.data.order_boreholes(nn);
                no2 = obj.data.order_boreholes(nn+1);
                
                % determination de l'azimuth du plan
                origin(1) = obj.data.boreholes(no1).X;
                origin(2) = obj.data.boreholes(no1).Y;
                origin(3) = obj.data.boreholes(no1).Z;
                
                % eq de la droite du borehole no2
                d=sqrt( sum([obj.data.boreholes(no2).Xmax-obj.data.boreholes(no2).X ...
                    obj.data.boreholes(no2).Ymax-obj.data.boreholes(no2).Y ...
                    obj.data.boreholes(no2).Zmax-obj.data.boreholes(no2).Z].^2) );
                % cosinus directeurs
                l=(obj.data.boreholes(no2).Xmax-obj.data.boreholes(no2).X)/d;
                m=(obj.data.boreholes(no2).Ymax-obj.data.boreholes(no2).Y)/d;
                n=(obj.data.boreholes(no2).Zmax-obj.data.boreholes(no2).Z)/d;
                % intersection  Z = origine(3)
                x2 = obj.data.boreholes(no2).X + l*(origin(3)-obj.data.boreholes(no2).Z)/n;
                y2 = obj.data.boreholes(no2).Y + m*(origin(3)-obj.data.boreholes(no2).Z)/n;
                %     az(nn) = atan2(y2-origine(2), x2-origine(1));
                
                no_plan = obj.data.order_planes(nn);
                d = sum( obj.data.planes(no_plan).x0 .* obj.data.planes(no_plan).a );
                x = d/obj.data.planes(no_plan).a(1);
                y = d/obj.data.planes(no_plan).a(2);
                az(no_plan) = atan2(y,x);
                
                rot = [cos(az(no_plan)) -sin(az(no_plan)); sin(az(no_plan)) cos(az(no_plan))];
                xr = [x2-origin(1) y2-origin(2)]*rot';
                if ( xr(1) < 0 )
                    az(no_plan) = az(no_plan) + pi;
                end
                
                % boreholes verticaux
                dip(nn) = 0;
            end
        end
        function showFit(obj,varargin)
            [Tx_p, Tx_no_plan] = Grid.proj_planes(obj.data.Tx, obj.data.planes);
            [Rx_p, Rx_no_plan] = Grid.proj_planes(obj.data.Rx, obj.data.planes);
            for n=1:obj.data.n_planes
                figure
                
                no_plan = obj.data.order_planes(n);
                ind = Tx_no_plan==no_plan;
                
                dTx = sqrt( sum( (obj.data.Tx(ind,:)-Tx_p(ind,:)).^2, 2) );
            
                subplot(221)
                plot(dTx,'o')
                title('Distance between original and projected Tx')
                
                % Cosinus directeurs des antennes après rotation
                
                subplot(223)
                plot(obj.grid.TxCosDir(ind,1),'b')
                hold on
                plot(obj.grid.TxCosDir(ind,2),'r')
                plot(obj.grid.TxCosDir(ind,3),'g')
                hold off
                title('Tx direction cosines after rotation')
                
                ind = Rx_no_plan==no_plan;
                dRx = sqrt( sum( (obj.data.Rx(ind,:)-Rx_p(ind,:)).^2, 2) );
                
                subplot(222)
                plot(dRx,'o')
                title('Distance between original and projected Rx')
                
                subplot(224)
                plot(obj.grid.RxCosDir(ind,1),'b')
                hold on
                plot(obj.grid.RxCosDir(ind,2),'r')
                plot(obj.grid.RxCosDir(ind,3),'g')
                hold off
                title('Rx direction cosines after rotation')
                
            end
        end
    end
end

