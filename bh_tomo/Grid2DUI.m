classdef Grid2DUI < handle
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
    end
    events
        gridEdited
    end

    methods
        function obj = Grid2DUI(f,type,fs,g,data)

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
    end

    methods (Static)
        function [grid,data] = buildGrid(data)
            grid = Grid2D();
            grid.bord = [1 1 1 1];
            grid.flip = 0;
            grid.borehole_x0 = 1;
            grid.x0 = [data.boreholes(1).X data.boreholes(1).Y data.boreholes(1).Z];
            grid.type = '2D';
            uTx = unique(data.Tx(data.in,:), 'rows');
            uRx = unique(data.Rx(data.in,:), 'rows');

            [data.x0,data.a]=Grid.lsplane([uTx; uRx]);
            % data.x0 : Centroid of the data = point on the best-fit plane
            % data.a  : Direction cosines of the normal to the best-fit plane
            if data.a(3)<0, data.a=-data.a; end

            data.Tx_p = Grid.proj_plane(data.Tx, data.x0, data.a);
            data.Rx_p = Grid.proj_plane(data.Rx, data.x0, data.a);
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

            bhNames = cell(length(obj.data.boreholes),1);
            for n=1:length(obj.data.boreholes)
                bhNames{n} = obj.data.boreholes(n).name;
            end
            obj.handles.bhList = uicontrol('Style','popupmenu',...
                'String',bhNames,...
                'Value',obj.grid.borehole_x0,...
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

            obj.handles.flip.Position = [1.5*hSize 2*vBorder+vSize 2*hSize vSize];
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
            bhNo = obj.handles.bhList.Value;
            obj.grid.x0(1) = obj.data.boreholes(bhNo).X;
            obj.grid.x0(2) = obj.data.boreholes(bhNo).Y;
            obj.grid.x0(3) = obj.data.boreholes(bhNo).Z;
            obj.handles.x0x.String = num2str(obj.grid.x0(1));
            obj.handles.x0y.String = num2str(obj.grid.x0(2));
            obj.handles.x0z.String = num2str(obj.grid.x0(3));

            obj.updateProj();

            obj.notify('gridEdited')
        end
        function updateProj(obj,varargin)
            [az,dip] = obj.get_azimuth_dip();
            obj.grid.Tx = Grid.transl_rotat(obj.data.Tx_p, obj.grid.x0, az, dip);
            obj.grid.Rx = Grid.transl_rotat(obj.data.Rx_p, obj.grid.x0, az, dip);
            obj.grid.TxCosDir = Grid.transl_rotat(obj.data.TxCosDir, [0 0 0], az, dip);
            obj.grid.RxCosDir = Grid.transl_rotat(obj.data.RxCosDir, [0 0 0], az, dip);
            if ~isnan(obj.data.Tx_Z_water(1,1))
                obj.grid.Tx_Z_water = Grid.transl_rotat(obj.data.Tx_Z_water, origine, az, dip);
            end
            if ~isnan(obj.data.Rx_Z_water(1,1))
                obj.grid.Rx_Z_water = Grid.transl_rotat(obj.data.Rx_Z_water, origine, az, dip);
            end
            obj.grid.in = obj.data.in;

            xmin = min([obj.grid.Tx(obj.grid.in,1);obj.grid.Rx(obj.grid.in,1)]) - 0.5*obj.dx;
            xmax = max([obj.grid.Tx(obj.grid.in,1);obj.grid.Rx(obj.grid.in,1)]) + 0.5*obj.dx;
            nx = ceil((xmax-xmin)/obj.dx);

            zmin = min([obj.grid.Tx(obj.grid.in,3);obj.grid.Rx(obj.grid.in,3)]) - 0.5*obj.dz;
            zmax = max([obj.grid.Tx(obj.grid.in,3);obj.grid.Rx(obj.grid.in,3)]) + 0.5*obj.dz;
            nz = ceil((zmax-zmin)/obj.dz);

            nxm = obj.grid.bord(1);
            nxp = obj.grid.bord(2);
            nzm = obj.grid.bord(3);
            nzp = obj.grid.bord(4);

            obj.grid.grx = xmin+obj.dx*((0-nxm):(nx+nxp));
            obj.grid.grz = zmin+obj.dz*((0-nzm):(nz+nzp));
            obj.grid.grx = obj.grid.grx';
            obj.grid.grz = obj.grid.grz';

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
        function [az,dip] = get_azimuth_dip(obj,varargin)
            d = sum( obj.data.x0 .* obj.data.a );
            x = d/obj.data.a(1);
            y = d/obj.data.a(2);
            az = atan2(y,x);
            dip = asin(obj.data.a(3));
            flip = obj.handles.flip.Value;
            az = az + flip*pi;
        end
        function showFit(obj,varargin)
            figure

            % distance entre les pts originaux et les pts projetes dans le plan

            dTx = sqrt( sum( (obj.data.Tx-obj.data.Tx_p).^2, 2) );

            subplot(321)
            plot(dTx,'o')
            title('Distance between original and projected Tx')

            dRx = sqrt( sum( (obj.data.Rx-obj.data.Rx_p).^2, 2) );

            subplot(322)
            plot(dRx,'o')
            title('Distance between original and projected Rx')

            l_orig = sqrt( sum( (obj.data.Tx-obj.data.Rx).^2,2 ) );
            l_new = sqrt( sum( (obj.data.Tx_p-obj.data.Rx_p).^2,2 ) );

            % Cosinus directeurs des antennes apr�s rotation

            subplot(323)
            plot(obj.grid.TxCosDir(:,1),'b')
            hold on
            plot(obj.grid.TxCosDir(:,2),'r')
            plot(obj.grid.TxCosDir(:,3),'g')
            hold off
            title('Tx direction cosines after rotation')

            subplot(324)
            plot(obj.grid.RxCosDir(:,1),'b')
            hold on
            plot(obj.grid.RxCosDir(:,2),'r')
            plot(obj.grid.RxCosDir(:,3),'g')
            hold off
            title('Rx direction cosines after rotation')

            % Erreur relative sur la longueur des rais après projection
            diff = 100*(l_orig - l_new)./l_orig;
            subplot(325)
            plot(diff)
            title('Relative error on ray length after projection [%]')
        end
    end
end
