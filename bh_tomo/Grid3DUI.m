classdef Grid3DUI < handle
    %UNTITLED Summary of this class goes here
    %   Detailed explanation goes here

    properties (Dependent)
        Position
    end
    properties (SetAccess=private)
        type
        grid
        dx
        dy
        dz
        axPanel
    end
    properties (Access=private)
        handles
        data
        haxes
        haxView
    end
    events
        gridEdited
    end

    methods
        function obj = Grid3DUI(f,type,fs,g,data)

            obj.type = type;
            if isempty(g)
                obj.grid = Grid3D();
            else
                if isa(g,'Grid3D')
                    obj.grid = g;
                else
                    error('grid not an instance of Grid3D')
                end
            end
            obj.data = data;

            if isempty(g.grx) % grid not initialized yet
                obj.dx = 1;
                obj.dy = 1;
                obj.dz = 1;
            else
                obj.dx = g.grx(2)-g.grx(1);
                obj.dy = g.gry(2)-g.gry(1);
                obj.dz = g.grz(2)-g.grz(1);
            end

            obj.addComponents(f,fs)
            obj.handles.params.Visible = 'on';
            obj.axPanel.Visible = 'on';

            obj.updateGrid()

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
        function set.dy(obj,dy)
            if isnumeric(dy)
                obj.dy = dy;
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
        function [grid] = buildGrid(data)
            grid = Grid3D();
            grid.bord = [1 1 1 1 1 1];
            grid.Tx = data.Tx;
            grid.Rx = data.Rx;
            grid.in = data.in;
            grid.TxCosDir = data.TxCosDir;
            grid.RxCosDir = data.RxCosDir;
            grid.type = '3D';
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
            obj.handles.y = uicontrol('Style','text',...
                'String','Y',...
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
            obj.handles.nym = uicontrol('Style','edit',...
                'String',num2str(obj.grid.bord(3)),...
                'FontSize',fs,...
                'Units','points',...
                'Callback',@obj.inputEdit,...
                'Parent',obj.handles.params);
            obj.handles.nyp = uicontrol('Style','edit',...
                'String',num2str(obj.grid.bord(4)),...
                'FontSize',fs,...
                'Units','points',...
                'Callback',@obj.inputEdit,...
                'Parent',obj.handles.params);
            obj.handles.nzm = uicontrol('Style','edit',...
                'String',num2str(obj.grid.bord(5)),...
                'FontSize',fs,...
                'Units','points',...
                'Callback',@obj.inputEdit,...
                'Parent',obj.handles.params);
            obj.handles.nzp = uicontrol('Style','edit',...
                'String',num2str(obj.grid.bord(6)),...
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
            obj.handles.dy = uicontrol('Style','edit',...
                'String',num2str(obj.dy),...
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


            obj.axPanel = uipanel(f,'Title','Grid View',...
                'FontSize',fs+1,...
                'Units','points',...
                'Visible','off',...
                'SizeChangedFcn',@obj.resizeAx);
            obj.haxes = axes('Units','points','Parent',obj.axPanel);
            obj.haxView = uicontrol('Style','popupmenu',...
                'String',{'XY Plane','XZ Plane','YZ Plane'},...
                'Value',2,...
                'Units','points',...
                'FontSize',fs,...
                'Callback',@obj.axView,...
                'Parent',obj.axPanel);


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

            obj.handles.y.Position = [hBorder height-vBorderTop-2*vSize-vSpace hSize vSize];
            obj.handles.nym.Position = [hBorder+hSize+hSpace height-vBorderTop-2*vSize-vSpace hSize vSize];
            obj.handles.nyp.Position = [hBorder+2*hSize+2*hSpace height-vBorderTop-2*vSize-vSpace hSize vSize];
            obj.handles.dy.Position = [hBorder+3*hSize+3*hSpace height-vBorderTop-2*vSize-vSpace hSize vSize];

            obj.handles.z.Position = [hBorder height-vBorderTop-3*vSize-2*vSpace hSize vSize];
            obj.handles.nzm.Position = [hBorder+hSize+hSpace height-vBorderTop-3*vSize-2*vSpace hSize vSize];
            obj.handles.nzp.Position = [hBorder+2*hSize+2*hSpace height-vBorderTop-3*vSize-2*vSpace hSize vSize];
            obj.handles.dz.Position = [hBorder+3*hSize+3*hSpace height-vBorderTop-3*vSize-2*vSpace hSize vSize];

            obj.handles.params.Visible = 'on';
        end
        function resizeAx(obj,varargin)
            obj.axPanel.Visible = 'off';
            width = obj.axPanel.Position(3);
            height= obj.axPanel.Position(4);

            axBorder = obj.haxes.Position(2);

            obj.haxes.Position = [axBorder axBorder width-1.5*axBorder height-2*axBorder];
            hSize = 2*obj.haxView.Extent(3);
            vSize = 1.5*obj.haxView.Extent(4);
            obj.haxView.Position = [width/2-hSize/2 height-3*axBorder/4 hSize vSize];
            obj.axPanel.Visible = 'on';
        end
        function inputEdit(obj,varargin)
            obj.dx = str2double( obj.handles.dx.String );
            obj.dy = str2double( obj.handles.dy.String );
            obj.dz = str2double( obj.handles.dz.String );
            obj.grid.bord(1) = str2double( obj.handles.nxm.String );
            obj.grid.bord(2) = str2double( obj.handles.nxp.String );
            obj.grid.bord(3) = str2double( obj.handles.nym.String );
            obj.grid.bord(4) = str2double( obj.handles.nyp.String );
            obj.grid.bord(5) = str2double( obj.handles.nzm.String );
            obj.grid.bord(6) = str2double( obj.handles.nzp.String );

            obj.updateGrid()

            obj.notify('gridEdited')
        end
        function updateGrid(obj)
            xmin = min([obj.grid.Tx(obj.grid.in,1);obj.grid.Rx(obj.grid.in,1)]) - 0.5*obj.dx;
            xmax = max([obj.grid.Tx(obj.grid.in,1);obj.grid.Rx(obj.grid.in,1)]) + 0.5*obj.dx;
            nx = ceil((xmax-xmin)/obj.dx);

            ymin = min([obj.grid.Tx(obj.grid.in,2);obj.grid.Rx(obj.grid.in,2)]) - 0.5*obj.dx;
            ymax = max([obj.grid.Tx(obj.grid.in,2);obj.grid.Rx(obj.grid.in,2)]) + 0.5*obj.dx;
            ny = ceil((ymax-ymin)/obj.dy);

            zmin = min([obj.grid.Tx(obj.grid.in,3);obj.grid.Rx(obj.grid.in,3)]) - 0.5*obj.dz;
            zmax = max([obj.grid.Tx(obj.grid.in,3);obj.grid.Rx(obj.grid.in,3)]) + 0.5*obj.dz;
            nz = ceil((zmax-zmin)/obj.dz);

            nxm = obj.grid.bord(1);
            nxp = obj.grid.bord(2);
            nym = obj.grid.bord(3);
            nyp = obj.grid.bord(4);
            nzm = obj.grid.bord(5);
            nzp = obj.grid.bord(6);

            obj.grid.grx = xmin+obj.dx*((0-nxm):(nx+nxp));
            obj.grid.gry = ymin+obj.dy*((0-nym):(ny+nyp));
            obj.grid.grz = zmin+obj.dz*((0-nzm):(nz+nzp));
            obj.grid.grx = obj.grid.grx';
            obj.grid.gry = obj.grid.gry';
            obj.grid.grz = obj.grid.grz';

            obj.axView()
        end
        function axView(obj,varargin)
            switch obj.haxView.Value
                case 1 % XY
                    xlim = [min(obj.grid.grx) max(obj.grid.grx)];
                    zlim = [min(obj.grid.gry) max(obj.grid.gry)];
                    xx1 = repmat(obj.grid.grx(:).',3,1) ;
                    zz1 = repmat([zlim(:) ; nan],1,numel(obj.grid.grx)) ;
                    xx2 = repmat([xlim(:) ; nan],1,numel(obj.grid.gry)) ;
                    zz2 = repmat(obj.grid.gry(:).',3,1) ;
                    xx1 = [xx1 xx2];
                    zz1 = [zz1 zz2];
                    xl = 'X';
                    yl = 'Y';
                case 2 % XZ
                    xlim = [min(obj.grid.grx) max(obj.grid.grx)];
                    zlim = [min(obj.grid.grz) max(obj.grid.grz)];
                    xx1 = repmat(obj.grid.grx(:).',3,1) ;
                    zz1 = repmat([zlim(:) ; nan],1,numel(obj.grid.grx)) ;
                    xx2 = repmat([xlim(:) ; nan],1,numel(obj.grid.grz)) ;
                    zz2 = repmat(obj.grid.grz(:).',3,1) ;
                    xx1 = [xx1 xx2];
                    zz1 = [zz1 zz2];
                    xl = 'X';
                    yl = 'Z';
                case 3 % YZ
                    xlim = [min(obj.grid.gry) max(obj.grid.gry)];
                    zlim = [min(obj.grid.grz) max(obj.grid.grz)];
                    xx1 = repmat(obj.grid.gry(:).',3,1) ;
                    zz1 = repmat([zlim(:) ; nan],1,numel(obj.grid.gry)) ;
                    xx2 = repmat([xlim(:) ; nan],1,numel(obj.grid.grz)) ;
                    zz2 = repmat(obj.grid.grz(:).',3,1) ;
                    xx1 = [xx1 xx2];
                    zz1 = [zz1 zz2];
                    xl = 'Y';
                    yl = 'Z';
            end

            plot(obj.haxes,xx1, zz1,'Color',[0.5 0.5 0.5])
            xlabel(obj.haxes,xl)
            ylabel(obj.haxes,yl)
            axis(obj.haxes,'equal')
            axis(obj.haxes,'tight')
        end
    end
end
