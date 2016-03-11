classdef Grid < handle
    %GRID Base class for 2D and 3D grids
    
    properties
        grx
        gry
        grz
        cont
        Tx
        Rx
        TxCosDir
        RxCosDir
        x0
        bord
        Tx_Z_water
        Rx_Z_water
        in
    end
    methods
        function set.grx(obj,g)
            if ~isnumeric(g)
                error('Grid coordinates must be numeric')
            end
            if any(abs(diff(g)-(g(2)-g(1)))>10*eps)
                error('Grid step size should be constant')
            end
            obj.grx = g;
        end
        function set.gry(obj,g)
            if ~isnumeric(g)
                error('Grid coordinates must be numeric')
            end
            if any(abs(diff(g)-(g(2)-g(1)))>10*eps)
                error('Grid step size should be constant')
            end
            obj.gry = g;
        end
        function set.grz(obj,g)
            if ~isnumeric(g)
                error('Grid coordinates must be numeric')
            end
            if any(abs(diff(g)-(g(2)-g(1)))>10*eps)
                error('Grid step size should be constant')
            end
            obj.grz = g;
        end
    end
    methods (Abstract)
        varargout = raytrace(obj, varargin)
    end
    
end

