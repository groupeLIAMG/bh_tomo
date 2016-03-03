classdef Borehole < handle
    %BOREHOLE Class to hold borehole parameters
    
    properties
        name
        X
        Y
        Z
        Xmax
        Ymax
        Zmax
        Z_surf
        Z_water
        scont
        acont
        diam
        fdata
    end
    
    methods
        function obj = Borehole(n)
            obj.name = n;
            obj.X = 0;
            obj.Y = 0;
            obj.Z = 0;
            obj.Xmax = 0;
            obj.Ymax = 0;
            obj.Zmax = 0;
            obj.Z_surf = 0;
            obj.Z_water = NaN;
            obj.diam = 0;
            obj.scont = [];
            obj.acont = [];
            obj.fdata = [0 0 0;0 0 0];
        end
        function set.name(obj, n)
            if ischar(n)
                obj.name = n;
            else
                error('Borehole name must be a string')
            end
        end
        function set.X(obj,value)
            if isnumeric(value)
                obj.X = value;
            else
                error('Value must be numeric')
            end
        end
        function set.Y(obj,value)
            if isnumeric(value)
                obj.Y = value;
            else
                error('Value must be numeric')
            end
        end
        function set.Z(obj,value)
            if isnumeric(value)
                obj.Z = value;
            else
                error('Value must be numeric')
            end
        end
        function set.Xmax(obj,value)
            if isnumeric(value)
                obj.Xmax = value;
            else
                error('Value must be numeric')
            end
        end
        function set.Ymax(obj,value)
            if isnumeric(value)
                obj.Ymax = value;
            else
                error('Value must be numeric')
            end
        end
        function set.Zmax(obj,value)
            if isnumeric(value)
                obj.Zmax = value;
            else
                error('Value must be numeric')
            end
        end
        function set.Z_surf(obj,value)
            if isnumeric(value)
                obj.Z_surf = value;
            else
                error('Value must be numeric')
            end
        end
        function set.Z_water(obj,value)
            if isnumeric(value)
                obj.Z_water = value;
            else
                error('Value must be numeric')
            end
        end
        function set.diam(obj,value)
            if isnumeric(value)
                obj.diam = value;
            else
                error('Value must be numeric')
            end
        end
        function set.fdata(obj,value)
            if isnumeric(value)
                if size(value,2)~=3
                    error('Borehole data should be nData x 3')
                end
                obj.fdata = value;
            else
                error('Value must be numeric')
            end
        end
    end
    
end
