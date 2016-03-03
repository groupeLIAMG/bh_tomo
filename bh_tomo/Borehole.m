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
        scont
        acont
        diam
        Z_water
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
            obj.scont = [];
            obj.acont = [];
            obj.diam = 0;
            obj.Z_water = NaN;
            obj.fdata = [];
        end
        function set.name(obj, n)
            if ischar(n)
                obj.name = n;
            else
                error('Borehole name must be a string')
            end
        end
    end
    
end
