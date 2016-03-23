classdef Model < handle
    %UNTITLED Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        name
        mogs         %
        boreholes    % list of boreholes
        type
        grid
        tt_covar     % covariance model - travel time
        amp_covar    % covariance model - amplitudes
        inv_res      % inversion results
    end
    
    methods
        function obj = Model(n)
            obj.name = n;
            obj.grid = Grid.empty;
        end
        function set.name(obj, n)
            if ischar(n)
                obj.name = n;
            else
                error('Model name must be a string')
            end
        end
    end
    
end

