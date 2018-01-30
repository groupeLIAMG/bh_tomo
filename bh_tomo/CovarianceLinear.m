classdef CovarianceLinear < Covariance & matlab.mixin.Copyable
    %UNTITLED2 Summary of this class goes here
    %   Detailed explanation goes here
    
    % B. Giroux
    % INRS-ETE
    % 2016-04-26
    methods
        function obj = CovarianceLinear(r,a,s)
            obj@Covariance(r,a,s);
            obj.type = CovarianceModels.Linear;
        end
        function k = compute(obj, x, x0)
            h = obj.compute_h(x, x0);
            k = obj.sill * (1-h);
        end
    end
    
end
