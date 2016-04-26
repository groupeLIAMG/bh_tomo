classdef CovarianceExponential < Covariance & matlab.mixin.Copyable
    %UNTITLED2 Summary of this class goes here
    %   Detailed explanation goes here
    
    % B. Giroux
    % INRS-ETE
    % 2016-04-26
    methods
        function obj = CovarianceExponential(r,a,s)
            obj@Covariance(r,a,s);
            obj.type = CovarianceModels.Exponential;
        end
        function k = compute(obj, x, x0)
            h = obj.compute_h(x, x0);
            k = obj.sill * exp(-h);
        end
    end
    
end
