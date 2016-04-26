classdef CovarianceCubic < Covariance & matlab.mixin.Copyable
    %UNTITLED2 Summary of this class goes here
    %   Detailed explanation goes here
    
    % B. Giroux
    % INRS-ETE
    % 2016-04-26
    methods
        function obj = CovarianceCubic(r,a,s)
            obj@Covariance(r,a,s);
            obj.type = CovarianceModels.Cubic;
        end
        function k = compute(obj, x, x0)
            h = obj.compute_h(x, x0);
            k = obj.sill * (1-3*min(h,1).^2+2*min(h,1).^3);
        end
    end
    
end
