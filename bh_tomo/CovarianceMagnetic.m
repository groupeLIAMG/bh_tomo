classdef CovarianceMagnetic < Covariance & matlab.mixin.Copyable
    %UNTITLED2 Summary of this class goes here
    %   Cauchy with b=3/2
    
    methods
        function obj = CovarianceMagnetic(r,a,s)
            obj@Covariance(r,a,s);
            obj.type = CovarianceModels.Magnetic;
        end
        function k = compute(obj, x, x0)
            h = obj.compute_h(x, x0);
            k = obj.sill * ((h.^2+1).^(-1.5));
        end
    end
    
end
