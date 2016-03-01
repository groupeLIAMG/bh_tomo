classdef CovarianceGravimetric < Covariance
    %UNTITLED2 Summary of this class goes here
    %   Cauchy with b=1/2
    
    properties (SetAccess=private)
        type
    end
    
    methods
        function obj = CovarianceGravimetric(r,a,s)
            obj@Covariance(r,a,s);
            obj.type = CovarianceModels.gravimetric;
        end
        function k = compute(obj, x, x0)
            h = obj.compute_h(x, x0);
            k = obj.sill * ((h.^2+1).^(-0.5));
        end
    end
    
end
