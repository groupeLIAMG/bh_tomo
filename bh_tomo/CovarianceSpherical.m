classdef CovarianceSpherical < Covariance
    %UNTITLED2 Summary of this class goes here
    %   Detailed explanation goes here
    
    properties (SetAccess=private)
        type
    end
    
    methods
        function obj = CovarianceSpherical(r,a,s)
            obj@Covariance(r,a,s);
            obj.type = CovarianceModels.spherical;
        end
        function k = compute(obj, x, x0)
            h = obj.compute_h(x, x0);
            k = obj.sill * (1-(1.5*min(h,1)/1-.5*(min(h,1)/1).^3));
        end
    end
    
end
