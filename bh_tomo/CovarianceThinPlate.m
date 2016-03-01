classdef CovarianceThinPlate < Covariance
    %UNTITLED2 Summary of this class goes here
    %   Detailed explanation goes here
    
    properties (SetAccess=private)
        type
    end
    
    methods
        function obj = CovarianceThinPlate(r,a,s)
            obj@Covariance(r,a,s);
            obj.type = CovarianceModels.thin_plate;
        end
        function k = compute(obj, x, x0)
            h = obj.compute_h(x, x0);
            k = obj.sill * ((h.^2).*log(max(h,eps)));
        end
    end
    
end
