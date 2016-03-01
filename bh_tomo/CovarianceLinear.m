classdef CovarianceLinear < Covariance
    %UNTITLED2 Summary of this class goes here
    %   Detailed explanation goes here
    
    properties (SetAccess=private)
        type
    end
    
    methods
        function obj = CovarianceLinear(r,a,s)
            obj@Covariance(r,a,s);
            obj.type = CovarianceModels.linear;
        end
        function k = compute(obj, x, x0)
            h = obj.compute_h(x, x0);
            k = obj.sill * (1-h);
        end
    end
    
end
