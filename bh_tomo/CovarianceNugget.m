classdef CovarianceNugget < Covariance
    %UNTITLED2 Summary of this class goes here
    %   Detailed explanation goes here
    
    properties (SetAccess=private)
        type
    end
    
    methods
        function obj = CovarianceNugget(s)
            obj@Covariance([],[],s);
            obj.type = CovarianceModels.nugget;
        end
        function k = compute(obj, x, x0)
            d=size(x,2);
            obj.range = ones(1,d);
            if d==3
                obj.angle = zeros(1,d);
            else
                obj.angle = 0;
            end
            h = obj.compute_h(x, x0);
            k = obj.sill * (h==0);
        end
    end
    
end
