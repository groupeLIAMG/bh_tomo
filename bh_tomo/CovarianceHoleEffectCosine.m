classdef CovarianceHoleEffectCosine < Covariance
    %UNTITLED2 Summary of this class goes here
        
    properties (SetAccess=private)
        type
    end
    
    methods
        function obj = CovarianceHoleEffectCosine(r,a,s)
            obj@Covariance(r,a,s);
            obj.type = CovarianceModels.hole_effect_cosine;
        end
        function k = compute(obj, x, x0)
            h = obj.compute_h(x, x0);
            k = obj.sill * cos(h*2*pi);
        end
    end
    
end
