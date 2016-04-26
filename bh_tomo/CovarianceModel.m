classdef CovarianceModel
    %UNTITLED Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        covar
        covar_xi    % anisotropy ratio
        covar_tilt  % anisotropy tilt angle
        nugget_d    % nugget -> data
        nugget_m    % nugget -> model
        nugget_xi
        nugget_tilt
        use_c0
    end
    
    methods
        function obj = CovarianceModel()
            obj.covar = Covariance.empty;
            obj.covar_xi = Covariance.empty;
            obj.covar_tilt = Covariance.empty;
            obj.nugget_d = 0;
            obj.nugget_m = 0;
            obj.nugget_xi = 0;
            obj.nugget_tilt = 0;
            obj.use_c0 = 0;
        end
    end
    
end

