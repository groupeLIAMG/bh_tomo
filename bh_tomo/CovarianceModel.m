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
        use_xi
        use_tilt
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
            obj.use_xi = 0;
            obj.use_tilt = 0;
        end
        
        function Cm = computeCm(obj,x,x0)
            Cm = obj.covar(1).compute(x,x0);
            for n=2:numel(obj.covar)
                Cm = Cm + obj.covar(n).compute(x,x0);
            end
            if obj.nugget_m~=0
                Cm = Cm + obj.nugget_m*speye(size(Cm,1));
            end
            if obj.use_xi==1
                Cx = obj.covar_xi(1).compute(x,x0);
                for n=2:numel(obj.covar_xi)
                    Cx = Cx + obj.covar_xi(n).compute(x,x0);
                end
                if obj.nugget_xi~=0
                    Cx = Cx + obj.nugget_xi*speye(size(Cx,1));
                end
                if obj.use_tilt==1
                    Ct = obj.covar_tilt(1).compute(x,x0);
                    for n=2:numel(obj.covar_tilt)
                        Ct = Ct + obj.covar_tilt(n).compute(x,x0);
                    end
                    if obj.nugget_tilt~=0
                        Ct = Ct + obj.nugget_tilt*speye(size(Ct,1));
                    end
                    Cm = sparse([Cm zeros(size(Cx)) zeros(size(Ct));
                        zeros(size(Cm)) Cx zeros(size(Ct));
                        zeros(size(Cm)) zeros(size(Cx)) Ct]);
                else
                    Cm = sparse([Cm zeros(size(Cx));zeros(size(Cm)) Cx]);
                end
            end
        end
    end
    
end

