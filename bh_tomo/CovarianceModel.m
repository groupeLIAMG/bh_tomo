classdef CovarianceModel < matlab.mixin.Copyable
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
        function obj = CovarianceModel(varargin)
            obj.covar = Covariance.empty;
            obj.covar_xi = Covariance.empty;
            obj.covar_tilt = Covariance.empty;
            if nargin == 0
                obj.nugget_d = 0;
                obj.nugget_m = 0;
                obj.nugget_xi = 0;
                obj.nugget_tilt = 0;
                obj.use_c0 = 0;
                obj.use_xi = 0;
                obj.use_tilt = 0;
            else
                if ~isstruct(varargin{1})
                    error('Invalid input')
                end
                s = varargin{1};
                % old indices in covardm are
                % 1 -> nugget
                % 2 -> exponential
                % 3 -> gaussian
                % 4 -> spherical
                % 5 -> linear
                % 6 -> cubic
                % 7 -> thin plate
                % 8 -> gravimetric
                % 9 -> magnetic
                % 10 -> hole effect sine
                % 11 -> hole effect cosine
                %
                % new indices are (CovarianceModels enum)
                % Cubic (1)
                % Spherical (2)
                % Gaussian (3)
                % Exponential (4)
                % Linear (5)
                % Thin_Plate (6)
                % Gravimetric (7)
                % Magnetic (8)
                % Hole_Effect_Sine (9)
                % Hole_Effect_Cosine (10)
                % Nugget (11)
                conv = int8([11 4 3 2 5 1 6 7 8 9 10]);
                if isfield(s, 'model')
                    for n=1:size(s.model,1)
                        % we have 2d models
                        obj.covar(n) = CovarianceModels.buildCov(conv(s.model(n,1)),...
                            s.model(n,3:-1:2), s.model(n,4), s.c(n));
                    end
                elseif isfield(s, 'modele')
                    for n=1:size(s.modele,1)
                        % we have 2d models
                        obj.covar(n) = CovarianceModels.buildCov(conv(s.modele(n,1)),...
                            s.modele(n,3:-1:2), s.modele(n,4), s.c(n));
                    end
                else
                    error('Invalid input')
                end
                if isfield(s, 'model_xi')
                    for n=1:size(s.model_xi,1)
                        obj.covar_xi(n) = CovarianceModels.buildCov(conv(s.model_xi(n,1)),...
                            s.model_xi(n,3:-1:2), s.model_xi(n,4), s.c_xi(n));
                    end
                elseif isfield(s, 'modele_xi')
                    for n=1:size(s.modele_xi,1)
                        obj.covar_xi(n) = CovarianceModels.buildCov(conv(s.modele_xi(n,1)),...
                            s.modele_xi(n,3:-1:2), s.modele_xi(n,4), s.c_xi(n));
                    end
                else
                    obj.covar_xi = Covariance.empty;
                end
                if isfield(s, 'model_th')
                    for n=1:size(s.model_th,1)
                        obj.covar_tilt(n) = CovarianceModels.buildCov(conv(s.model_th(n,1)),...
                            s.model_th(n,3:-1:2), s.model_th(n,4), s.c_th(n));
                    end
                elseif isfield(s, 'modele_th')
                    for n=1:size(s.modele_th,1)
                        obj.covar_tilt(n) = CovarianceModels.buildCov(conv(s.modele_th(n,1)),...
                            s.modele_th(n,3:-1:2), s.modele_th(n,4), s.c_th(n));
                    end
                else
                    obj.covar_tilt = Covariance.empty;
                end
                if isfield(s, 'nugget_t')
                    obj.nugget_d = s.nugget_t;
                elseif isfield(s, 'pepite_t')
                    obj.nugget_d = s.pepite_t;
                else
                    error('Invalid input')
                end
                if isfield(s, 'nugget_l')
                    obj.nugget_m = s.nugget_l;
                elseif isfield(s, 'pepite_l')
                    obj.nugget_m = s.pepite_l;
                else
                    error('Invalid input')
                end
                if isfield(s, 'nugget_xi')
                    obj.nugget_xi = s.nugget_xi;
                else
                    obj.nugget_xi = 0;
                end
                if isfield(s, 'nugget_th')
                    obj.nugget_tilt = s.nugget_th;
                else
                    obj.nugget_tilt = 0;
                end
                obj.use_c0 = s.use_c0;
                if isfield(s, 'aniso')
                    obj.use_xi = s.aniso;
                else
                    obj.use_xi = 0;
                end
                obj.use_tilt = 0;
                
            end
        end
        
        function Cm = compute(obj,x,x0)
            Cm = obj.covar(1).compute(x,x0);
            for n=2:numel(obj.covar)
                Cm = Cm + obj.covar(n).compute(x,x0);
            end
            if obj.nugget_m~=0
                Cm = Cm + obj.nugget_m*eye(size(Cm,1));
            end
            if obj.use_xi==1
                Cx = obj.covar_xi(1).compute(x,x0);
                for n=2:numel(obj.covar_xi)
                    Cx = Cx + obj.covar_xi(n).compute(x,x0);
                end
                if obj.nugget_xi~=0
                    Cx = Cx + obj.nugget_xi*eye(size(Cx,1));
                end
                if obj.use_tilt==1
                    Ct = obj.covar_tilt(1).compute(x,x0);
                    for n=2:numel(obj.covar_tilt)
                        Ct = Ct + obj.covar_tilt(n).compute(x,x0);
                    end
                    if obj.nugget_tilt~=0
                        Ct = Ct + obj.nugget_tilt*eye(size(Ct,1));
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

