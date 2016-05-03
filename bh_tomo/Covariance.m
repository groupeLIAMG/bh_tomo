classdef Covariance < matlab.mixin.Heterogeneous & matlab.mixin.Copyable
    %COVARIANCE Base class for covariance models
    
    % Code in part from covardm from Denis Marcotte
    %
    
    % B. Giroux
    % INRS-ETE
    % 2016-04-26
    properties
        range
        angle   % angle in degrees
        sill
        type
    end
    methods
        function obj = Covariance(r,a,s)
            if nargin > 0
                
                if numel(r)==2 && numel(a)~=1
                    error('angle must be a scalar for 2D media')
                elseif numel(r)==3 && numel(a)~=3
                    error('3 angle values needed for 3D media')
                end
                
                obj.range = r;
                obj.angle = a;
                obj.sill = s;
            end
            obj.type = CovarianceModels.empty;
        end
        function set.range(obj, r)
            if isnumeric(r)
                if numel(r)~=2 && numel(r)~=3
                    error('Covariance should be defined for 2D or 3D models')
                end
                obj.range = max(r, 100*eps); % range of 0 not allowed
            else
                error('range must be numeric')
            end
        end
        function set.angle(obj, a)
            if isnumeric(a)
                obj.angle = a;
            else
                error('angle must be numeric')
            end
        end
        function set.sill(obj, s)
            if isscalar(s)
                obj.sill = s;
            else
                error('sill must be a scalar')
            end
        end
    end
    
    methods (Access=protected)
        function [cx,rot] = trans(obj, cx)
            d = size(cx,2);
            rot = [];
            if d~=numel(obj.range)
                error('Dimensionality of input data inconsistent')
            end
            if d>1
                
                if d==2,
                    cang=cos(obj.angle/180*pi); sang=sin(obj.angle/180*pi);
                    rot=[cang,-sang;sang,cang];
                else
                     % rotation matrix in 3-D is computed around z, y and x
                     % in that order

                     cangz=cos(obj.angle(3)/180*pi);
                     sangz=sin(obj.angle(3)/180*pi);
                     cangy=cos(obj.angle(2)/180*pi);
                     sangy=sin(obj.angle(2)/180*pi);
                     cangx=cos(obj.angle(1)/180*pi);
                     sangx=sin(obj.angle(1)/180*pi);
                     rotz=[cangz,-sangz,0;sangz,cangz,0;0 0 1];
                     roty=[cangy,0,sangy;0 1 0;-sangy,0,cangy];
                     rotx=[1 0 0;0 cangx -sangx;0 sangx cangx];
                     rot=rotz*roty*rotx;
                end

                  % rotation is performed around z, y and x in that order,
                  %   the other coordinates are left unchanged.

                  cx=cx*rot;
                  t=diag(obj.range);

            else
                t = eye(d)*obj.range;
            end
            % perform contractions or dilatations (reduced h)
            cx=cx/t;
            
        end
        function h = compute_h(obj, x, x0)
            [n1,d]=size(x); % d dimension de l'espace
            [n2,d2]=size(x0);
            if d~=d2
                disp('Error: x and x0 should have the same dimensionality')
                return
            end            
            t1 = obj.trans(x);
            t2 = obj.trans(x0);
            h=0;
            for id=1:d
                h=h+(repmat(t1(:,id),1,n2)-repmat(t2(:,id)',n1,1)).^2;
            end
            h=sqrt(h);
        end
    end
    
    methods (Abstract)
        k = compute(obj, x, x0);
    end
    
end

