classdef Grid2D < Grid
    properties (Access = private, Hidden = true)
        mexObj   % Handle to the underlying C++ class instance
        nthreads
    end
    methods
        % Constructor
        function obj = Grid3D(grx,grz,varargin)
            obj.grx = grx(:);
            obj.grz = grz(:);
            if nargin>=3
                obj.nthreads = varargin{1};
            else
                obj.nthreads = 1;
            end
            
        end
    end
    
end

