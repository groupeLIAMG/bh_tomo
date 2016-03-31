classdef Grid2D < Grid
    properties
        flip
        borehole_x0
        x0
    end
    properties (Access = private, Hidden = true)
        mexObj   % Handle to the underlying C++ class instance
        nthreads
        nsnx
        nsnz
    end
    methods
        % Constructor
        function obj = Grid2D(grx,grz,varargin)
            if nargin>0
                obj.grx = grx(:);
                obj.grz = grz(:);
            end
            if nargin>=3
                obj.nthreads = varargin{1};
            else
                obj.nthreads = 1;
            end
            obj.nsnx = 10;
            obj.nsnz = 10;
        end
    
        % Destructor - Destroy the C++ class instance
        function delete(obj)
            if ~isempty(obj.mexObj)
                grid2d_mex('delete', obj.mexObj);
            end
        end
        
        % raytrace
        function varargout = raytrace(obj, varargin)
            if isempty(obj.mexObj)
                s.xmin = obj.grx(1);
                s.zmin = obj.grz(1);
                s.dx = obj.grx(2)-obj.grx(1);
                s.dz = obj.grz(2)-obj.grz(1);
                s.nx = length(obj.grx)-1;
                s.nz = length(obj.grz)-1;
                s.nsnx = obj.nsnx;
                s.nsnz = obj.nsnz;
                obj.mexObj = grid2d_mex('new', s, obj.nthreads);
            else
                % check that mexObj is consistent with current obj values
                if abs(grid2d_mex('get_xmin', obj.mexObj)-obj.grx(1))>10*eps || ...
                        abs(grid2d_mex('get_zmin', obj.mexObj)-obj.grz(1))>10*eps || ...
                        abs(grid2d_mex('get_dx', obj.mexObj)-(obj.grx(2)-obj.grx(1)))>10*eps || ...
                        abs(grid2d_mex('get_dz', obj.mexObj)-(obj.grz(2)-obj.grz(1)))>10*eps || ...
                        abs(grid2d_mex('get_nx', obj.mexObj)-(length(obj.grx)-1))>10*eps || ...
                        abs(grid2d_mex('get_nz', obj.mexObj)-(length(obj.grz)-1))>10*eps
                    
                    % delete old instance
                    grid2d_mex('delete', obj.mexObj);
                    
                    % create new instance with right dimensions
                    s.xmin = obj.grx(1);
                    s.zmin = obj.grz(1);
                    s.dx = obj.grx(2)-obj.grx(1);
                    s.dz = obj.grz(2)-obj.grz(1);
                    s.nx = length(obj.grx)-1;
                    s.nz = length(obj.grz)-1;
                    s.nsnx = obj.nsnx;
                    s.nsnz = obj.nsnz;
                    obj.mexObj = grid2d_mex('new', s, obj.nthreads);
                end
            end
            [varargout{1:nargout}] = grid2d_mex('raytrace', obj.mexObj, varargin{:});
        end
        
        % for saving in mat-files
        function s = saveobj(obj)
            s.nthreads = obj.nthreads;
            s.nsnx = obj.nsnx;
            s.nsnz = obj.nsnz;
            s.grx = obj.grx;
            s.grz = obj.grz;
            s.cont = obj.cont;
            s.Tx = obj.Tx;
            s.Rx = obj.Rx;
            s.TxCosDir = obj.TxCosDir;
            s.RxCosDir = obj.RxCosDir;
            s.x0 = obj.x0;
            s.flip = obj.flip;
            s.borehole_x0 = obj.borehole_x0;
            s.bord = obj.bord;
            s.Tx_Z_water = obj.Tx_Z_water;
            s.Rx_Z_water = obj.Rx_Z_water;
            s.in = obj.in;
            s.type = obj.type;
        end
    end
    methods(Static)
        % for loading from mat-files
        function obj = loadobj(s)
            if isstruct(s)
                obj = Grid2D(s.grx, s.grz, s.nthreads);
                
                obj.nsnx = s.nsnx;
                obj.nsnz = s.nsnz;
                obj.cont = s.cont;
                obj.Tx = s.Tx;
                obj.Rx = s.Rx;
                obj.TxCosDir = s.TxCosDir;
                obj.RxCosDir = s.RxCosDir;
                obj.x0 = s.x0;
                obj.flip = s.flip;
                obj.borehole_x0 = s.borehole_x0;
                obj.bord = s.bord;
                obj.Tx_Z_water = s.Tx_Z_water;
                obj.Rx_Z_water = s.Rx_Z_water;
                obj.in = s.in;
                obj.type = s.type;
            else
                error('Wrong input arguments')
            end
        end
    end    
end

