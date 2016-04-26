%GRID3D class to perform raytracing in 3D with the fast sweeping method
%       on a regular 3D rectilinear grid
%
%  Usage:
%
%  Create and destroy instance of class
%
%    g = Grid3D(grx, gry, grz, nthreads)
%    clear g
%
%   Input for instantiation
%    grx: coordinates of nodes along X
%    gry: coordinates of nodes along Y
%    grz: coordinates of nodes along Z
%    nthreads: number of threads (optional, default = 1)
%
%  Raytracing
%    [tt] = g.raytrace(s, Tx, Rx, t0)
%    [tt, rays] = g.raytrace(s, Tx, Rx, t0)
%    [tt, rays, L] = g.raytrace(s, Tx, Rx, t0)
%
%   Input
%    g: grid instance
%    s: slowness vector ( nSlowness by 1 )
%    Tx: source coordinates, nTx by 3
%          1st column contains X coordinates, 2nd contains Y coordinates,
%          3rd contains Z coordinates
%    Rx: receiver coordinates, nRx by 3
%          1st column contains X coordinates, 2nd contains Y coordinates,
%          3rd contains Z coordinates
%    t0: source epoch, nTx by 1
%          t0 is optional (0 if not given)
%
%    *** IMPORTANT: Tx or Rx should _not_ lie on (or close to) an external
%                   face of the grid when rays are needed ***
%    *** nTx must be equal to nRx, i.e. each row define one Tx-Rx pair ***
%    *** nSlowness must equal the number of _voxels_ of the grid ***
%    *** Indexing of slowness values is done by "vectorizing" a 3D array,
%        i.e. if slowness field s is of size (nx,ny,nz), enter s(:) as
%        first argument
%
%
%   Output
%    tt:   vector of traveltimes, nRx by 1
%    rays: cell object containing the matrices of coordinates of the ray
%          paths, nRx by 1.  Each matrix is nPts by 3
%    L:    data kernel matrix (tt = L*s)
%
% -----------
%
% Bernard Giroux
% INRS-ETE
% 2016-02-14

% Copyright (C) 2016 Bernard Giroux
%
% This program is free software; you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation; either version 3 of the License, or
% (at your option) any later version.
%
% This program is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details.
%
% You should have received a copy of the GNU General Public License
% along with this program.  If not, see <http://www.gnu.org/licenses/>.
%
%


classdef Grid3D < Grid
    properties (Access = private, Hidden = true)
        mexObj   % Handle to the underlying C++ class instance
        nthreads
    end
    methods
        % Constructor - Create a new C++ class instance
        function obj = Grid3D(grx,gry,grz,varargin)
            if nargin>0
                obj.grx = grx(:);
                obj.gry = gry(:);
                obj.grz = grz(:);
                dx=grx(2)-grx(1);
                dy=gry(2)-gry(1);
                dz=grz(2)-grz(1);
                if dx~=dy || dx~=dz || dy~=dz
                    error('Grid step size should be equal in all directions')
                end
            end
            if nargin>=4
                obj.nthreads = varargin{1};
            else
                obj.nthreads = 1;
            end
            
        end
        
        % Destructor - Destroy the C++ class instance
        function delete(obj)
            if ~isempty(obj.mexObj)
                grid3d_mex('delete', obj.mexObj);
            end
        end
        
        % raytrace
        function varargout = raytrace(obj, varargin)
            if isempty(obj.mexObj)
                s.xmin = obj.grx(1);
                s.ymin = obj.gry(1);
                s.zmin = obj.grz(1);
                s.dx = obj.grx(2)-obj.grx(1);
                s.dy = obj.gry(2)-obj.gry(1);
                s.dz = obj.grz(2)-obj.grz(1);
                s.nx = length(obj.grx)-1;
                s.ny = length(obj.gry)-1;
                s.nz = length(obj.grz)-1;
                obj.mexObj = grid3d_mex('new', s, obj.nthreads);
            else
                % check that mexObj is consistent with current obj values
                if abs(grid3d_mex('get_xmin', obj.mexObj)-obj.grx(1))>10*eps || ...
                        abs(grid3d_mex('get_ymin', obj.mexObj)-obj.gry(1))>10*eps || ...
                        abs(grid3d_mex('get_zmin', obj.mexObj)-obj.grz(1))>10*eps || ...
                        abs(grid3d_mex('get_dx', obj.mexObj)-(obj.grx(2)-obj.grx(1)))>10*eps || ...
                        abs(grid3d_mex('get_dy', obj.mexObj)-(obj.gry(2)-obj.gry(1)))>10*eps || ...
                        abs(grid3d_mex('get_dz', obj.mexObj)-(obj.grz(2)-obj.grz(1)))>10*eps || ...
                        abs(grid3d_mex('get_nx', obj.mexObj)-(length(obj.grx)-1))>10*eps || ...
                        abs(grid3d_mex('get_ny', obj.mexObj)-(length(obj.gry)-1))>10*eps || ...
                        abs(grid3d_mex('get_nz', obj.mexObj)-(length(obj.grz)-1))>10*eps
                    
                    % delete old instance
                    grid3d_mex('delete', obj.mexObj);
                    
                    % create new instance with right dimensions
                    s.xmin = obj.grx(1);
                    s.ymin = obj.gry(1);
                    s.zmin = obj.grz(1);
                    s.dx = obj.grx(2)-obj.grx(1);
                    s.dy = obj.gry(2)-obj.gry(1);
                    s.dz = obj.grz(2)-obj.grz(1);
                    s.nx = length(obj.grx)-1;
                    s.ny = length(obj.gry)-1;
                    s.nz = length(obj.grz)-1;
                    obj.mexObj = grid3d_mex('new', s, obj.nthreads);
                end
            end
            [varargout{1:nargout}] = grid3d_mex('raytrace', obj.mexObj, varargin{:});
        end
        
        function L = getForwardStraightRays(obj,varargin)
            ind = true(size(obj.Tx,1),1);
            dx = obj.grx(2)-obj.grx(1);
            dy = obj.gry(2)-obj.gry(1);
            dz = obj.grz(2)-obj.grz(1);
            if nargin>=2
                ind = varargin{1};
            end
            if nargin>=3
                if ~isempty(varargin{2}) && isfinite(varargin{2})
                    dx = varargin{2};
                end
            end
            if nargin>=4
                if ~isempty(varargin{3}) && isfinite(varargin{3})
                    dy = varargin{3};
                end
            end
            if nargin>=5
                if ~isempty(varargin{4}) && isfinite(varargin{4})
                    dz = varargin{4};
                end
            end
            grx = obj.grx(1):dx:obj.grx(end);
            gry = obj.gry(1):dy:obj.gry(end);
            grz = obj.grz(1):dz:obj.grz(end);
            L = Lsr3d(obj.Tx(ind,:),obj.Rx(ind,:),grx,gry,grz);
        end
        function c = getCellCenter(obj,varargin)
            dx = obj.grx(2)-obj.grx(1);
            dy = obj.gry(2)-obj.gry(1);
            dz = obj.grz(2)-obj.grz(1);
            if nargin==4
                dx = varargin{1};
                dy = varargin{2};
                dz = varargin{3};
            end
            xmin = obj.grx(1)+dx/2;
            ymin = obj.gry(1)+dy/2;
            zmin = obj.grz(1)+dz/2;
            xmax = obj.grx(end)-dx/3;  % divide by 3 to avoid truncation error
            ymax = obj.gry(end)-dy/3;
            zmax = obj.grz(end)-dz/3;
            nx = ceil((xmax-xmin)/dx);
            ny = ceil((ymax-ymin)/dy);
            nz = ceil((zmax-zmin)/dz);
            c=[kron(ones(ny*nz,1),(1:nx)'*dx) ...
                kron(kron(ones(nz,1),(1:ny)'*dy),ones(nx,1)) ...
                kron((1:nz)'*dz,ones(nx*ny,1))];
            c(:,1)=xmin+c(:,1)-dx;
            c(:,2)=ymin+c(:,2)-dy;
            c(:,3)=zmin+c(:,3)-dz;
        end
        % for saving in mat-files
        function s = saveobj(obj)
            s.nthreads = obj.nthreads;
            s.grx = obj.grx;
            s.gry = obj.gry;
            s.grz = obj.grz;
            s.cont = obj.cont;
            s.Tx = obj.Tx;
            s.Rx = obj.Rx;
            s.TxCosDir = obj.TxCosDir;
            s.RxCosDir = obj.RxCosDir;
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
                obj = Grid3D(s.grx, s.gry, s.grz, s.nthreads);
                
                obj.cont = s.cont;
                obj.Tx = s.Tx;
                obj.Rx = s.Rx;
                obj.TxCosDir = s.TxCosDir;
                obj.RxCosDir = s.RxCosDir;
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
