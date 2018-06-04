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
    properties
        rotation
    end
    properties (Access = private, Hidden = true)
        mexObj   % Handle to the underlying C++ class instance
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
                if abs(grid3d_mex('get_xmin', obj.mexObj)-obj.grx(1))>100*eps || ...
                        abs(grid3d_mex('get_ymin', obj.mexObj)-obj.gry(1))>100*eps || ...
                        abs(grid3d_mex('get_zmin', obj.mexObj)-obj.grz(1))>100*eps || ...
                        abs(grid3d_mex('get_dx', obj.mexObj)-(obj.grx(2)-obj.grx(1)))>100*eps || ...
                        abs(grid3d_mex('get_dy', obj.mexObj)-(obj.gry(2)-obj.gry(1)))>100*eps || ...
                        abs(grid3d_mex('get_dz', obj.mexObj)-(obj.grz(2)-obj.grz(1)))>100*eps || ...
                        abs(grid3d_mex('get_nx', obj.mexObj)-(length(obj.grx)-1))>100*eps || ...
                        abs(grid3d_mex('get_ny', obj.mexObj)-(length(obj.gry)-1))>100*eps || ...
                        abs(grid3d_mex('get_nz', obj.mexObj)-(length(obj.grz)-1))>100*eps || ...
                        grid3d_mex('get_nthreads', obj.mexObj) ~= obj.nthreads

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
            grx = obj.grx;
            gry = obj.gry;
            grz = obj.grz;
            if nargin>=2
                ind = varargin{1};
            end
            if nargin>=3
                if ~isempty(varargin{2}) && isfinite(varargin{2})
                    dx = varargin{2};
                    grx = obj.grx(1):dx:obj.grx(end);
                end
            end
            if nargin>=4
                if ~isempty(varargin{3}) && isfinite(varargin{3})
                    dy = varargin{3};
                    gry = obj.gry(1):dy:obj.gry(end);
                end
            end
            if nargin>=5
                if ~isempty(varargin{4}) && isfinite(varargin{4})
                    dz = varargin{4};
                    grz = obj.grz(1):dz:obj.grz(end);
                end
            end
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
        function c = checkCenter(obj,x,y,z)
            dx = obj.grx(2)-obj.grx(1);
            xmin = obj.grx(1)+dx/2;
            xmax = obj.grx(end)-dx/3;  % divide by 3 to avoid truncation error
            xx=xmin:dx:xmax;
            if numel(x)~=numel(xx)
                c=0;
                return
            end
            if any(abs(x-reshape(xx,size(x)))>1000*eps)
                c=0;
                return
            end

            dy = obj.gry(2)-obj.gry(1);
            ymin = obj.gry(1)+dy/2;
            ymax = obj.gry(end)-dy/3;  % divide by 3 to avoid truncation error
            yy=ymin:dy:ymax;
            if numel(y)~=numel(yy)
                c=0;
                return
            end
            if any(abs(y-reshape(yy,size(y)))>1000*eps)
                c=0;
                return
            end

            dz = obj.grz(2)-obj.grz(1);
            zmin = obj.grz(1)+dz/2;
            zmax = obj.grz(end)-dz/3;
            zz=zmin:dz:zmax;
            if numel(z)~=numel(zz)
                c=0;
                return
            end
            if any(abs(z-reshape(zz,size(z)))>1000*eps)
                c=0;
                return
            end
            c=1;
        end

        function indc = getContIndices(obj,cont,varargin)
            if nargin==3
                x = varargin{1};
            else
                x = obj.getCellCenter();
            end
            indc = zeros(length(cont(:,1)), 1);
            for i=1:length(cont(:,1))
                ind11=findnear(cont(i,1),x(:,1));
                ind22=findnear(cont(i,2),x(ind11,2));
                ind33=findnear(cont(i,3),x(ind11(ind22),3));
                indc(i)=min(ind11(ind22(ind33)));
            end
        end
        function [Dx,Dy,Dz] = derivative(obj,order,varargin)
            % compute derivative operators for grid _cells_

            % x is "fastest" dimension, z is "slowest"

            normalize=false;
            if nargin>=3
                % normalize to true cell size
                normalize=varargin{1};
            end
            dx=1;
            dy=1;
            dz=1;
            if normalize
                dx=obj.dx;
                dy=obj.dy;
                dz=obj.dz;
            end
            nx=length(obj.grx)-1;
            ny=length(obj.gry)-1;
            nz=length(obj.grz)-1;
            if order==1
                idx = 1/dx;
                idy = 1/dy;
                idz = 1/dz;

                i = kron(1:nx*ny*nz,ones(1,2));
                jj = [[1 1:nx-1];[2 3:nx nx]];
                jj = reshape(jj,1,numel(jj));
                j = zeros(1,nz*ny*nx*2);
                for n=1:nz*ny
                    j((1:2*nx)+(n-1)*2*nx) = (n-1)*nx+jj;
                end
                v = [-idx idx repmat(0.5*[-idx idx],1,nx-2) -idx idx];
                v = kron(ones(1,ny*nz),v);
                Dx = sparse(i,j,v);


                nxy=nx*ny;
                i = zeros(1,nz*nxy*2);
                j = zeros(1,nz*nxy*2);
                v = zeros(1,nz*nxy*2);
                i(1:2*nx) = [1:nx 1:nx];
                j(1:2*nx) = [1:nx nx+(1:nx)];
                v(1:2*nx) = [-idy+zeros(1,nx) idy+zeros(1,nx)];
                for n=2:ny-1
                    i((1:2*nx)+(n-1)*2*nx) = (n-1)*nx+[1:nx 1:nx];
                    j((1:2*nx)+(n-1)*2*nx) = (n-1)*nx+[-nx+(1:nx) nx+(1:nx)];
                    v((1:2*nx)+(n-1)*2*nx) = 0.5*[-idy+zeros(1,nx) idy+zeros(1,nx)];
                end
                i((1:2*nx)+(ny-1)*2*nx) = (ny-1)*nx+[1:nx 1:nx];
                j((1:2*nx)+(ny-1)*2*nx) = (ny-1)*nx+[-nx+(1:nx) 1:nx];
                v((1:2*nx)+(ny-1)*2*nx) = [-idy+zeros(1,nx) idy+zeros(1,nx)];
                for n=2:nz
                    i((1:2*nxy)+(n-1)*2*nxy) = (n-1)*nxy+i(1:2*nxy);
                    j((1:2*nxy)+(n-1)*2*nxy) = (n-1)*nxy+j(1:2*nxy);
                    v((1:2*nxy)+(n-1)*2*nxy) = v(1:2*nxy);
                end
                Dy = sparse(i,j,v);

                i = zeros(1,nz*nxy*2);
                j = zeros(1,nz*nxy*2);
                v = zeros(1,nz*nxy*2);
                % forward operator is (u_{i+1} - u_i)/dx
                i(1:2*nxy) = [1:nxy 1:nxy];
                j(1:2*nxy) = [1:nxy nxy+(1:nxy)];
                v(1:2*nxy) = [-idz+zeros(1,nxy) idz+zeros(1,nxy)];
                % centered operator is (u_{i+1} - u_{i-1})/(2dx)
                for n=2:nz-1
                    i((1:2*nxy)+(n-1)*2*nxy) = (n-1)*nxy+[1:nxy 1:nxy];
                    j((1:2*nxy)+(n-1)*2*nxy) = (n-1)*nxy+[-nxy+(1:nxy) nxy+(1:nxy)];
                    v((1:2*nxy)+(n-1)*2*nxy) = 0.5*[-idz+zeros(1,nxy) idz+zeros(1,nxy)];
                end
                % backward operator is (u_i - u_{i-1})/dx
                i((1:2*nxy)+(nz-1)*2*nxy) = (nz-1)*nxy+[1:nxy 1:nxy];
                j((1:2*nxy)+(nz-1)*2*nxy) = (nz-1)*nxy+[-nxy+(1:nxy) 1:nxy];
                v((1:2*nxy)+(nz-1)*2*nxy) = [-idz+zeros(1,nxy) idz+zeros(1,nxy)];
                Dz = sparse(i,j,v);
            else
                % order = 2
                idx2 = 1/(dx*dx);
                idy2 = 1/(dy*dy);
                idz2 = 1/(dz*dz);

                i = kron(1:nx*ny*nz,ones(1,3));
                jj = [[1 1:nx-2 nx-2];[2 2:nx-1 nx-1];[3 3:nx nx]];
                jj = reshape(jj,1,numel(jj));
                j = zeros(1,nz*ny*nx*3);
                for n=1:nz*ny
                    j((1:3*nx)+(n-1)*3*nx) = (n-1)*nx+jj;
                end
                v = kron(ones(1,nx*ny*nz),idx2*[1 -2 1]);
                Dx = sparse(i,j,v);


                nxy=nx*ny;
                i = zeros(1,nz*nxy*3);
                j = zeros(1,nz*nxy*3);
                v = zeros(1,nz*nxy*3);
                i(1:3*nx) = [1:nx 1:nx 1:nx];
                j(1:3*nx) = [1:nx nx+(1:nx) 2*nx+(1:nx)];
                v(1:3*nx) = [idy2+zeros(1,nx) -2*idy2+zeros(1,nx) idy2+zeros(1,nx)];
                for n=2:ny-1
                    i((1:3*nx)+(n-1)*3*nx) = (n-1)*nx+[1:nx 1:nx 1:nx];
                    j((1:3*nx)+(n-1)*3*nx) = (n-1)*nx+[-nx+(1:nx) 1:nx nx+(1:nx)];
                    v((1:3*nx)+(n-1)*3*nx) = [idy2+zeros(1,nx) -2*idy2+zeros(1,nx) idy2+zeros(1,nx)];
                end
                i((1:3*nx)+(ny-1)*3*nx) = (ny-1)*nx+[1:nx 1:nx 1:nx];
                j((1:3*nx)+(ny-1)*3*nx) = (ny-1)*nx+[-2*nx+(1:nx) -nx+(1:nx) 1:nx];
                v((1:3*nx)+(ny-1)*3*nx) = [idy2+zeros(1,nx) -2*idy2+zeros(1,nx) idy2+zeros(1,nx)];
                for n=2:nz
                    i((1:3*nxy)+(n-1)*3*nxy) = (n-1)*nxy+i(1:3*nxy);
                    j((1:3*nxy)+(n-1)*3*nxy) = (n-1)*nxy+j(1:3*nxy);
                    v((1:3*nxy)+(n-1)*3*nxy) = v(1:3*nxy);
                end
                Dy = sparse(i,j,v);


                i = zeros(1,nz*nxy*3);
                j = zeros(1,nz*nxy*3);
                v = zeros(1,nz*nxy*3);
                % forward operator is (u_{i+1} - u_i)/dx
                i(1:3*nxy)=[1:nxy 1:nxy 1:nxy];
                j(1:3*nxy)=[1:nxy nxy+(1:nxy) 2*nxy+(1:nxy)];
                v(1:3*nxy)=[idz2+zeros(1,nxy) -2*idz2+zeros(1,nxy) idz2+zeros(1,nxy)];
                % centered operator is (u_{i+1} - u_{i-1})/(2dx)
                for n=2:nz-1
                    i((1:3*nxy)+(n-1)*3*nxy) = (n-1)*nxy+[1:nxy 1:nxy 1:nxy];
                    j((1:3*nxy)+(n-1)*3*nxy) = (n-1)*nxy+[-nxy+(1:nxy) 1:nxy nxy+(1:nxy)];
                    v((1:3*nxy)+(n-1)*3*nxy) = [idz2+zeros(1,nxy) -2*idz2+zeros(1,nxy) idz2+zeros(1,nxy)];
                end
                i((1:3*nxy)+(nz-1)*3*nxy) = (nz-1)*nxy+[1:nxy 1:nxy 1:nxy];
                j((1:3*nxy)+(nz-1)*3*nxy) = (nz-1)*nxy+[-2*nxy+(1:nxy) -nxy+(1:nxy) 1:nxy];
                v((1:3*nxy)+(nz-1)*3*nxy) = [idz2+zeros(1,nxy) -2*idz2+zeros(1,nxy) idz2+zeros(1,nxy)];
                Dz = sparse(i,j,v);
            end
        end
        function G = preFFTMA(obj,cm)
            small = 1e-6;
            Nx = 2*length(obj.grx);
            Ny = 2*length(obj.gry);
            Nz = 2*length(obj.grz);

            Nx2 = Nx/2;
            Ny2 = Ny/2;
            Nz2 = Nz/2;

            x = obj.dx*(0:Nx2-1);
            x = [x fliplr(-x)]';
            y = obj.dy*(0:Ny2-1);
            y = [y fliplr(-y)]';
            z = obj.dz*(0:Nz2-1);
            z = [z fliplr(-z)]';

            x = kron(ones(Ny*Nz,1), x);
            y = kron(kron(ones(Nz,1), y),ones(Nx,1));
            z = kron(z, ones(Nx*Ny,1));

            d = cm(1).compute([x y z],[0 0 0]);
            for n=2:numel(cm)
                d = d + cm(n).compute([x y z],[0 0 0]);
            end
            K = reshape(d,Nx,Ny,Nz);

            mk=0;
            tmp = K(:,:,1);
            if min(tmp(:))>small
                Nz=5*Nz;
                mk=1;
            end
            tmp = K(:,1,:);
            if min(tmp(:))>small
                Ny=5*Ny;
                mk=1;
            end
            tmp = K(1,:,:);
            if min(tmp)>small
                Nx=5*Nx;
                mk=1;
            end

            if mk==1;
                Nx2 = Nx/2;
                Ny2 = Ny/2;
                Nz2 = Nz/2;

                x = obj.dx*(0:Nx2-1);
                x = [x fliplr(-x)]';
                y = obj.dy*(0:Ny2-1);
                y = [y fliplr(-y)]';
                z = obj.dz*(0:Nz2-1);
                z = [z fliplr(-z)]';

                x = kron(ones(Ny*Nz,1), x);
                y = kron(kron(ones(Nz,1), y),ones(Nx,1));
                z = kron(z, ones(Nx*Ny,1));

                d = cm(1).compute([x y z],[0 0 0]);
                for n=2:numel(cm)
                    d = d + cm(n).compute([x y z],[0 0 0]);
                end
                K = reshape(d,Nx,Ny,Nz);
            end
            G=fftn(K).^0.5;
        end
        function ms = FFTMA(obj,G)
            [Nx,Ny,Nz] = size(G);
            U=fftn(randn(size(G)));
            GU=G.*U;
            % Transformation de Fourier inverse donnant g*u et z
            Z=real(ifftn(GU));
            ms = Z((Nx+2)/2+1:(Nx+2)/2+length(obj.grx)-1,...
                (Ny+2)/2+1:(Ny+2)/2+length(obj.gry)-1,...
                (Nz+2)/2+1:(Nz+2)/2+length(obj.grz)-1);
        end
        function toXdmf(obj,field,fieldname,filename)
            % for some reason we have to permute x & z (which h5create and
            % h5write do already)

            nx=length(obj.grx)-1;
            ny=length(obj.gry)-1;
            nz=length(obj.grz)-1;
            ox=obj.grx(1)+obj.dx/2;
            oy=obj.gry(1)+obj.dy/2;
            oz=obj.grz(1)+obj.dz/2;

            fid = fopen([filename,'.xmf'],'wt');
            fprintf(fid,'<?xml version="1.0" ?>\n');
            fprintf(fid,'<!DOCTYPE Xdmf SYSTEM "Xdmf.dtd" []>\n');
            fprintf(fid,'<Xdmf xmlns:xi="http://www.w3.org/2003/XInclude" Version="2.2">\n');
%            fprintf(fid,' <Information Name="SampleLocation" Value="4"/>\n');
            fprintf(fid,' <Domain>\n');
            fprintf(fid,'   <Grid Name="Structured Grid" GridType="Uniform">\n');
            fprintf(fid,'     <Topology TopologyType="3DCORECTMesh" NumberOfElements="%d %d %d "/>\n',nz+1,ny+1,nx+1);
            fprintf(fid,'     <Geometry GeometryType="ORIGIN_DXDYDZ">\n');
            fprintf(fid,'       <DataItem Dimensions="3 " NumberType="Float" Precision="4" Format="XML">\n');
            fprintf(fid,'          %f %f %f\n',oz,oy,ox);
            fprintf(fid,'       </DataItem>\n');
            fprintf(fid,'       <DataItem Dimensions="3 " NumberType="Float" Precision="4" Format="XML">\n');
            fprintf(fid,'        %f %f %f\n',obj.dz,obj.dy,obj.dx);
            fprintf(fid,'       </DataItem>\n');
            fprintf(fid,'     </Geometry>\n');
            fprintf(fid,'     <Attribute Name="%s" AttributeType="Scalar" Center="Cell">\n',fieldname);
            fprintf(fid,'       <DataItem Dimensions="%d %d %d" NumberType="Float" Precision="4" Format="HDF">%s.h5:/%s</DataItem>\n',nz,ny,nx,filename,fieldname);
            fprintf(fid,'     </Attribute>\n');
            fprintf(fid,'   </Grid>\n');
            fprintf(fid,' </Domain>\n');
            fprintf(fid,'</Xdmf>\n');
            fclose(fid);

            if exist([filename,'.h5'],'file')
                delete([filename,'.h5'])
            end
            field=reshape(single(field),nx,ny,nz);
            h5create([filename,'.h5'],['/',fieldname],size(field),'Datatype','single');
            h5write([filename,'.h5'],['/',fieldname],field);

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
            s.x0 = obj.x0;
            s.borehole_x0 = obj.borehole_x0;
            s.rotation = obj.rotation;
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
                if isfield(s, 'x0')
                    obj.x0 = s.x0;
                else
                    obj.x0 = [0 0 0];
                end
                if isfield(s, 'borehole_x0')
                    obj.borehole_x0 = s.borehole_x0;
                else
                    obj.borehole_x0 = 0;
                end
                if isfield(s, 'rotation')
                    obj.rotation = s.rotation;
                else
                    obj.rotation = 0;
                end
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
