classdef Grid2D < Grid

    % when building vectors from 2D grids, Z is the "fast" axis, i.e.
    % indices are incremented for Z first: ind = (ix-1)*nz + iz
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

            % check if we have to raytrace in anisotropic media
            % possible arguments are:
            %   slowness, Tx, Rx, t0              -> isotropic
            %   slowness, xi, Tx, Rx, t0          -> elliptical anisotropy
            %   slowness, xi, theta, Tx, Rx, t0   -> tilted elliptical
            % in all cases t0 is optional
            if size(varargin{1}) ~= size(varargin{2})
                % second argument not same size as slowness: must be Tx
                type = 'iso';
            elseif size(varargin{1}) ~= size(varargin{3})
                % second arg same size as slowness but third not equal
                type = 'elliptical';
            else
                type = 'tilted';
            end


            if isempty(obj.mexObj)
                s.xmin = obj.grx(1);
                s.zmin = obj.grz(1);
                s.dx = obj.grx(2)-obj.grx(1);
                s.dz = obj.grz(2)-obj.grz(1);
                s.nx = length(obj.grx)-1;
                s.nz = length(obj.grz)-1;
                s.nsnx = obj.nsnx;
                s.nsnz = obj.nsnz;
                obj.mexObj = grid2d_mex('new', s, type, obj.nthreads);
            else
                % check that mexObj is consistent with current obj values
                % (in bh_tomo, this should not happen)
                if abs(grid2d_mex('get_xmin', obj.mexObj)-obj.grx(1))>100*eps || ...
                        abs(grid2d_mex('get_zmin', obj.mexObj)-obj.grz(1))>100*eps || ...
                        abs(grid2d_mex('get_dx', obj.mexObj)-(obj.grx(2)-obj.grx(1)))>100*eps || ...
                        abs(grid2d_mex('get_dz', obj.mexObj)-(obj.grz(2)-obj.grz(1)))>100*eps || ...
                        abs(grid2d_mex('get_nx', obj.mexObj)-(length(obj.grx)-1))>100*eps || ...
                        abs(grid2d_mex('get_nz', obj.mexObj)-(length(obj.grz)-1))>100*eps || ...
                        strcmp(grid2d_mex('get_type', obj.mexObj), type)~=1

                    % (in bh_tomo, we should not get here)

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
                    obj.mexObj = grid2d_mex('new', s, type, obj.nthreads);
                end
            end
            [varargout{1:nargout}] = grid2d_mex('raytrace', obj.mexObj, varargin{:});
        end

        function L = getForwardStraightRays(obj,varargin)
            % L = obj.getForwardStraightRays(ind,dx,dy,dz,aniso)
            aniso=false;
            ind = true(size(obj.Tx,1),1);
            grx = obj.grx;
            grz = obj.grz;
            if nargin>=2
                if ~isempty(varargin{1})
                    ind = varargin{1};
                end
            end
            if nargin>=3
                if ~isempty(varargin{2}) && isfinite(varargin{2})
                    dx = varargin{2};
                    grx = obj.grx(1):dx:obj.grx(end);
                end
            end
            % dy is ignored
            if nargin>=5
                if ~isempty(varargin{4}) && isfinite(varargin{4})
                    dz = varargin{4};
                    grz = obj.grz(1):dz:obj.grz(end);
                end
            end
            if nargin>=6
                aniso=varargin{5};
            end
            if aniso
                L = Lsr2da(obj.Tx(ind,[1 3]),obj.Rx(ind,[1 3]),grx,grz);
            else
                L = Lsr2d(obj.Tx(ind,[1 3]),obj.Rx(ind,[1 3]),grx,grz);
            end
        end
        function c = getCellCenter(obj,varargin)
            dx = obj.grx(2)-obj.grx(1);
            dz = obj.grz(2)-obj.grz(1);
            if nargin==4
                dx = varargin{1};
                % dy is ignored
                dz = varargin{3};
            end
            xmin = obj.grx(1)+dx/2;
            zmin = obj.grz(1)+dz/2;
            xmax = obj.grx(end)-dx/3;  % divide by 3 to avoid truncation error
            zmax = obj.grz(end)-dz/3;
            nx = ceil((xmax-xmin)/dx);
            nz = ceil((zmax-zmin)/dz);
            % z is the "fast" axis
            c=[kron((1:nx)',ones(nz,1)*dx), kron(ones(nx,1),(1:nz)'*dz)];
            c(:,1)=xmin+c(:,1)-dx;
            c(:,2)=zmin+c(:,2)-dz;
        end
        function c = checkCenter(obj,x,y,z)
            if ~isempty(y)
                c=0;
                return
            end
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

            if isstruct(cont) % check if we have anisotropy
                nc = obj.getNumberOfCells();
                nc1 = size(cont.data,1);
                nc2 = size(cont.data_xi,1);
                nc3 = size(cont.data_theta,1);
                indc = zeros(nc1+nc2+nc3,1);
                for i=1:size(cont.data,1)
                    ind11=findnear(cont.data(i,1),x(:,1));
                    ind22=findnear(cont.data(i,2),x(ind11,2));
                    indc(i)=min(ind11(ind22));
                end
                for i=1:size(cont.data_xi,1)
                    ind11=findnear(cont.data_xi(i,1),x(:,1));
                    ind22=findnear(cont.data_xi(i,2),x(ind11,2));
                    indc(nc1+i)=nc+min(ind11(ind22));
                end
                for i=1:size(cont.data_theta,1)
                    ind11=findnear(cont.data_theta(i,1),x(:,1));
                    ind22=findnear(cont.data_theta(i,2),x(ind11,2));
                    indc(nc1+nc2+i)=2*nc+min(ind11(ind22));
                end
            else
                indc = zeros(size(cont,1),1);
                for i=1:size(cont,1)
                    ind11=findnear(cont(i,1),x(:,1));
                    ind22=findnear(cont(i,2),x(ind11,2));
                    indc(i)=min(ind11(ind22));
                end
            end
        end
        function [Dx,Dy,Dz] = derivative(obj,order,varargin)
            % compute derivative operators for grid _cells_
            normalize=false;
            if nargin>=3
                % normalize to true cell size
                normalize=varargin{1};
            end
            dx=1;
            dz=1;
            if normalize
                dx=obj.dx;
                dz=obj.dz;
            end
            nx=length(obj.grx)-1;
            nz=length(obj.grz)-1;
            Dy = [];
            if order==1
                idx = 1/dx;
                idz = 1/dz;

                i = zeros(1,nz*nx*2);
                j = zeros(1,nz*nx*2);
                v = zeros(1,nz*nx*2);
                % forward operator is (u_{i+1} - u_i)/dx
                i(1:2*nz)=[1:nz 1:nz];
                j(1:2*nz)=[1:nz nz+(1:nz)];
                v(1:2*nz)=[-idx+zeros(1,nz) idx+zeros(1,nz)];
                % centered operator is (u_{i+1} - u_{i-1})/(2dx)
                for n=2:nx-1
                    i((1:2*nz)+(n-1)*2*nz) = (n-1)*nz+[1:nz 1:nz];
                    j((1:2*nz)+(n-1)*2*nz) = (n-1)*nz+[-nz+(1:nz) nz+(1:nz)];
                    v((1:2*nz)+(n-1)*2*nz) = 0.5*[-idx+zeros(1,nz) idx+zeros(1,nz)];
                end
                % backward operator is (u_i - u_{i-1})/dx
                i((1:2*nz)+(nx-1)*2*nz) = (nx-1)*nz+[1:nz 1:nz];
                j((1:2*nz)+(nx-1)*2*nz) = (nx-1)*nz+[-nz+(1:nz) 1:nz];
                v((1:2*nz)+(nx-1)*2*nz) = [-idx+zeros(1,nz) idx+zeros(1,nz)];
                Dx = sparse(i,j,v);

                i = kron(1:nx*nz,ones(1,2));
                jj = [[1 1:nz-1];[2 3:nz nz]];
                jj = reshape(jj,1,numel(jj));
                j = zeros(1,nz*nx*2);
                for n=1:nx
                    j((1:2*nz)+(n-1)*2*nz) = (n-1)*nz+jj;
                end
                v = [-idz idz repmat(0.5*[-idz idz],1,nz-2) -idz idz];
                v = kron(ones(1,nx),v);
                Dz = sparse(i,j,v);
            else
                % order==2
                idx2 = 1/(dx*dx);
                idz2 = 1/(dz*dz);

                i = zeros(1,nz*nx*3);
                j = zeros(1,nz*nx*3);
                v = zeros(1,nz*nx*3);
                % forward operator is (u_i - 2u_{i+1} + u_{i+2})/dx^2
                i(1:3*nz)=[1:nz 1:nz 1:nz];
                j(1:3*nz)=[1:nz nz+(1:nz) 2*nz+(1:nz)];
                v(1:3*nz)=[idx2+zeros(1,nz) -2*idx2+zeros(1,nz) idx2+zeros(1,nz)];
                % centered operator is (u_{i-1} - 2u_i + u_{i+1})/dx^2
                for n=2:nx-1
                    i((1:3*nz)+(n-1)*3*nz) = (n-1)*nz+[1:nz 1:nz 1:nz];
                    j((1:3*nz)+(n-1)*3*nz) = (n-1)*nz+[-nz+(1:nz) 1:nz nz+(1:nz)];
                    v((1:3*nz)+(n-1)*3*nz)=[idx2+zeros(1,nz) -2*idx2+zeros(1,nz) idx2+zeros(1,nz)];
                end
                % backward operator is (u_{i-2} - 2u_{i-1} + u_i)/dx^2
                i((1:3*nz)+(nx-1)*3*nz) = (nx-1)*nz+[1:nz 1:nz 1:nz];
                j((1:3*nz)+(nx-1)*3*nz) = (nx-1)*nz+[-2*nz+(1:nz) -nz+(1:nz) 1:nz];
                v((1:3*nz)+(nx-1)*3*nz)=[idx2+zeros(1,nz) -2*idx2+zeros(1,nz) idx2+zeros(1,nz)];
                Dx = sparse(i,j,v);

                i=kron(1:nx*nz,ones(1,3));
                jj=[[1 1:nz-2 nz-2];[2 2:nz-1 nz-1];[3 3:nz nz]];
                jj = reshape(jj,1,numel(jj));
                j = zeros(1,nz*nx*3);
                for n=1:nx
                    j((1:3*nz)+(n-1)*3*nz) = (n-1)*nz+jj;
                end
                v=kron(ones(1,nx*nz),idz2*[1 -2 1]);
                Dz = sparse(i,j,v);
            end
%             figure; imagesc(Dx); axis tight; axis equal; colorbar
%             figure; imagesc(Dz); axis tight; axis equal; colorbar
        end
        function G = preFFTMA(obj,cm)
            small = 1e-6;
            Nx = 2*length(obj.grx);
            Nz = 2*length(obj.grz);

            Nx2 = Nx/2;
            Nz2 = Nz/2;

            x = obj.dx*(0:Nx2-1);
            x = [x fliplr(-x)]';
            z = obj.dz*(0:Nz2-1);
            z = [z fliplr(-z)]';

            x = kron(x, ones(Nz,1));
            z = kron(ones(Nx,1), z);

            d = cm(1).compute([x z],[0 0]);
            for n=2:numel(cm)
                d = d + cm(n).compute([x z],[0 0]);
            end
            K = reshape(d,Nx,Nz)';  % transpose so that we get a Nz x Nx field

            mk=0;
            if min(K(:,1))>small
                % Enlarge grid to make sure that covariance falls to zero
                Nz=5*Nz;
                mk=1;
            end
            if min(K(1,:))>small
                Nx=5*Nx;
                mk=1;
            end
            if mk==1
                Nx2 = Nx/2;
                Nz2 = Nz/2;

                x = obj.dx*(0:Nx2-1);
                x = [x fliplr(-x)]';
                z = obj.dz*(0:Nz2-1);
                z = [z fliplr(-z)]';

                x = kron(ones(Nz,1), x);
                z = kron(z, ones(Nx,1));

                d = cm(1).compute([x z],[0 0]);
                for n=2:numel(cm)
                    d = d + cm(n).compute([x z],[0 0]);
                end
                K = reshape(d,Nx,Nz)';
            end
            G=fft2(K).^0.5;
        end
        function ms = FFTMA(obj,G)
            [Nz,Nx] = size(G);
            U=fft2(randn(size(G)));
            GU=G.*U;
            % Transformation de Fourier inverse donnant g*u et z
            Z=real(ifft2(GU));
            ms = Z((Nz+2)/2+1:(Nz+2)/2+length(obj.grz)-1,(Nx+2)/2+1:(Nx+2)/2+length(obj.grx)-1);
        end
        function toXdmf(obj,field,fieldname,filename)
            % for some reason we have to permute x & z (which h5create and
            % h5write do already)

            nx=length(obj.grx)-1;
            nz=length(obj.grz)-1;
            ox=obj.grx(1)+obj.dx/2;
            oz=obj.grz(1)+obj.dz/2;

            fid = fopen(filename,'wt');
            fprintf(fid,'<?xml version="1.0" ?>\n');
            fprintf(fid,'<!DOCTYPE Xdmf SYSTEM "Xdmf.dtd" []>\n');
            fprintf(fid,'<Xdmf xmlns:xi="http://www.w3.org/2003/XInclude" Version="2.2">\n');
%            fprintf(fid,' <Information Name="SampleLocation" Value="4"/>\n');
            fprintf(fid,' <Domain>\n');
            fprintf(fid,'   <Grid Name="Structured Grid" GridType="Uniform">\n');
            fprintf(fid,'     <Topology TopologyType="2DCORECTMesh" NumberOfElements="%d %d"/>\n',nz+1,nx+1);
            fprintf(fid,'     <Geometry GeometryType="ORIGIN_DXDY">\n');
            fprintf(fid,'       <DataItem Dimensions="2 " NumberType="Float" Precision="4" Format="XML">\n');
            fprintf(fid,'          %f %f\n',oz,ox);
            fprintf(fid,'       </DataItem>\n');
            fprintf(fid,'       <DataItem Dimensions="2 " NumberType="Float" Precision="4" Format="XML">\n');
            fprintf(fid,'        %f %f\n',obj.dz,obj.dx);
            fprintf(fid,'       </DataItem>\n');
            fprintf(fid,'     </Geometry>\n');
            fprintf(fid,'     <Attribute Name="%s" AttributeType="Scalar" Center="Cell">\n',fieldname);
            fprintf(fid,'       <DataItem Dimensions="%d %d" NumberType="Float" Precision="4" Format="HDF">%s.h5:/%s</DataItem>\n',nz,nx,filename,fieldname);
            fprintf(fid,'     </Attribute>\n');
            fprintf(fid,'   </Grid>\n');
            fprintf(fid,' </Domain>\n');
            fprintf(fid,'</Xdmf>\n');
            fclose(fid);

            if exist([filename,'.h5'],'file')
                delete([filename,'.h5'])
            end
            field=reshape(single(field),nz,nx)'; % z is fast dim, we have
            % to transpose to get x,z rather than z,x (which is needed by
            % h5create & h5write)
            h5create([filename,'.h5'],['/',fieldname],size(field),'Datatype','single');
            h5write([filename,'.h5'],['/',fieldname],field);

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
