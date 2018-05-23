classdef Grid < matlab.mixin.Copyable
    %GRID Base class for 2D and 3D grids

    properties
        grx
        gry
        grz
        cont
        Tx
        Rx
        TxCosDir
        RxCosDir
        bord
        Tx_Z_water
        Rx_Z_water
        in
        type
    end
    methods
        function obj = Grid()
            obj.cont = Constraints();
        end
        function set.grx(obj,g)
            if ~isnumeric(g)
                error('Grid coordinates must be numeric')
            end
            if any(abs(diff(g)-(g(2)-g(1)))>100*eps)
                error('Grid step size should be constant')
            end
            obj.grx = g;
        end
        function set.gry(obj,g)
            if ~isnumeric(g)
                error('Grid coordinates must be numeric')
            end
            if any(abs(diff(g)-(g(2)-g(1)))>100*eps)
                error('Grid step size should be constant')
            end
            obj.gry = g;
        end
        function set.grz(obj,g)
            if ~isnumeric(g)
                error('Grid coordinates must be numeric')
            end
            if any(abs(diff(g)-(g(2)-g(1)))>100*eps)
                error('Grid step size should be constant')
            end
            obj.grz = g;
        end
        function nc = getNumberOfCells(obj)
            nc = (length(obj.grx)-1)*(length(obj.grz)-1);
            if ~isempty(obj.gry)
                nc = nc*(length(obj.gry)-1);
            end
        end
        function s = getNcell(obj)
            if ~isempty(obj.gry)
                s = [numel(obj.grx)-1 numel(obj.gry)-1 numel(obj.grz)-1];
            else
                s = [numel(obj.grz)-1 numel(obj.grx)-1];
            end
        end
        function d = dx(obj)
            d=obj.grx(2)-obj.grx(1);
        end
        function d = dy(obj)
            if isempty(obj.gry)
                d = 0;
            else
                d=obj.gry(2)-obj.gry(1);
            end
        end
        function d = dz(obj)
            d=obj.grz(2)-obj.grz(1);
        end
    end
    methods (Static)
        function [x0, a, d, normd] = lsplane(X)
            % ---------------------------------------------------------------------
            % LSPLANE.M Least-squares plane (orthogonal distance
            % regression).
            %
            % Version 1.0
            % Last amended I M Smith 27 May 2002.
            % Created I M Smith 08 Mar 2002
            % ---------------------------------------------------------------------
            % Input
            % X Array [x y z] where x = vector of x-coordinates,
            % y = vector of y-coordinates and z = vector of
            % z-coordinates.
            % Dimension: m x 3.
            %
            % Output
            % x0 Centroid of the data = point on the best-fit plane.
            % Dimension: 3 x 1.
            %
            % a Direction cosines of the normal to the best-fit
            % plane.
            % Dimension: 3 x 1.
            %
            % <Optional...
            % d Residuals.
            % Dimension: m x 1.
            %
            % normd Norm of residual errors.
            % Dimension: 1 x 1.
            % ...>
            %
            % [x0, a <, d, normd >] = lsplane(X)
            % ---------------------------------------------------------------------
            % check number of data points
            m = size(X, 1);
            if m < 3
                error('At least 3 data points required: ' )
            end
            %
            % calculate centroid
            x0 = mean(X);
            %
            % form matrix A of translated points
            A = [(X(:, 1) - x0(1)) (X(:, 2) - x0(2)) (X(:, 3) - x0(3))];
            %
            % calculate the SVD of A
            [U, S, V] = svd(A, 0);
            %
            % find the smallest singular value in S and extract from V the
            % corresponding right singular vector
            [s, i] = min(diag(S));
            a = V(:, i)';
            %
            % calculate residual distances, if required
            if nargout > 2
                d = U(:, i)*s;
                normd = norm(d);
            end
            % ---------------------------------------------------------------------
            % End of LSPLANE.M.
        end
        function order = boreholes_order( bh )
            % order = boreholes_order( bh )
            %
            %

            nd = length(bh);
            x = zeros(1,nd);
            y = x;

            for n=1:nd
                x(n) = bh(n).X;
                y(n) = bh(n).Y;
            end

            dx = max(x)-min(x);
            dy = max(y)-min(y);

            if ( dx > dy )
                v1 = x;
                v2 = y;
            else
                v1 = y;
                v2 = x;
            end

            % premier forage a x min
            order = 1:nd;
            [v1, iv1] = sort( v1 );
            order = order(iv1);
            v2 = v2(iv1);

            ind = v1==v1(1);
            % si x_min est repete
            if sum( ind ) > 1
                [~, iv2] = sort( v2(ind) );
                v2(ind) = v2( iv2 );
                order(ind) = order( iv2 );
            end

            for n=1:nd-2
                dist = sqrt( ( v1(n)-v1((n+1):nd) ).^2 + ( v2(n)-v2((n+1):nd) ).^2 );
                [~, ind] = sort(dist);
                tmp = v1((n+1):nd);
                v1((n+1):nd) = tmp(ind);
                tmp = v2((n+1):nd);
                v2((n+1):nd) = tmp(ind);
                tmp = order((n+1):nd);
                order((n+1):nd) = tmp(ind);
            end
        end
        function  p_data = proj_plane(data, x0, a)
            %  p_data = proj_plane(data, x0, a)
            %
            % data : coordonnees a projeter (nx3)
            % x0   : pt sur le plan (1x3)
            % a    : cosinus directeur de la normale au plan (1x3)
            %
            % p_data : coordonnees projetees
            %
            %
            p_data = data;
            for n=1:size(data,1)
                r = x0-data(n,:);     % vecteur pointant de data vers x0
                p = dot(a,r);         % distance entre data et le plan
                p_data(n,:) = data(n,:) + p*a;   % coord de data projete sur le plan
            end
        end
        function [p_data,no_plane] = proj_planes(data, planes)
            %  p_data = proj_planes(data, planes)
            %
            % data : coordonnees a projeter (nx3)
            % plans: vecteur des plans, contenant struct
            %        x0   : pt sur le plan (1x3)
            %        a    : cosinus directeur de la normale au plan (1x3)
            %
            % p_data : coordonnees projetees
            % no_plan : no du plan de projection
            %
            %
            p_data = data;
            no_plane = zeros(1,size(data,1));
            p = zeros(1,length(planes));
            d = zeros(1,length(planes));
            for n=1:size(data,1)
                for nn=1:length(planes)
                    r = planes(nn).x0-data(n,:);           % vecteur pointant de data vers x0
                    d(nn) = sqrt(sum(r.*r));               % distance entre data et le centroid
                    p(nn) = dot(planes(nn).a, r);          % distance entre data et le plan
                end
                [~,no] = min(abs(d));
                % on va garder le plan pour lequel la distance au centroid est la plus faible
                p_data(n,:) = data(n,:) + p(no)*planes(no).a;  % coord de data projete sur le plan
                no_plane(n) = no;
            end
        end
        function m_data = transl_rotat(data, origine, az, dip)
            % m_data = transl_rotat(data, origine, az, dip)

            % translation p/r origine
            m_data = data-repmat(origine,size(data,1),1);
            % rotation p/r azimuth
            if abs(az) > (pi/720)  % si plus grand que 1/4 degre
                rot = [cos(az) -sin(az); sin(az) cos(az)];
                for n=1:size(m_data,1)
                    m_data(n,1:2) = m_data(n,1:2)*rot';
                end
                % data(n,1:2)*rot' est egal a (rot*data(n,1:2)')'
            end

            % rotation p/r pendage
            if abs(dip) > (pi/720)
                rot = [cos(dip) -sin(dip); sin(dip) cos(dip)];
                for n=1:size(m_data,1)
                    m_data(n,2:3) = m_data(n,2:3)*rot';
                end
                % data(n,2:3)*rot' est egal a (rot*data(n,2:3)')'
            end
        end

    end

end
