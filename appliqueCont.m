function s = appliqueCont(s, cont, varargin)
%
%
%

if nargin == 3
    g = varargin{1};
    grx = g.xmin + (0:g.nx)*g.dx;
    grz = g.zmin + (0:g.nz)*g.dz;
    xx = (grx(1)+g.dx/2):g.dx:(grx(end)-g.dx/2);
    zz = (grz(1)+g.dz/2):g.dz:(grz(end)-g.dz/2);

    if ~isempty(cont)
        if ischar(cont)  % cont est un fichier
            data = load(cont);
        else
            data = cont.data;
        end
        for n=1:size(data,1)
            xc=data(n,2);
            if xc<xx(1) || xc>xx(end)
                continue
            end
            zc=data(n,1);
            if zc<zz(1) || zc>zz(end)
                continue
            end
            
            [~, ic] = min( abs( xc - xx ) );
            [~, jc] = min( abs( zc - zz ) );
            s((ic-1)*g.nz + jc) = data(n,3);
        end
    end


elseif nargin == 5
    g = varargin{1};
    grx = g.grx;
    grz = g.grz;
    gridx = varargin{2};
    gridz = varargin{3};

    DX=grx(2)-grx(1);
    DZ=grz(2)-grz(1);

    xx = (grx(1)+DX/2):DX:(grx(end)-DX/2);
    zz = (grz(1)+DZ/2):DZ:(grz(end)-DZ/2);

    s=interp2(gridx, gridz', reshape(s,length(gridz),length(gridx)),xx,zz');
    if ~isempty(cont)
        if ischar(cont)  % cont est un fichier
            data = load(cont);
        else
            data = cont.data;
        end
        for n=1:size(data,1)
            xc=data(n,2);
            if xc<xx(1) || xc>xx(end)
                continue
            end
            zc=data(n,1);
            if zc<zz(1) || zc>zz(end)
                continue
            end
            ic=findnear(xc,xx);
            jc=findnear(zc,zz);
            s(jc,ic) = data(n,3);
        end
    end

else
    s = [];
    return
end

