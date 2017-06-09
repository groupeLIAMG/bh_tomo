function indc = calcul_kss_aniso(x, covar, cont, use_cont, nc, varargin)
%

% Copyright (C) 2010 Bernard Giroux
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

global kss ks0 Css

handle = [];
if nargin>5
    handle = varargin{1};
end
if isempty(handle)
    disp('Calcul direct de la matrice de Covariance')
else
    str = get_str_locale();
    set(handle, 'String', [str.s225,' - ',str.s235]); drawnow
end
drawnow

kss_s  = covardm(x, x, covar.model,    covar.c);
kss_xi = covardm(x, x, covar.model_xi, covar.c_xi);
kss = [kss_s zeros(size(kss_s)); zeros(size(kss_xi)) kss_xi];

if covar.nugget_l ~= 0
    kss = kss + covar.nugget_l*eye(size(kss,1));
end

indc = [];
if use_cont
    if isfield(cont,'data_xi')  %%YH
    if ~isempty( cont.data ) && ~isempty( cont.data_xi )   %%YH || -> &&  
        %
        % selon x
        %
        ncx = size(cont.data,1);
        ncz = size(cont.data_xi,1);
        indc = zeros(ncx+ncz,1);
        for i=1:size(cont.data,1)
            ind11=findnear(cont.data(i,1),x(:,1));
            ind22=findnear(cont.data(i,2),x(ind11,2));
            indc(i)=min(ind11(ind22));
        end
        for i=1:size(cont.data_xi,1)
            ind11=findnear(cont.data_xi(i,1),x(:,1));
            ind22=findnear(cont.data_xi(i,2),x(ind11,2));
            indc(ncx+i)=nc+min(ind11(ind22));
        end
        clear ind11 ind22
    
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
        %Terme de gauche
        ks0=kss(indc,:);
        %Covariance des donnees de contrainte
    
        % contraintes dures
        Css = kss(indc,indc);
    
        % contraintes souples
        % la variance est contenue ds la 4e colonne
        if size(cont.data,2)==4 && size(cont.data_xi,2)==4
            Css = Css + diag([cont.data(:,4); cont.data_xi(:,4)]);
        end
    end
    end
end
