function indc = calcul_kss_aniso2(x, covar, cont, use_cont, ncell, varargin)
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
    set(handle, 'String', [str.s225,' - ',str.s235])
end
drawnow

kss_s  = covardm(x, x, covar.model,    covar.c);
kss_xi = covardm(x, x, covar.model_xi, covar.c_xi);
kss_th = covardm(x, x, covar.model_th, covar.c_th);
kss = [kss_s zeros(size(kss_xi)) zeros(size(kss_th)) ;
       zeros(size(kss_s)) kss_xi zeros(size(kss_th)) ;
       zeros(size(kss_s)) zeros(size(kss_xi)) kss_th] ;

if covar.nugget_l ~= 0
    kss = kss + covar.nugget_l*eye(size(kss,1));
end

indc = [];
if use_cont && isfield(cont, 'data_xi') && isfield(cont, 'data_th')
    if ~isempty( cont.data ) || ~isempty( cont.data_xi ) || ...
            ~isempty( cont.data_th )
        
        %
        % selon x
        %
        ncx = size(cont.data,1);
        ncz = size(cont.data_xi,1);
        nct = size(cont.data_th,1);
        indc = zeros(ncx+ncz+nct,1);
        for i=1:ncx
            ind11=findnear(cont.data(i,1),x(:,1));
            ind22=findnear(cont.data(i,2),x(ind11,2));
            indc(i)=min(ind11(ind22));
        end
        for i=1:ncz
            ind11=findnear(cont.data_xi(i,1),x(:,1));
            ind22=findnear(cont.data_xi(i,2),x(ind11,2));
            indc(ncx+i)=ncell+min(ind11(ind22));
        end
        for i=1:nct
            ind11=findnear(cont.data_th(i,1),x(:,1));
            ind22=findnear(cont.data_th(i,2),x(ind11,2));
            indc(ncx+ncz+i)=2*ncell+min(ind11(ind22));
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
        tmp = [];
        if size(cont.data,2)==4 
            tmp = [tmp; cont.data(:,4)];
        else
            tmp = [tmp; zeros(ncx,1)];
        end
        if size(cont.data_xi,2)==4
            tmp = [tmp; cont.data_xi(:,4)];
        else
            tmp = [tmp; zeros(ncz,1)];
        end
        if size(cont.data_th,2)==4
            tmp = [tmp; cont.data_th(:,4)];
        else
            tmp = [tmp; zeros(nct)];
        end
        Css = Css + diag(tmp);
        
    end
end
