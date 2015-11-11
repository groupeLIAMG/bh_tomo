function [np]=calcul_kss_global3D(x,indi,covar,varargin)
% [np]=calcul_kss_global3D(x,indi,covar,varargin)
% kss is global
  
% Copyright (C) 2005 Erwan Gloaguen
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


global kss
handle = [];
if nargin>3
	handle = varargin{1};
end

if isempty(handle)
	disp('Calcul direct de la matrice de Covariance')
else
	str = get_str_locale();
	set(handle, 'String', [str.s225,' - ',str.s235]); drawnow
end
x=x(indi,:);
np=length(x(:,1));
kss=covardm(x,x,covar.model,covar.c);

if covar.nugget_l ~= 0
	kss = kss + covar.nugget_l*eye(np);
end
