function p_data = proj_plan(data, x0, a)
%  p_data = proj_plan(data, x0, a)
%
% data : coordonnees a projeter (nx3)
% x0   : pt sur le plan (1x3)
% a    : cosinus directeur de la normale au plan (1x3)
%
% p_data : coordonnees projetees
%
%
  
% Copyright (C) 2005 Bernard Giroux
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


p_data = data;
for n=1:size(data,1)
    r = x0-data(n,:);     % vecteur pointant de data vers x0
    p = dot(a,r);         % distance entre data et le plan
    p_data(n,:) = data(n,:) + p*a;   % coord de data projete sur le plan
end
