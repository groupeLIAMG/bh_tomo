function s = getAntSlown(c, k_form, diam, water)
%function s = getAntSlown(c, k_form, diam, water)
%
% retourne la lenteur le long de l'antenne
%  input
%       c : objet corrAntL
%       k_form : permttivite de la formation
%       diam : diametre du trou de forage
%       water : 1 si antenne dans l'eau, 0 si dans l'air
%

% Copyright (C) 2008 Bernard Giroux
%
% This program is free software; you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation; either version 2 of the License, or
% (at your option) any later version.
%
% This program is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details.
% 
% You should have received a copy of the GNU General Public License
% along with this program; if not, write to the Free Software
% Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA 02110-1301 USA
% 

if water
	s = interp2(c.k_formation', c.hole_diam, c.water_data, k_form, diam);
else
	s = interp2(c.k_formation', c.hole_diam, c.air_data, k_form, diam);
end
