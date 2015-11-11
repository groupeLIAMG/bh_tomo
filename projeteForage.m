function [x, y, z, c] = projeteForage(fdata, prof, nom)
% function [x, y, z, c] = projeteForage(fdata, prof, nom)
%
% fdata: tableau [n x 3]
%        coordonnées (x,y,elev) de la trajectoire du
%        forage, classées du collet vers le fond du trou
% prof:  tableau [1 x n2]
%        profondeur des n2 pts de mesure le long du trou à partir du collet
% nom:   nom du forage
%
% x:     coordonnées en x des pts de mesure [1 x n2]
% y:     coordonnées en y des pts de mesure [1 x n2]
% z:     élévation des pts de mesure [1 x n2]
% c:     cosinus directeurs du forage aux pts de mesure [n2 x 3], pointent
%           vers le fond
%

% Copyright (C) 2007 Bernard Giroux
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
npts = length(prof);
x = zeros(1,npts);
y = zeros(1,npts);
z = zeros(1,npts);
c = zeros(npts,3);

profForage = [0; cumsum(sqrt(sum( (diff(fdata,1)).^2, 2)))];

for n=1:npts
	% indice du début du segment
	i1 = find(prof(n)>=profForage);
	if isempty(i1)
		str=get_str_locale();
		errordlg([str.s271, ': ',nom])
		x = zeros(1,npts);
		y = zeros(1,npts);
		z = zeros(1,npts);
		c = zeros(npts,3);
		return
	end
	i1 = i1(end);
	% indice de fin du segment
	i2 = find(prof(n)<profForage);
	if isempty(i2)
		str=get_str_locale();
		errordlg([str.s271, ': ',nom])	
		x = zeros(1,npts);
		y = zeros(1,npts);
		z = zeros(1,npts);
		c = zeros(npts,3);
		return
	end
	i2 = i2(1);
	
	d = sqrt(sum((fdata(i2,1:3)-fdata(i1,1:3)).^2));
	% cosinus directeurs
	l = (fdata(i2,1:3)-fdata(i1,1:3))./d;
	
	d2 = prof(n) - profForage(i1);
	
	x(n) = fdata(i1,1) + d2*l(1);
	y(n) = fdata(i1,2) + d2*l(2);
	z(n) = fdata(i1,3) + d2*l(3);
	c(n,:) = l;
end
