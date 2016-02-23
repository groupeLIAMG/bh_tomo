function order = boreholes_order( bh )
% order = boreholes_order( bh )
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
	[tmp, iv2] = sort( v2(ind) );
	v2(ind) = v2( iv2 );
	order(ind) = order( iv2 );
end

for n=1:nd-2
    dist = sqrt( ( v1(n)-v1((n+1):nd) ).^2 + ( v2(n)-v2((n+1):nd) ).^2 );
	[tmp, ind] = sort(dist);
	tmp = v1((n+1):nd);
	v1((n+1):nd) = tmp(ind);
	tmp = v2((n+1):nd);
	v2((n+1):nd) = tmp(ind);
	tmp = order((n+1):nd);
	order((n+1):nd) = tmp(ind);
end
