function [dif]=fitModeliKss3d(x0libre,covar,id,L,x,iktt,afi,lclas,c0,ax1,ax2)
% [dif]=fitModeliKss3d(x0libre,covar,id,L,x,iktt,afi,lclas,c0,ax1,ax2)

% Copyright (C) 2013 Bernard Giroux
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


[g,gt] = modeliKss3d(x0libre,covar,id,L,x,iktt,afi,lclas,c0,ax1,ax2);

q = fliplr((1:length(g)).^2);
q = 1 + q./max(q);
%q=round(q(1:floor(length(q)/afi)));
%q=q(1:floor(length(q)/afi));
g = g(1:length(q));
gt = gt(1:length(q));

dif = sum((((gt-g).*q).^2));

