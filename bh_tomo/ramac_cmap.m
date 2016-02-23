function c = ramac_cmap(m)
% A colormap based on the RAMAC GPR program colormap
%   by Bernard Giroux, Ecole Polytechnique de Montreal

% Copyright (C) 2006 Bernard Giroux
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

if nargin < 1
   m = size(get(gcf,'colormap'),1);
end
CMRmap_v=[0.00 0.00 0.00;
    0.00 0.00 0.68;
    0.00 0.67 0.00;
    0.00 0.67 0.68;
    0.68 0.00 0.00;
    0.68 0.00 0.68;
    0.68 0.33 0.00;
    0.68 0.67 0.68;
    0.32 0.33 0.32;
    0.32 0.33 1.00;
    0.32 1.00 0.32;
    0.32 1.00 1.00;
    1.00 0.33 0.32;
    1.00 0.33 1.00;
    1.00 1.00 0.32;
    1.00 1.00 1.00];
x = 1:15/(m-1):16;    % m color levels instead of 16
x1 = 1:16;
for i = 1:3
  c(:,i) = spline(x1,CMRmap_v(:,i),x)';
  % spline fit intermediate values
end
c = abs(c/max(max(c)));
