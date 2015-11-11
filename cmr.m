function c = cmr(m)
% A colormap for effective black and white rendering of color scale images
%   cmr(m) is a m-by-3 matrix containing the colormap.
%   cmr, by itself, is the same length as the current figure's
%   colormap. If no figure exists, MATLAB creates one.
%
% Adapted from the CMRmap.m of Carey Rappaport
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
CMRmap_v=[0 0 0;.15 .15 .5;.3 .15 .75;.6 .2 .50;1 .25 .15;.9 .5 0;.9 .75 .1;.9 .9 .5;1 1 1];
x = 1:8/(m-1):9;    % m color levels instead of 9
x1 = 1:9;
for i = 1:3
  c(:,i) = spline(x1,CMRmap_v(:,i),x)';
  % spline fit intermediate values
end
c = abs(c/max(max(c)));
