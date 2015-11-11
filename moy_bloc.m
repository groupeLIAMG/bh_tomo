function m = moy_bloc(xy,lclas)
% m = moy_bloc(xy,lclas)

% Copyright (C) 2005 Erwan Gloaguen, Bernard Giroux
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


k = floor(length(xy)/lclas);

%m = zeros(1,k);
%for i=1:(k-1)
%	m(i) = mean(xy((i-1)*lclas+1:i*lclas));
%end
%m(k) = mean(xy((k-1)*lclas+1:length(xy)));

m = mean(reshape(xy(1:k*lclas),lclas,k));
m(k) = mean(xy((k-1)*lclas+1:end));
