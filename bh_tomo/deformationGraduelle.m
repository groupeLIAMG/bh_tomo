function T = deformationGraduelle(S, L, c0, tobs, varargin)
% function T = deformationGraduelle(S, L, c0, tobs, t_handle)
%
%
  
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


if nargin>=5
    set(varargin{1}, 'String','Geostatistical Inversion - Gradual Deformation ...')
end
n=log2(size(S,2));

if isempty( c0 )
	w = 1;
else
	w = max(c0)./c0;
end

while n>0
    T = zeros(size(S,1), size(S,2)/2);
    ind=1;
    for k=1:2:2^n
		banana= @(x)sum(w.*(L*(cos(x)*S(:,k)+sin(x)*S(:,k+1))-tobs).^2);
        x = fminbnd(banana,0,2*pi);
        T(:,ind)=(cos(x)*S(:,k)+sin(x)*S(:,k+1));
        ind=ind+1;
    end
    S = T;
    n=n-1;
end
