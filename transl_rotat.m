function m_data = transl_rotat(data, origine, az, dip)
% m_data = transl_rotat(data, origine, az, dip)
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

% translation p/r origine
m_data = data-repmat(origine,size(data,1),1);
% rotation p/r azimuth
if abs(az) > (pi/720)  % si plus grand que 1/4 degre
    rot = [cos(az) -sin(az); sin(az) cos(az)];
    for n=1:size(m_data,1)
        m_data(n,1:2) = m_data(n,1:2)*rot';
    end
    % data(n,1:2)*rot' est egal a (rot*data(n,1:2)')'
end

% rotation p/r pendage
if abs(dip) > (pi/720)
    rot = [cos(dip) -sin(dip); sin(dip) cos(dip)];
    for n=1:size(m_data,1)
        m_data(n,2:3) = m_data(n,2:3)*rot';
    end
    % data(n,2:3)*rot' est egal a (rot*data(n,2:3)')'
end
