function d = lisRAMAC2( basename )
% function d = lisRAMAC2( basename )

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

d = lisRAD( basename );
if ~isstruct(d)
  return
end

d.rdata = lisRD3( basename, d.nptsptrc, d.ntrace );
if d.rdata==-1
  return
end

pos = lisTLF( basename, d.ntrace );
if isstruct(pos)
  d.Tx_z = pos.Tx_z;
  d.Rx_z = pos.Rx_z;
else
    d = [];
    return
end

% Offset du pt milieu des antennes forage p/r au colet du trou
d.TxOffset = 0;
d.RxOffset = 0;
if d.synthetique == 0
	if d.rnomfreq == 100
		d.TxOffset = 0.665;
		d.RxOffset = 0.665;
	elseif d.rnomfreq == 250
		d.TxOffset = 0.325;
		d.RxOffset = 0.365;
	end
end

d.Tx_z = d.Tx_z(1:d.ntrace);
d.Rx_z = d.Rx_z(1:d.ntrace);

d.Rx_x = zeros(1,d.ntrace);
d.Rx_y = zeros(1,d.ntrace);
d.Tx_x = zeros(1,d.ntrace);
d.Tx_y = zeros(1,d.ntrace);
