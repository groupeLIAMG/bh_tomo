function d = lisSU( fichier )
% function d = lisSU( fichier )
%
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

[data, tr_header, header] = ReadSu(fichier);

if isempty(data)
	disp(['Fichier SU vide ou absent: ',fichier])
	d=-1;
	return
end

d.cdate = '';
d.ntrace = length(tr_header);
d.nptsptrc = header.ns;
d.ntzpt = 0;
d.nttwin = 0;
d.rstpos = 0;
d.rfpos = 0;
d.rstepsz = 0;
d.cunits = 'm';
d.rnomfreq = 0;
d.rantsep = 0;
d.npulsev = 0;
d.nnstack = 0;
d.csurvmod = '';
d.chstrngs = '';
d.nchstrngs = 0;
d.timec = 1e-3*header.dt;
d.rdata = data;
d.ltime = 0;
d.timestp = 0;
d.xloc = 0;
d.Tx_x = zeros(1,d.ntrace);
d.Tx_y = zeros(1,d.ntrace);
d.Tx_z = zeros(1,d.ntrace);
d.Rx_x = zeros(1,d.ntrace);
d.Rx_y = zeros(1,d.ntrace);
d.Rx_z = zeros(1,d.ntrace);
d.sigpos = 0;
d.rsigpos = 0;
d.antennas = 'Seismic Un*x';
d.synthetique = false;
d.comment='true positions';
d.TxOffset = 0;
d.RxOffset = 0;
d.tunits = 'ms';

d.timestp = d.timec*(1:d.nptsptrc);

for n=1:d.ntrace
	d.Tx_x(n) = tr_header(n).SourceX;
	d.Tx_y(n) = tr_header(n).SourceY;
	d.Tx_z(n) = tr_header(n).SourceSurfaceElevation;
	d.Rx_x(n) = tr_header(n).GroupX;
	d.Rx_y(n) = tr_header(n).GroupY;
	d.Rx_z(n) = tr_header(n).ReceiverGroupElevation;
end
