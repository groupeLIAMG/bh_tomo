function d = lisSEGY(fichier)
% function d = lisSEGY(fichier)


% Copyright (C) 2010 Bernard Giroux
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


s = read_segy(fichier, [], 'scalco,sx,sy,gx,gy,selev,gelev');


d.cdate = '';
d.ntrace = size( s.data, 2 );
d.nptsptrc = double(s.bh.hns);
d.ntzpt = 0;
d.nttwin = 0;
d.rstpos = 0;
d.rfpos = 0;
d.rstepsz = 0;
d.cunits = 'm';
if s.bh.mfeet == 2
    d.cunits = 'ft';
end
d.rnomfreq = 0;
d.rantsep = 0;
d.npulsev = 0;
d.nnstack = 0;
d.csurvmod = '';
d.chstrngs = '';
d.nchstrngs = 0;
d.timec = 0.001 * double(s.bh.hdt);
d.rdata = s.data;
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
d.antennas = 'SEG Y';
d.synthetique = false;
d.comment='true positions';
d.TxOffset = 0;
d.RxOffset = 0;
d.tunits = 'ms';

d.timestp = d.timec*(1:d.nptsptrc);

ind = s.th.scalco>0;
if any(ind)
    d.Tx_x(ind) = double(s.th.scalco(ind)) .* double(s.th.sx(ind));
    d.Tx_y(ind) = double(s.th.scalco(ind)) .* double(s.th.sy(ind));
    d.Tx_z(ind) = double(s.th.scalco(ind)) .* double(s.th.selev(ind));
    d.Rx_x(ind) = double(s.th.scalco(ind)) .* double(s.th.gx(ind));
    d.Rx_y(ind) = double(s.th.scalco(ind)) .* double(s.th.gy(ind));
    d.Rx_z(ind) = double(s.th.scalco(ind)) .* double(s.th.gelev(ind));
end

ind = ~ind;
if any(ind)
    d.Tx_x(ind) = double(s.th.sx(ind)) ./ abs( double(s.th.scalco(ind)) );
    d.Tx_y(ind) = double(s.th.sy(ind)) ./ abs( double(s.th.scalco(ind)) );
    d.Tx_z(ind) = double(s.th.selev(ind)) ./ abs( double(s.th.scalco(ind)) );
    d.Rx_x(ind) = double(s.th.gx(ind)) ./ abs( double(s.th.scalco(ind)) );
    d.Rx_y(ind) = double(s.th.gy(ind)) ./ abs( double(s.th.scalco(ind)) );
    d.Rx_z(ind) = double(s.th.gelev(ind)) ./ abs( double(s.th.scalco(ind)) );
end
