function d = lisTLF( basename )
% d = lisTLF( basename )

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

[fid,mesg]=fopen(basename,'rt');
if fid < 2
  fichier = [basename,'.tlf'];
  [fid,mesg]=fopen(fichier,'rt');
  if fid < 2
	fichier = [basename,'.TLF'];
	[fid,mesg]=fopen(fichier,'rt');
	if fid < 2
%	  disp(['Fichier TLF absent: ',basename])
	  d = -1;
	  return
	end
  end
end
tline = fgetl(fid);
ind = 1;

d.Tx_z = [];
d.Rx_z = [];
while ~feof(fid)
  tnd = fscanf(fid,'%d', 1);
  tnf = fscanf(fid,'%d', 1);
  nt = tnf-tnd+1;

  Rxd = fscanf(fid,'%f', 1);
  Rxf = fscanf(fid,'%f', 1);
  if nt == 1
	dRx = 1;
	if Rxd>Rxf, Rxd = Rxf; end
  else
	dRx = (Rxf-Rxd)/(nt-1);
  end
  Tx  = fscanf(fid,'%f', 1);
  
  if nt > 0, d.Tx_z = [d.Tx_z Tx*ones(1, nt)]; end
  d.Rx_z = [d.Rx_z Rxd:dRx:Rxf];
  fgetl(fid);
end

fclose(fid);
