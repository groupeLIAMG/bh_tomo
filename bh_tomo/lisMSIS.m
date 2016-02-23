function d = lisMSIS( basename )
% function d = lisMSIS( basename )
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

fid=fopen(basename,'rt');
if fid < 2
  fichier = [basename,'.dat'];
  fid=fopen(fichier,'rt');
  if fid < 2
	fichier = [basename,'.DAT'];
	[fid,mesg]=fopen(fichier,'rt');
	if fid < 2
	  disp(['Impossible d''ouvrir le fichier MSIS: ',basename,' - ',mesg])
	  d=-1;
	  return
	end
  end
end

d.cdate = '';
d.ntrace = 0;
d.nptsptrc = 0;
d.ntzpt = 0;
d.nttwin = 0;
d.rstpos = 0;
d.rfpos = 0;
d.rstepsz = 0;
d.cunits = '';
d.rnomfreq = 0;
d.rantsep = 0;
d.npulsev = 0;
d.nnstack = 0;
d.csurvmod = '';
d.chstrngs = '';
d.nchstrngs = 0;
d.timec = 0;
d.rdata = 0;
d.ltime = 0;
d.timestp = 0;
d.xloc = 0;
d.Tx_x = 0;
d.Tx_y = 0;
d.Tx_z = 0;
d.Rx_x = 0;
d.Rx_y = 0;
d.Rx_z = 0;
d.sigpos = 0;
d.rsigpos = 0;
d.antennas = 'Microsis';
d.synthetique = false;
d.comment='';
d.TxOffset = 0;
d.RxOffset = 0;
d.tunits = 'ms';
d.cunits = 'm';

while 1
  tline = fgetl(fid);
  if ~ischar(tline), break, end
  if ~isempty(findstr(tline,'Number of channel     :	 '))
	[debut,fin] = regexp(tline,'\d+');
   	d.ntrace = str2double( tline(debut:fin) );
    %d.ntrace = sscanf(tline,'%*25c%d',1);
  elseif ~isempty(findstr(tline,'Sampling rate         :	 '))
    [debut,fin] = regexp(tline,'\d+');
   	d.timec = str2double( tline(debut:fin) );  % MHz
  elseif ~isempty(findstr(tline,'Number of saved pts  :	 '))  
      [debut,fin] = regexp(tline,'\d+');
      d.nptsptrc = str2double( tline(debut:fin) );
      %d.nptsptrc = sscanf(tline,'%*24c%d',1);
  elseif ~isempty(findstr(tline,'***************************'))
    position=ftell(fid);
  end
end
fseek (fid, position, 'bof');
d.rdata = fscanf(fid, '%f',[d.ntrace inf]);
fclose(fid);

d.timec = 1e-3/d.timec;
d.timestp = d.timec*(1:d.nptsptrc);

pos = lisTLF( basename );
if isstruct(pos)
  d.Tx_z = pos.Tx_z;
  d.Rx_z = pos.Rx_z;
else
    d = [];
    return
end

d.Tx_z = d.Tx_z(1:d.ntrace);
d.Rx_z = d.Rx_z(1:d.ntrace);

d.Rx_x = zeros(1,d.ntrace);
d.Rx_y = zeros(1,d.ntrace);
d.Tx_x = zeros(1,d.ntrace);
d.Tx_y = zeros(1,d.ntrace);


d.rdata = d.rdata';
