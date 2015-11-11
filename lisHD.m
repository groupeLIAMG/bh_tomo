function d = lisHD( basename )
% function d = lisHD( basename )

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
  fichier = [basename,'.hd'];
  [fid,mesg]=fopen(fichier,'rt');
  if fid < 2
	fichier = [basename,'.HD'];
	[fid,mesg]=fopen(fichier,'rt');
	if fid < 2
	  disp(['Impossible d''ouvrir le fichier HD: ',basename])
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
d.Tx_z = 0;
d.Rx_x = 0;
d.Rx_z = 0;
d.sigpos = 0;
d.rsigpos = 0;
d.antennas = '';
d.synthetique = 0;

tline = fgetl(fid);
d.comment = fgetl(fid);
tline = fgetl(fid);
if strcmp(tline(1:3),'VRP'),
  vrp = 1;
  ic = findstr(tline,'=');
  count = length(tline);
  Tx_x = sscanf(tline(ic(1)+1:count),'%f');
  Rx_x = sscanf(tline(ic(2)+1:count),'%f');
  Rx_z = sscanf(tline(ic(3)+1:count),'%f');
elseif strcmp(tline(1:3),'ZOP'),
elseif strcmp(tline(1:3),'MOG'),
  mog = 1;
  ic = findstr(tline,'=');
  count = length(tline);
  Tx_z = sscanf(tline(ic(1)+1:count),'%f');
  Rx_z = sscanf(tline(ic(2)+1:count),'%f');
end

d.cdate = fgetl(fid);

while 1
  tline = fgetl(fid);
  if ~ischar(tline), break, end
  if ~isempty(findstr(tline,'NUMBER OF TRACES'))
	  d.ntrace = sscanf(tline,'%*20c %d',1);
  elseif ~isempty(findstr(tline,'NUMBER OF PTS/TRC'))
	  d.nptsptrc = sscanf(tline,'%*20c %d',1);
  elseif ~isempty(findstr(tline,'TIMEZERO AT POINT'))
	  d.ntzpt = sscanf(tline,'%*20c %f',1);
  elseif ~isempty(findstr(tline,'TOTAL TIME WINDOW'))
	  d.nttwin = sscanf(tline,'%*20c %f',1);
  elseif ~isempty(findstr(tline,'STARTING POSITION'))
	  d.rstpos = sscanf(tline,'%*20c %f',1);
  elseif ~isempty(findstr(tline,'FINAL POSITION'))
	  d.rfpos = sscanf(tline,'%*20c %f',1);
  elseif ~isempty(findstr(tline,'STEP SIZE USED'))
	  d.rstepsz = sscanf(tline,'%*20c %f',1);
  elseif ~isempty(findstr(tline,'POSITION UNITS'))
	  d.cunits = sscanf(tline,'%*20c %s',1);
  elseif ~isempty(findstr(tline,'NOMINAL FREQUENCY'))
	  d.rnomfreq = sscanf(tline,'%*20c %f',1);
  elseif ~isempty(findstr(tline,'ANTENNA SEPARATION'))
	  d.rantsep = sscanf(tline,'%*20c %f',1);
  elseif ~isempty(findstr(tline,'PULSER VOLTAGE (V)'))
	  d.npulsev = sscanf(tline,'%*20c %f',1);
  elseif ~isempty(findstr(tline,'NUMBER OF STACKS'))
	  d.nnstack = sscanf(tline,'%*20c %d',1);
  elseif ~isempty(findstr(tline,'SURVEY MODE'))
	  d.csurvmod = tline(21:end);
  end
end

fclose(fid);

d.ltime = d.nttwin;
d.timec = d.nttwin/d.nptsptrc;
d.timestp = d.timec*(1:d.nptsptrc);


d.Tx_x = zeros(1,d.ntrace);
d.Tx_y = zeros(1,d.ntrace);
d.Tx_z = Tx_z*ones(1,d.ntrace);
d.Rx_x = d.rantsep*ones(1,d.ntrace);
d.Rx_y = zeros(1,d.ntrace);
d.Rx_z = Rx_z + d.rstepsz*(0:(d.ntrace-1));
