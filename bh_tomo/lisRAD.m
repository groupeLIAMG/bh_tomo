function d = lisRAD( basename )
% function d = lisRAD( basename )

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


fid=fopen(basename,'rt');
if fid < 2
  fichier = [basename,'.rad'];
  fid=fopen(fichier,'rt');
  if fid < 2
	fichier = [basename,'.RAD'];
	fid=fopen(fichier,'rt');
	if fid < 2
	  disp(['Impossible d''ouvrir le fichier RAD: ',basename])
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
d.cunits = 'm';
d.rnomfreq = 0;
d.rantsep = 0;
d.npulsev = 300;
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
d.tunits = 'ns';

while 1
  tline = fgetl(fid);
  if ~ischar(tline), break, end
  if ~isempty(findstr(tline,'SAMPLES:'))
		d.nptsptrc = sscanf(tline,'%*8c%d',1);
  elseif ~isempty(findstr(tline,'FREQUENCY:'))
		d.timec = sscanf(tline,'%*10c%f',1);
  elseif ~isempty(findstr(tline,'FREQUENCY STEPS:45'))
  elseif ~isempty(findstr(tline,'SIGNAL POSITION:')) && tline(1) == 'S'
		d.sigpos = sscanf(tline, '%*16c%f',1);
  elseif ~isempty(findstr(tline,'RAW SIGNAL POSITION:'))
		d.rsigpos = sscanf(tline, '%*20c%d',1);
  elseif ~isempty(findstr(tline,'DISTANCE FLAG:1'))
  elseif ~isempty(findstr(tline,'TIME FLAG:0'))
  elseif ~isempty(findstr(tline,'PROGRAM FLAG:0'))
  elseif ~isempty(findstr(tline,'EXTERNAL FLAG:0'))
  elseif ~isempty(findstr(tline,'TIME INTERVAL:0.100000'))
  elseif ~isempty(findstr(tline,'DISTANCE INTERVAL:'))
		d.rstepsz = sscanf(tline, '%*18c%f',1);
  elseif ~isempty(findstr(tline,'OPERATOR:'))
		d.synthetique = ~isempty(regexp(tline,'MoRad', 'once')) || ...
				~isempty(regexp(tline,'synthetic', 'once'));
  elseif ~isempty(findstr(tline,'CUSTOMER:'))
  elseif ~isempty(findstr(tline,'SITE:'))
  elseif ~isempty(findstr(tline,'ANTENNAS:'))
		[debut,fin] = regexp(tline,'\d+');
		d.rnomfreq = str2double( tline(debut:fin) );
		d.antennas = tline(10:length(tline));
  elseif ~isempty(findstr(tline,'ANTENNA ORIENTATION:NOT VALID FIELD'))
  elseif ~isempty(findstr(tline,'ANTENNA SEPARATION'))
		d.rantsep = sscanf(tline, '%*s %*11c%f', 1);
  elseif ~isempty(findstr(tline,'COMMENT:'))
		d.comment = tline(9:length(tline));
  elseif ~isempty(findstr(tline,'TIMEWINDOW'))
		d.ltime = sscanf(tline, '%*11c%f',1);
  elseif ~isempty(findstr(tline,'STACKS'))
		d.nnstack = sscanf(tline,'%*7c%d',1);
  elseif ~isempty(findstr(tline,'STACK EXPONENT:5'))
  elseif ~isempty(findstr(tline,'STACKING TIME:0.163840'))
  elseif ~isempty(findstr(tline,'LAST TRACE'))
		d.ntrace = sscanf(tline,'%*s %*6c%d',1);
  elseif ~isempty(findstr(tline,'STOP POSITION:'))
		d.rfpos = sscanf(tline, '%*14c%f',1);
  elseif ~isempty(findstr(tline,'SYSTEM CALIBRATION:0.0000233154'))
  elseif ~isempty(findstr(tline,'START POSITION:2.0000000000'))
  elseif ~isempty(findstr(tline,'SHORT FLAG:0'))
  elseif ~isempty(findstr(tline,'INTERMEDIATE FLAG:1'))
  elseif ~isempty(findstr(tline,'LONG FLAG:0'))
  elseif ~isempty(findstr(tline,'PREPROCESSING:0'))
  elseif ~isempty(findstr(tline,'HIGH:0'))
  elseif ~isempty(findstr(tline,'LOW:0'))
  elseif ~isempty(findstr(tline,'FIXED INCREMENT:0.3000000000'))
  elseif ~isempty(findstr(tline,'FIXED MOVES UP:0'))
  elseif ~isempty(findstr(tline,'FIXED MOVES DOWN:1'))
  elseif ~isempty(findstr(tline,'FIXED POSITION:3.9000000000'))
  elseif ~isempty(findstr(tline,'WHEEL CALIBRATION:0.0000000000'))
  elseif ~isempty(findstr(tline,'POSITIVE DIRECTION:1'))
  end
end
d.timec = 1000/d.timec;
d.timestp = d.timec*(0:d.nptsptrc-1);
if d.synthetique==0
	d.antennas = strcat(d.antennas, ' - Ramac');
end

fclose(fid);
