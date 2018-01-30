function d = lisRD3( basename, npts, ntrc )
% function d = lisRD3( basename, npts, ntrc )

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

[fid,mesg]=fopen(basename,'r');
if fid < 2
  fichier = [basename,'.rd3'];
  [fid,mesg]=fopen(fichier,'r');
  if fid < 2
	fichier = [basename,'.RD3'];
	[fid,mesg]=fopen(fichier,'r');
	if fid < 2
	  disp(['Impossible d''ouvrir le fichier RD3: ',basename])
	  d=-1;
	  return
	end
  end
end

[d, c] = fread(fid, [npts ntrc], 'int16','ieee-le');

fclose(fid);
