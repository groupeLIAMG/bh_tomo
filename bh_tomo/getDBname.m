function name = getDBname(db_file)
% name = getDBname(db_file)

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


if isempty( db_file )
	name = '';
	return
end
tmp = findstr('/',db_file);
if ~isempty(tmp)
	name = db_file((tmp(length(tmp))+1):length(db_file));
else
	tmp = findstr('\',db_file); %windoze
	if ~isempty(tmp)
		name = db_file((tmp(length(tmp))+1):length(db_file));
	else
		name = dbfile;
	end
end
