function c = updateObj(c)
%
%

% Copyright (C) 2008 Bernard Giroux
%
% This program is free software; you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation; either version 2 of the License, or
% (at your option) any later version.
%
% This program is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details.
% 
% You should have received a copy of the GNU General Public License
% along with this program; if not, write to the Free Software
% Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA 02110-1301 USA
% 

if strcmpi(c.type, 'ramac 100')
    
	disp('*** Correction des antennes: mise à jour pour Ramac 100 ***')
	c.length = 2*0.48;
	c.air_data = load('usgs_air_ant_slown.dat');
    c.water_data = load('usgs_water_ant_slown.dat');
    c.k_formation = 4:2:30;
    c.hole_diam = 0.1:0.02:0.2;

elseif strcmpi(c.type, 'ramac 250')

	disp('*** Correction des antennes: mise à jour pour Ramac 250 ***')
	c.length = 2*0.5;
	c.air_data = load('ramac_250_air_ant_slown.dat');
	c.water_data = load('ramac_250_water_ant_slown.dat');
	c.k_formation = 2:2:36;
	c.hole_diam = 0.05:0.01:0.16;
   
elseif strcmpi(c.type, 'usgs')
	disp('*** Correction des antennes: mise à jour pour USGS ***')
    c.length = 2*0.48;
    c.air_data = load('usgs_air_ant_slown.dat');
    c.water_data = load('usgs_water_ant_slown.dat');
    c.k_formation = 4:2:30;
    c.hole_diam = 0.1:0.02:0.2;
else
    ME = MException('AntennaLengthCorrection:InvalidInput','Antenna type not recognized');
    throw(ME);
end