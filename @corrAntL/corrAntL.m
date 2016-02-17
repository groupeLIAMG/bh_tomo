function c = corrAntL(arg)
% corrAntL - Constructeur d'une classe permettant d'effectuer une
% correction relative à la longueur finie des antennes radar
%
%   c = corrAntL(arg) 
%         arg est le type d'antenne: 'Ramac 100', 'Ramac 250' ou 'USGS'
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

if nargin == 0
    c.type = 'Ramac 100';
    c.length = [];
    c.air_data = [];
    c.water_data = [];
    c.k_formation = [];
    c.hole_diam = [];
    c = class(c, 'corrAntL');
    c = updateObj(c);
elseif isa(arg, 'corrAntL')
    c = arg;
elseif isa(arg, 'char')
    c.type = arg;
    c.length = [];
    c.air_data = [];
    c.water_data = [];
    c.k_formation = [];
    c.hole_diam = [];
    c = class(c, 'corrAntL');
    c = updateObj(c);
else
    ME = MException('AntennaLengthCorrection:InvalidInput','Give antenna type in char');
    throw(ME);
end
