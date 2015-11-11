function Ow = topp(k,varargin)
% function Ow = topp(k, dir)
%
% Calcul la teneur en eau Ow a partir de la constante dielectrique k
%
%  dir (optionel) : si 1, retourne Ow
%                   si -1, retourne k (l'input doit alors etre Ow)
%
%

% Copyright (C) 2008 Bernard Giroux
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

if nargin == 2 && varargin{1}==-1
  Ow = 3.03 + 9.3*k + 146*k.^2 - 76.7*k.^3;
else
  Ow = -5.3e-2 + 2.92e-2*k - 5.5e-4*k.^2 + 4.3e-6*k.^3;
end
