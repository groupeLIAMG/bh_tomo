function e_eff = self_similar( e_m, e_f, phi, psi )
%SELF_SIMILAR Compute the effective property of a two-phase mixture
%
%   e_eff = self_similar( e_m, e_f, phi, psi )
%             e_m : property of the host medium
%             e_f : property of the inclusion
%             phi : volumetric proportion of e_f
%             psi : depolarization factor
%
%   Reference:
%
%@Article{sen81a,
%  Title   = {A self-similar model for sedimentary rocks with application to the dielectric constant of fused glass beads},
%  Author  = {P. N. Sen and C. Scala and M. H. Cohen},
%  Journal = {Geophysics},
%  Year    = {1981},
%  Number  = {5},
%  Pages   = {781-795},
%  Volume  = {46},
%  Doi     = {10.1190/1.1441215}
%}
%
%  B. Giroux
%  INRS-ETE
%  2014-06-10

% Copyright (C) 2014 Bernard Giroux
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

% starting value is the average
e0 = (e_m+e_f)/2;

if exist('optimoptions')
	options = optimoptions('fsolve','Display','off');
else
	options = optimset('Display','off');
end

% objective function
f = @(x)func(x, e_m, e_f, phi, psi);

e_eff = fsolve(f, e0, options);
end

function y = func(x, e_m, e_f, phi, psi)
y = (x-e_m)/(e_f - e_m)*(e_f/x)^psi - phi;
end
