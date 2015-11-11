function w = warcone(t, f, t0, m, n)
% function w = warcone(t, f, t0, m, n)
%
% voir papier
%
%@ARTICLE{arcone91,
%  author = {Steven A. Arcone},
%  title = {Dielectric constant and layer-thickness interpretation
%  ofhelicopter-borne short-pulse radar waveforms reflected from wet
%  and dryriver-ice sheets},
%  journal = {IEEE Transactions on Geoscience and Remote Sensing},
%  year = {1991},
%  volume = {29},
%  pages = {768--777},
%  number = {5},
%  doi = {10.1109/36.83992}
%}

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

omega = 2*pi*f;
w = sin(omega*(t-t0)) .* (sin(omega/(2*n)*(t-t0))).^m;
ind = t<t0 | t>(n/f);
w(ind) = 0;
