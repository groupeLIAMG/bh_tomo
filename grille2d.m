function out=grille2d(xmin,pasx,nx,ymin,pasy,ny)
% function out=grille2d(xmin,pasx,nx,ymin,pasy,ny)
% grille2d : define a 2d grid
%

% Copyright (C) 2005 Erwan Gloaguen
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

out=[kron(ones(ny,1),(1:nx)'*pasx), kron((1:ny)',ones(nx,1)*pasy)];
out(:,1)=xmin+out(:,1)-pasx;
out(:,2)=ymin+out(:,2)-pasy;
