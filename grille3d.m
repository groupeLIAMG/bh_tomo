function out=grille3d(x0,dx,nx,y0,dy,ny,z0,dz,nz)
% function out=grille3d(x0,dx,nx,y0,dy,ny,z0,dz,nz)
% grille3d : define a 3d grid
%

% Copyright (C) 2013 Bernard Giroux
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

out=[kron(ones(ny*nz,1),(1:nx)'*dx) ...
		 kron(kron(ones(nz,1),(1:ny)'*dy),ones(nx,1)) ...
		 kron((1:nz)'*dz,ones(nx*ny,1))];

out(:,1)=x0+out(:,1)-dx;
out(:,2)=y0+out(:,2)-dy;
out(:,3)=z0+out(:,3)-dz;
