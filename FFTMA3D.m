function zs=FFTMA3D(xmin,dx,xmax,ymin,dy,ymax,zmin,dz,zmax,G,mk)
% function zs=FFTMA3D(xmin,dx,xmax,ymin,dy,ymax,zmin,dz,zmax,G,mk)
%
%
				
% Copyright (C) 2006 Erwan Gloaguen, Bernard Giroux
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

small=1e-6;

[Nx,Ny,Nz] = size(G);

U=fftn(randn(size(G)));

GU=G.*U;

% Transformation de Fourier inverse donnant g*u et z
Z=real(ifftn(GU));
%reconstruction de la simulation sur la taille du champ
if mk==0
  zs=Z((Nx+2)/2+1:end, (Ny+2)/2+1:end, (Nz+2)/2+1:end);
elseif mk==1
  zs=Z((Nx+2)/2+1:(Nx+2)/2+length(xmin:dx:xmax+small),...
			 (Ny+2)/2+1:(Ny+2)/2+length(ymin:dy:ymax+small),...
			 (Nz+2)/2+1:(Nz+2)/2+length(zmin:dz:zmax+small));
end
