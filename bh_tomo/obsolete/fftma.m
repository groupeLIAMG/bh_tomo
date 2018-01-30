function [zs]=fftma(xmin,dx,xmax,ymin,dy,ymax,covar)
% [zs]=fftma(xmin,dx,xmax,ymin,dy,ymax,covar)
  
% Copyright (C) 2005 Erwan Gloaguen, Bernard Giroux
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

Nx=2*(1+length(xmin:dx:xmax+small));
Ny=2*(1+length(ymin:dy:ymax+small));

Nx2 = Nx/2;
Ny2 = Ny/2;

x = dx*(0:Nx2-1);
x = [x fliplr(-x)]';
y = dy*(0:Ny2-1);
y = [y fliplr(-y)]';

x = kron(ones(Ny,1), x);
y = kron(y, ones(Nx,1));
d = covardm([y x],[0 0],covar.modele,covar.c);

K = reshape(d,Nx,Ny)';
%imagesc(K); axis equal
%%%%%%%%%%%%%%%%%%%%%%%%%%
%On s'assure que la covariance tombe bien à zéro
mk=0;
if min(K(:,1))>small
  Ny=2*Ny;
  mk=1;
end
if min(K(1,:))>small
  Nx=2*Nx;
  mk=1;
end  

if mk==1;
  Nx2 = Nx/2;
  Ny2 = Ny/2;

  x = dx*(0:Nx2-1);
  x = [x fliplr(-x)]';
  y = dy*(0:Ny2-1);
  y = [y fliplr(-y)]';

  x = kron(ones(Ny,1), x);
  y = kron(y, ones(Nx,1));
  d = covardm([y x],[0 0],covar.modele,covar.c);

  K = reshape(d,Nx,Ny)';
end

% Calcul de G

G=fft2(K).^0.5;

U=fft2(randn(size(K)));

GU=G.*U;

% Transformation de Fourier inverse donnant g*u et z
Z=real(ifft2(GU));
%reconstruction de la simulation sur la taille du champ
if mk==0
  zs=Z((Ny+2)/2+1:end,(Nx+2)/2+1:end);
elseif mk==1
  zs=Z((Ny+2)/2+1:(Ny+2)/2+length(ymin:dy:ymax+small),(Nx+2)/2+1:(Nx+2)/2+length(xmin:dx:xmax+small));
end
