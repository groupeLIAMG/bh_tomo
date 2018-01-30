function [G, mk] = preFFTMA2D(xmin,dx,xmax,ymin,dy,ymax,covar)
% [G, mk] = preFFTMA2D(xmin,dx,xmax,ymin,dy,ymax,covar)
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
	
	
new = 0;

small=1e-6;
Nx=2*(1+length(xmin:dx:xmax+small));
Ny=2*(1+length(ymin:dy:ymax+small));

if new
	fact = 1;
	if covar.model(1)==2 || covar.model(1)==3
    long = xmax-xmin;
    temp = 3*covar.model(3)/long;
    fact = round(temp);
	elseif covar.model(1)==4
    long = xmax-xmin;
    temp = 2/sqrt(3)*covar.model(3)/long;
    fact = round(temp);
	end
	if fact/2 ~= round(fact/2)
    fact = fact+1;
	end
	Nx = Nx*(fact+4);
	Ny = Ny*(fact+4);
	
	Nx2 = Nx/2;
	Ny2 = Ny/2;
	
	%Erwan
	x = dx*(1:Nx2);
	x = [0 x fliplr(-x(1:end-1))]';
	y = dy*(1:Ny2);
	y = [0 y fliplr(-y(1:end-1))]';
	
	x = kron(ones(Ny,1), x);
	y = kron(y, ones(Nx,1));
	
	d = Covardm([y x],[0 0],covar.model,covar.c);
	K = reshape(d,Nx,Ny)';
	
	%%%%%%%%%%%%%%%%%%%%%%%%%%
	%On s'assure que la covariance tombe bien à zéro
	
	
	mk=0;
	if covar.model(1)>4 && min(K(:,1))>1e-6
		Ny=2*Ny;
		Nx=2*Nx;
		mk=1;
	end
	if covar.model(1)>4 && min(K(1,:))>1e-6
		Nx=2*Nx;
		Ny=2*Ny;
		mk=1;
	end  
	
	if mk==1;
		Nx2 = Nx/2;
		Ny2 = Ny/2;
		
		x = dx*(1:Nx2);
		x = [0 x fliplr(-x(1:end-1))]';
		y = dy*(1:Ny2);
		y = [0 y fliplr(-y(1:end-1))]';
		
		x = kron(ones(Ny,1), x);
		y = kron(y, ones(Nx,1));
		
		d = Covardm([y x],[0 0],covar.model,covar.c);
		
		K = reshape(d,Nx,Ny)';
	end
	
else
	
	Nx2 = Nx/2;
	Ny2 = Ny/2;
	
	x = dx*(0:Nx2-1);
	x = [x fliplr(-x)]';
	y = dy*(0:Ny2-1);
	y = [y fliplr(-y)]';
	
	x = kron(ones(Ny,1), x);
	y = kron(y, ones(Nx,1));
	d = covardm([y x],[0 0],covar.model,covar.c);
	
	K = reshape(d,Nx,Ny)';
	%imagesc(K); axis equal
	%%%%%%%%%%%%%%%%%%%%%%%%%%
	%On s'assure que la covariance tombe bien à zéro
	mk=0;
	if min(K(:,1))>small
		Ny=5*Ny;
		mk=1;
	end
	if min(K(1,:))>small
		Nx=5*Nx;
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
		d = covardm([y x],[0 0],covar.model,covar.c);
		
		K = reshape(d,Nx,Ny)';
	end
end
	
	
% Calcul de G

G=fft2(K).^0.5;
