function [kss]=kss_it(x,dx,dz,covar,m,n,indi);
%function kss=kss_it(x,dx,dz,covar,m,n,indi);
%  Fonction qui calcul une matrice de covariance de bloc
%  sur une grille 2d régulière 

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

h=1;

ind=1:m*n;

lindi=length(ind);
d=covardm(min(x),x,covar.modele,covar.c);

x(:,1)=round((x(:,1)-x(1,1))./dz+1);
x(:,2)=round((x(:,2)-x(1,2))./dx+1);

if lindi<1e4
	kss=zeros(n*m,m*n);%
else
%	kss = sparse([],[],[],lindi,lindi,lindi);
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if lindi<1e1
	k=2;
	for jj=1:m
		for ii=1:n
			kss(k-1,k:lindi)=d(int32((abs(jj-x(k:lindi,2)))*n+(abs(ii-(x(k:lindi,1)))+1)));
			k=k+1;
		end
	end
	
	kss=kss+kss';%+speye(lindi)*(ca(1)+ca(2));
	
	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%5
elseif lindi<1e3
	k=2;
	h = waitbar(0,'Calcul de la matrice de covariance des parametres (2)...');
	for i=1:length(x(:,1))
		waitbar(i/lindi,h)
		%%	kss(k-1,k:lindi)=d(abs(x(i,2)-(x(k:lindi,2))*n+abs(x(i,1)-(x(k:lindi,1)))+1));
		kss(k-1,k:lindi)=d((abs(x(i,2)-(x(k:lindi,2)))*n+abs(x(i,1)-(x(k:lindi,1)))+1));
		k=k+1;
	end
	kss=kss+kss';%+speye(lindi)*(ca(1)+ca(2));
else%%%%%%%%%%%%%%%nouvelle m�thode
	k=2;
	h = waitbar(0,'Calcul de la matrice de covariance des parametres (3)...');
	for i=1:length(x(:,1))
		waitbar(i/lindi,h)
		kss(k-1,1:lindi)=d((abs(x(i,2)-(x(1:lindi,2)))*n+abs(x(i,1)-(x(1:lindi,1)))+1));
		k=k+1;
	end
	
end

%%%%%
%on ferme le waitbar
close(h)
% si indi contient des zeros
if nnz(indi)~=numel(indi)
	kss=kss(indi,indi);
end
