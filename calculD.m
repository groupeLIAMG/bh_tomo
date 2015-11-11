function D = calculD(x,modele,c,m,n)
% function D = calculD(x,modele,c,m,n)
%
% modele de covariance d'une grille 2D
%
% x : vecteur positions des noeuds de la grille
% modele : modele de covariance
% m : n pts en x
% n : n pts en z

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


if m*n>1e9
  if m<=n
%	d=covardm(x,[0 0],modele,[0; 1]);
%	d=covardm(x,[0 0],modele(2,:),[1]);
	d=covardm(x,[0 0],modele,c);
	d=reshape(d,n,m)';
  else
%	d=covardm([0 0],x,modele,[0;1]);
%	d=covardm([0 0],x,modele(2,:),[1]);
	d=covardm([0 0],x,modele,c);
	d=reshape(d,n,m);
  end
  m1=min(m,n);
  m2=max(m,n);
  D = zeros( m1*m2 );
  l=1;
  ll=1;
  k=1;
  for jj=1:m1
	for ii=1:m2
	  for j=1:m1
		for i=1:m2
		  D(ll,l)=d(abs(jj-j)+1,abs(ii-i)+1);
		  k=k+1;
		  l=l+1;
		end
	  end
	  ll=ll+1;
	  l=1;
	end
  end
else
	%  D=covardm(x,x,modele,[0;1]);
	%  D=covardm(x,x,modele(2,:),[1]);
	D=covardm(x,x,modele,c);
end
