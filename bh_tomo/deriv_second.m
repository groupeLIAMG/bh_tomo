function [D]=deriv_second(m,n)
% [D]=deriv_second(m,n)
  
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


  
  
%Calcul de la matrice des dérivées secondes;
mn=m*n;
if mn<=1e6
  %disp('tt')
  
  x=grille2d(1,1,n,1,1,m); 
  x=fliplr(x);
  D=-4*speye(m*n,m*n);   
  for i=1:length(x)
	ind=find(abs(x(i,1)-x(:,1))+abs(x(i,2)-x(:,2))==1);
	D(i,ind)=1;
	% ind2=find(x(i,1)==x(:,1)&x(i,2)==x(:,2));
	%	 D(i,ind2)=-4;
  end 
elseif mn>1e6    
  if m>n %on parcours le champ en commencant par
		 %la valeur en haut à droite et on descend les z
		 %
		 %d=zeros(min(m,n),max(m,n)); 
	d=zeros(n,m);
	d(1,1)=-4;
	d(2,1)=1;
	d(1,2)=1;
	%d(2,2)=1;
    l=1;ll=1;
    k=1;
    for jj=1:m%max(m,n)
	  for ii=1:n%min(m,n)
		for j=1:m%max(m,n)
		  for i=1:n%min(m,n)
			D(ll,l)=d(abs(ii-i)+1,abs(jj-j)+1);
			k=k+1;
			l=l+1;
		  end
		end
		ll=ll+1;
		l=1;
	  end
	end
	
  else 	%d=zeros(min(m,n),max(m,n)); 
	d=zeros(n,m);
	d(1,1)=-4;
	d(2,1)=1;
	d(1,2)=1;
	%d(2,2)=1;
    l=1;ll=1;
    k=1;
    for jj=1:m%max(m,n)
	  for ii=1:n%min(m,n)
		for j=1:m%max(m,n)
		  for i=1:n%min(m,n)
			D(ll,l)=d(abs(ii-i)+1,abs(jj-j)+1);
			k=k+1;
			l=l+1;
		  end
		end
		ll=ll+1;
		l=1;
	  end
	end 
  end 
end %m*n>1e6