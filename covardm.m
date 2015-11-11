function k=covardm(x,x0,model,c)
% [k]=covardm(x,x0,model,c)
%
% Fonction pour calculer les covariances avec des modèles
% spécifiées comme avec cokri 
% La fonction calcule pour toutes les positions de x (n1) avec toutes les
% positions de x0 (n2)  K est donc n1xn2 
% auteur D. Marcotte avril 2002

% Copyright (C) 2005 Denis Marcotte
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

% here we define the equations for the various covariograms. Any new model
% can be added here.
k=[];
Gam=['h==0                                     '; %nugget
    'exp(-h)                                  '; %exponential
    'exp(-(h).^2)                             '; %gaussian
    '1-(1.5*min(h,1)/1-.5*(min(h,1)/1).^3)    '; %spherical
    '1-h                                      '; %linear
    '1-3*min(h,1).^2+2*min(h,1).^3            '; %modele cubique
    '(h.^2).*log(max(h,eps))                  '; %spline plaque mince
    '(h.^2+1).^(-0.5)                         '; %modèle gravimétrique (Cauchy avec b=0.5)
    '(h.^2+1).^(-1.5)                         '; %modele magnétique (Cauchy avec b=1.5) 
    'sin(max(eps,h*2*pi))./max(eps,h*2*pi)    '; %effet de trou sinusoidal
    'cos(h*2*pi)                              '; %effet de trou cosinusoidal
    '1-(1.5*min(h,1)/1-.5*(min(h,1)/1).^3)+1-h'];%spherique+lineaire


% definition of some constants


[n1,d]=size(x); % d dimension de l'espace
[n2,d]=size(x0);
[rp,p]=size(c);
r=rp/p;  % nombre de structures
cx=[x(:,1:d);x0];
nm=size(model,2);

% ne pas permettre des portées de 0 en input 
model(:,2)=max(model(:,2),100*eps);

if nm>2
   model(:,3:1+d)=max(model(:,3:1+d),100*eps);
end

% calculer les covariances

k=zeros(n1*p,n2*p);
for i=1:r,

  % calculation of matrix of reduced rotated distances H

  [t1]=trans(x(:,1:d),model,i);
  [t2]=trans(x0,model,i);
  h=0;
  for id=1:d
     h=h+(t1(:,id)*ones(1,n2)-ones(n1,1)*t2(:,id)').^2;
  end
  h=sqrt(h);
  ji=(i-1)*p+1; js=i*p ;

  % evaluation of the current basic structure
  g=eval(Gam(model(i,1),:));
  k=k+kron(g,c(ji:js,:));
end
