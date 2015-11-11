function [x,residu_moy,residu_ecar]=lsqr_erwan3d(L,dt,iter,epsilon,gradmin,correlmin,scont,x,l_moy,indi)
% [x,residu_moy,residu_ecar]=lsqr_erwan3d(L,dt,iter,epsilon,gradmin,correlmin)
% Resolution tomographique par methode LSQR
%
%        L        = matrice dite de tomographie 
%                    nombre de colonne = nombre d'inconnue
%                    nombre de ligne = nombre de rais
%        dt        = parametres a inverser observes
%        iter     = nombre d'iterations
%        epsilon  = tolerance sur la fonction objective
%        gradmin  = gradient minimum de la fonction objective
%        correlmin= correlation d'une iteration a l'autre
%
%        x           = resultat de l'inversion
%        residu_moy  = moyenne des residus de chaque iterartion
%        residu_ecar = ecart-type des residus de chaque iteration


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

L=sparse(L);
if ~isempty( scont )
	%Calcul des indices des donnees contraignantes
	indc = zeros(length(scont(:,1)), 1);
	for i=1:length(scont(:,1))
		ind11=findnear(scont(i,1),x(:,1));
		ind22=findnear(scont(i,2),x(ind11,2));
		ind33=findnear(scont(i,3),x(ind11(ind22),3));
		indc(i)=min(ind11(ind22(ind33)));
	end
	clear ind11 ind22 ind33
	
	npx=length(scont(:,1));
	K=zeros(npx,length(L(1,:)));
	for ndc=1:npx
		K(ndc,indc(ndc))=1;
	end
	L=[L;K];
	
	%	 W=diag(1./sqrt(sum(L,1)));
	%    L=L*W;
	dt=[dt;scont(:,4)-l_moy];
end


xx=zeros(length(L(1,:)),iter);
% (1) initialisation

beta=norm(dt);
u=dt/beta;

tmp=L'*u;
alpha=norm(tmp);
v=tmp/alpha;

w=v;

x=zeros(size(L,2),1);

phi_n=beta;

ro_n=alpha;

residu_moy = zeros(iter,1);
residu_ecar = zeros(iter,1);

% (2) Boucle
h = waitbar(0,'Iterations LSQR ...');
for n=1:iter
    waitbar(n/iter+1)
    
    % (3) continue bidiagonalisation
    
    % (3a)
    tmp=L*v-alpha*u;
    beta=norm(tmp);
    u=tmp/beta;
    
    % (3b)
    tmp=L'*u-beta*v;
    alpha=norm(tmp);
    v=tmp/alpha;
    
    % (4) construct and apply nexte orthogonal transformation
    
    % (4a)
    ro=sqrt(ro_n*ro_n+beta*beta);
    
    % (4b)
    c=ro_n/ro;
    
    % (4c)
    s=beta/ro;
    
    % (4d)
    theta=s*alpha;
    
    % (4e)
    ro_n=-1*c*alpha;
    
    % (4f)
    phi=c*phi_n;
    
    % (4g)
    phi_n=s*phi_n;
    
    % (5) Mise a jour du nouveau vecteur solution
    
    % (5a)
    x=x+(phi/ro)*w;
    
    % (5b)
    w=v-(theta/ro)*w;
    
    % Controle de la convergence
    
    res=L*x-dt;
    residu_moy(n)=sqrt((res'*res)./length(res));
    residu_ecar(n)=std(res);
    %Construction de la matrice de l'ensemble des resultats
    xx(:,n)=x;
    
    %Sortie sur differents criteres de convergences
    if n>=2
        correlminm=corrcoef(xx(:,n),xx(:,n-1));
        correl=correlminm(2,1);
        gradfo=abs(residu_moy(n)-residu_moy(n-1));
        if residu_moy(n)<epsilon
            disp('Sortie sur epsilon')
            break
        elseif gradfo<gradmin
            disp('Sortie sur gradfo')
            break
        elseif correl>correlmin
            disp('Sortie sur correl')
            break
        end
    end
end
close(h)

residu_moy = residu_moy(1:n);
residu_ecar = residu_ecar(1:n);
