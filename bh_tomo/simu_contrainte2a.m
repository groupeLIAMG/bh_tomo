function [S_sim,sco]=simu_contrainte2a(L,dt,covar,indc,nt,gridx,gridz,nsim,scont,mode,c0)
%%[S_sim,sco]=simu_contrainte(L,dt,covar,indc,nt,l_moy,gridx,gridz,nsim,scont,mode,c0);
%
% mode=1 si cokri seulement et 2 si cokri+simu

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


global kss ks0 Css


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%5
%Calcul de ktt theorique
ktt = L*kss*L';
if isempty(c0)==1
    ktt = ktt + covar.nugget_t*eye(nt);
else
    ktt = ktt + diag(covar.nugget_t.*c0(1:nt));
end
%Calcul de la covariance croisee
kts = L*kss;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%55
%Matrice de gauche
K=[Css kts(:,indc)';kts(:,indc) ktt];

%Regularisation
K=K+eye(size(K))*1e-6;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%5
%Cokrigeage primal
%tic
Lambda = (K\[ks0;kts])';
sco = Lambda*[scont;dt];
%toc

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%5
%Cokrigeage dual
%tic
%Gamma = (K\[scont;dt])';
%sco = (Gamma*[ks0;kts])';
%toc

S_sim=0;

if mode==2
    S_sim=zeros(length(sco),nsim);
    
    [Gs, mks] = preFFTMA2D(0,gridx(2)-gridx(1),max(gridx)-min(gridx),...
													 0,gridz(2)-gridz(1),max(gridz)-min(gridz), ...
													 covar);
		tmp.modele = covar.modele_xi;
		tmp.c = covar.c_xi;
    [Gxi, mkxi] = preFFTMA2D(0,gridx(2)-gridx(1),max(gridx)-min(gridx),...
														 0,gridz(2)-gridz(1),max(gridz)-min(gridz),tmp);
		
		
		nCells = size(L,2)/2;
		
    for i=1:nsim
        
			szs = FFTMA2D(0,gridx(2)-gridx(1),max(gridx)-min(gridx),...
										0,gridz(2)-gridz(1),max(gridz)-min(gridz),Gs,mks);
			xizs = FFTMA2D(0,gridx(2)-gridx(1),max(gridx)-min(gridx),...
										 0,gridz(2)-gridz(1),max(gridz)-min(gridz),Gxi,mkxi);
			szs = szs(:);
			xizs = xizs(:);
			
			%Calcul des temps pour chaque simu
%			Lx = L(:,1:nCells);
%			Lz = L(:,(1+nCells):end);
%			ts = Lx.^2 + Lz.^2.*kron(ones(nt,1),xizs'.^2);
%			ts = kron(ones(nt,1),szs') .* sqrt(ts);
%			ts = sum(ts,2);
%			ts = ts + randn(nt,1)*covar.nugget_t;
			
			szs = [szs; xizs];
			ts = L*szs + randn(nt,1)*covar.nugget_t;
			%simulation des donnees modelisees
			s_sim=Lambda*[szs(indc);ts];
			%Simu conditionnee
			S_sim(:,i)=sco+(szs-s_sim);
        
    end
end
