function       [S_sim,sco]=cokri_simua(L,dt,covar,nt,gridx,gridz,nsim,mode,c0)
% [S_sim,ssco]=cokri_simua(L,dt,covar,indi,nt,gridx,gridz,nsim,mode,ind,ind2,c0)

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


global kss

if isempty(c0)==1
   ktt = L*kss*L' + covar.nugget_t*eye(nt);
else
   ktt = L*kss*L' + diag(covar.nugget_t.*c0(1:nt));
end

kts = L*kss;

% cokrigeage primal
%tic
Lambda = (ktt\kts)';
sco = Lambda*dt;
%toc

% cokrigeage dual
%tic
%Gamma = (ktt\dt)';
%sco = (Gamma*kts)';
%toc



S_sim=0;
if mode==2
	S_sim=zeros(length(sco),nsim);
	
	[Gs, mks] = preFFTMA2D(0,gridx(2)-gridx(1),max(gridx)-min(gridx),...
												 0,gridz(2)-gridz(1),max(gridz)-min(gridz),covar);
	tmp.model = covar.model_xi;
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
%		Lx = L(:,1:nCells);
%		Lz = L(:,(1+nCells):end);
%		ts = Lx.^2 + Lz.^2.*kron(ones(nt,1),xizs'.^2);
%		ts = kron(ones(nt,1),szs') .* sqrt(ts);
%		ts = sum(ts,2);
%		ts = ts + randn(nt,1)*covar.pepite_t;
		
		szs = [szs; xizs];
		ts = L*szs + randn(nt,1)*covar.nugget_t;
				
		s_sim=Lambda*ts;
		
		S_sim(:,i)=sco+(szs-s_sim);
	end
end


