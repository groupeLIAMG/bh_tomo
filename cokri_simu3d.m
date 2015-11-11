function [S_sim,sco]=cokri_simu3d(L,dt,covar,nt,gridx,gridy,gridz,nsim,mode,c0)
% [S_sim,ssco]=cokri_simu3d(L,dt,covar,indi,nt,gridx,gridy,gridz,nsim,mode,ind,ind2,c0)

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




S_sim=0;
%if mode==2 && covar.aniso==0
if mode==2    
	% cokrigeage primal
	Lambda = (ktt\kts)';
	sco = Lambda*dt;
	
	
	S_sim=zeros(length(sco),nsim);
	
	[G, mk] = preFFTMA3D(0,gridx(2)-gridx(1),max(gridx)-min(gridx),...
											 0,gridy(2)-gridy(1),max(gridy)-min(gridy),...
											 0,gridz(2)-gridz(1),max(gridz)-min(gridz),covar);
	for i=1:nsim
		
		szs=FFTMA3D(0,gridx(2)-gridx(1),max(gridx)-min(gridx),...
								0,gridy(2)-gridy(1),max(gridy)-min(gridy),...
								0,gridz(2)-gridz(1),max(gridz)-min(gridz),G,mk);
		
		szs=reshape(szs,numel(szs),1);
		
		ts=L*szs + randn(nt,1)*covar.nugget_t;
		
		s_sim=Lambda*ts;
		
		S_sim(:,i)=sco+(szs-s_sim);
	end
else
    
    % cokrigeage dual
    %tic
    Gamma = (ktt\dt)';
    sco = (Gamma*kts)';
    %toc

end


