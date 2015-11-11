function [S_sim,sco]=simulation_cond(L,dt,covar,nt,gridx,gridz,nsim,mode,c0,varargin)
% [S_sim,ssco]=simulation_cond(L,dt,covar,indi,nt,gridx,gridz,nsim,mode,ind,ind2,c0,varargin)

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

handle = [];
if nargin==10
	handle = varargin{1};
end
if isempty(handle)
    disp('Inversion par cokrigeage - Cokrigeage et/ou simulation sans contraintes')
else
  str = get_str_locale();
  set(handle,'String',[str.s225,' - ',str.s234])
  drawnow
end

global kss

if isempty(c0)==1
   ktt2=L*kss*L'+covar.pepite_t*eye(nt);
else
%   ktt2=L*kss*L'+ca(3)*eye(nt);
   ktt2=L*kss*L'+diag(covar.pepite_t.*c0(1:nt));
end

kts=L*kss;

%Lambda=zeros(nt,np);
%Lambda=(inv(ktt2)*kts)';
Lambda=(ktt2\kts)';

sco=Lambda*dt;

S_sim=0;
if mode==2
	S_sim=zeros(length(sco),nsim);
	if isempty(handle)
		disp('Debut des simulations')
	else
		set(handle, 'String', [str.s225,' - ',str.s233])
	end
	
	[G, mk] = preFFTMA2D(0,gridx(2)-gridx(1),max(gridx)-min(gridx),...
				0,gridz(2)-gridz(1),max(gridz)-min(gridz),covar);
	for i=1:nsim

		szs=FFTMA2D(0,gridx(2)-gridx(1),max(gridx)-min(gridx),...
				0,gridz(2)-gridz(1),max(gridz)-min(gridz),G,mk);
%		[szs]=fftma(0,gridx(2)-gridx(1),max(gridx)-min(gridx),....
%			0,gridz(2)-gridz(1),max(gridz)-min(gridz),covar);
		
		szs=reshape(szs,numel(szs),1);

		ts=L*szs + randn(nt,1)*covar.pepite_t;
		
		s_sim=Lambda*ts;
		
		S_sim(:,i)=sco+(szs-s_sim);
	end
end


