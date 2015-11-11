function [diff1,diff2,diff1_min,diff2_min]=choixsimu(L,S_sim,dt,c0)
% [diff1,diff2,diff1_min,diff2_min]=choixsimu(L,S_sim,dt,c0)
  
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

diff1 = zeros(1, size(S_sim,2));
diff2 = zeros(1, size(S_sim,2));
if isempty( c0 )
	w = 1;
else
	w = max(c0)./c0;
end
for i=1:size(S_sim,2)
  diff1(i)=sum(w.*(abs(L*S_sim(:,i)-dt)));
  diff2(i)=sum(w.*(L*S_sim(:,i)-dt).^2);  
end

% figure
% subplot(2,1,1)
% plot(diff1)
% grid on
% subplot(2,1,2)
% plot(diff2)
% grid on

diff1_min=find(diff1==min(diff1));
diff2_min=find(diff2==min(diff2));
