function [diff1,diff2,diff1_min,diff2_min]=choixsimuaa(L,S_sim,dt,c0)
% [diff1,diff2,diff1_min,diff2_min]=choixsimua(L,S_sim,dt,c0)

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

nCells = size(L,2)/2;
nt = length(dt);
Lx = L(:,1:nCells);
Lz = L(:,(1+nCells):end);

for i=1:size(S_sim,2)
    th = S_sim((1+2*nCells):end,i);
    xi = S_sim((1+nCells):2*nCells,i);
    s = S_sim(1:nCells, i);
    co = kron(ones(nt,1),cos(th)');
    si = kron(ones(nt,1),sin(th)');
    t = (Lx.*co + Lz.*si).^2 + (Lz.*co - Lx.*si).^2.*kron(ones(nt,1),xi'.^2);
    t = kron(ones(nt,1),s') .* sqrt(t);
    t = sum(t,2);
    
    diff1(i)=sum(w.*(abs(t-dt)));
    diff2(i)=sum(w.*(t-dt).^2);
end

diff1_min=find(diff1==min(diff1));
diff2_min=find(diff2==min(diff2));
