function [g,gt] = modeliKss3d(x0libre,covar,id,L,x,iktt,afi,lclas,c0,ax1,ax2)
% [g,gt] = modeliKss3d(x0libre,covar,id,L,x,iktt,afi,lclas,c0,ax1,ax2)

% Copyright (C) 2013 Bernard Giroux
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

modelt = covar.model;
ct = covar.c;
nugget_tt = covar.nugget_t;
nugget_lt = covar.nugget_l;

n1 = sum(sum(id{1}));
n2 = sum(sum(id{2}));
n3 = sum(id{3});
n4 = sum(id{4});


modelt(id{1})    = x0libre(1:n1);
ct(id{2})        = x0libre(n1+1:n1+n2);
nugget_tt(id{3}) = x0libre(n1+n2+1:n1+n2+n3);
nugget_lt(id{4}) = x0libre(n1+n2+n3+1:n1+n2+n3+n4);

np = length(L(1,:));
nt = length(L(:,1));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

D=covardm(x,x,modelt,ct);

if covar.aniso==0 && nugget_lt ~= 0
	D = D + nugget_lt*speye(size(D,1));
end

D = L*D*L';
if isempty(c0)==1
	D = D + nugget_tt*speye(nt);
else
	D = D + nugget_tt*sparse(diag(c0(1:nt)));
end
%%%%%%%%%%%%%%%%%%%%%%%%%%


[D,ind] = sort(D(:),1,'descend');
D = full(D);

gt = moy_bloc(D,lclas);
ind0 = find(gt<Inf);
gt = gt(ind0);

iktt2 = iktt(ind);
g = moy_bloc(iktt2,lclas);
g = g(ind0);

%%%%%%%%%%%%%%%%%%%%%%%%%%%

N = round(length(g)/afi);
g = g(1:N);
gt = gt(1:N);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if ~isempty(ax1) && ~isempty(ax2)
	gmin = min([g gt]);
	gmax = max([g gt]);
	xl = get(get(ax1,'XLabel'),'String');
	fs = get(get(ax1,'XLabel'),'FontSize');
	yl = get(get(ax1,'YLabel'),'String');
	plot(ax1,g,gt,'o',[gmin gmax],[gmin gmax],'g-.')
	set(ax1,'DataAspectRatio',[1 1 1])
	set(get(ax1,'XLabel'),'String',xl,'FontSize',fs);
	set(get(ax1,'YLabel'),'String',yl,'FontSize',fs);

	n1 = 1:length(gt);
	n2 = 1:length(g);
	xl = get(get(ax2,'XLabel'),'String');
	fs = get(get(ax2,'XLabel'),'FontSize');
	yl = get(get(ax2,'YLabel'),'String');
	plot(ax2, n2,g,'b+', n1,gt,'ro')
	set(get(ax2,'XLabel'),'String',xl,'FontSize',fs);
	set(get(ax2,'YLabel'),'String',yl,'FontSize',fs);
end
