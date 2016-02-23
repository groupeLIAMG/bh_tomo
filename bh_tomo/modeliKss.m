function [g,gt] = modeliKss(x0libre,covar,id,L,gridx,gridz,x,iktt,afi,lclas,c0,ax1,ax2)
% [g,gt] = modeliKss(x0libre,covar,id,L,gridx,gridz,x,iktt,afi,lclas,c0,ax1,ax2)

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

modelt = covar.model;
ct = covar.c;
nugget_tt = covar.nugget_t;
nugget_lt = covar.nugget_l;

n1 = sum(sum(id{1}));
n2 = sum(sum(id{2}));
n3 = sum(id{3});
n4 = sum(id{4});

if covar.aniso > 0
    model_xi_t = covar.model_xi;
    c_xi_t = covar.c_xi;
    nugget_xi_t = covar.nugget_xi;
    n5 = sum(sum(id{5}));
    n6 = sum(sum(id{6}));
    n7 = sum(sum(id{7}));
end
if covar.aniso > 1
    model_theta_t = covar.model_th;
    c_theta_t = covar.c_th;
    nugget_th_t = covar.nugget_th;
    n8 = sum(sum(id{8}));
    n9 = sum(sum(id{9}));
    n10 = sum(sum(id{10}));
end

modelt(id{1})    = x0libre(1:n1);
ct(id{2})        = x0libre(n1+1:n1+n2);
nugget_tt(id{3}) = x0libre(n1+n2+1:n1+n2+n3);
nugget_lt(id{4}) = x0libre(n1+n2+n3+1:n1+n2+n3+n4);

if covar.aniso > 0
    model_xi_t(id{5})  = x0libre(n1+n2+n3+n4+1:n1+n2+n3+n4+n5);
    c_xi_t(id{6})      = x0libre(n1+n2+n3+n4+n5+1:n1+n2+n3+n4+n5+n6);
    nugget_xi_t(id{7}) = x0libre(n1+n2+n3+n4+n5+n6+1:n1+n2+n3+n4+n5+n6+n7);
end
if covar.aniso > 1
    model_theta_t(id{8}) = x0libre(n1+n2+n3+n4+n5+n6+n7+1:n1+n2+n3+n4+n5+n6+n7+n8);
    c_theta_t(id{9})     = x0libre(n1+n2+n3+n4+n5+n6+n7+n8+1:n1+n2+n3+n4+n5+n6+n7+n8+n9);
    nugget_th_t(id{10})  = x0libre(n1+n2+n3+n4+n5+n6+n7+n8+n9+1:n1+n2+n3+n4+n5+n6+n7+n8+n9+n10);
end


np = length(L(1,:));
nt = length(L(:,1));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

D = calculD(x, modelt, ct, length(gridx), length(gridz));
if covar.aniso==1
    Dxi = calculD(x, model_xi_t, c_xi_t, length(gridx), length(gridz));
    D = sparse([D zeros(size(Dxi)); zeros(size(D)) Dxi]);
elseif covar.aniso==2
    Dxi = calculD(x, model_xi_t, c_xi_t, length(gridx), length(gridz));
    Dth = calculD(x, model_theta_t, c_theta_t, length(gridx), length(gridz));
    D = sparse([D              zeros(size(Dxi)) zeros(size(Dth));
                zeros(size(D)) Dxi              zeros(size(Dth));
                zeros(size(D)) zeros(size(Dxi)) Dth]);
end

if covar.aniso==0 && nugget_lt ~= 0
	D = D + nugget_lt*speye(size(D,1));
elseif covar.aniso==1 && ( nugget_lt ~= 0 || nugget_xi_t ~= 0 )
    np = np/2;
    D = D + sparse([nugget_lt*speye(np)   zeros(np);
                    zeros(np) nugget_xi_t*speye(np)]); %%YH zeros(no) --> np
elseif covar.aniso==2 && ( nugget_lt ~= 0 || nugget_xi_t ~= 0 || nugget_th_t ~= 0 )
    np = np/3;
    D = D + sparse([nugget_lt*speye(np)   zeros(np) zeros(np);
                    zeros(np) nugget_xi_t*speye(np) zeros(np);  %%YH first zeros(no) --> np
                    zeros(np) zeros(np) nugget_th_t*speye(np)]);
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
