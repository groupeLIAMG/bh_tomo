function [covara,dif]=ajuster_Cov3d(options,covar,covar_code,L,...
	gridx,gridz,x,iktt,afi,lclas,c0,ax1,ax2)
%  fonction pour ajuster simultanement un modele dans plusieurs directions

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

% verifier si l'on est avec un cas 2D anisotrope ou 3D anisotrope
%[n,p] = size(covar.model);

id{1} = covar_code.model==0;
p1 = sum(sum(id{1}));
x0libre = covar.model(id{1});
x0libre = x0libre(:);  % on s'assure qu'on a un vecteur colonne

id{2} = covar_code.c==0;
p2 = sum(sum(id{2}));
x0libre = [x0libre; covar.c(id{2})];

id{3} = covar_code.nugget_t==0;
p3 = sum(id{3});
x0libre = [x0libre; covar.nugget_t(id{3})];

id{4} = covar_code.nugget_l==0;
p4 = sum(id{4});
x0libre = [x0libre; covar.nugget_l(id{4})];

if ~isempty(x0libre)
	[t1,dif,flag]=fminsearch('fitModeliKss3d',x0libre,options,covar,...
		id,L,x,iktt,afi,lclas,c0,ax1,ax2);
end

if p1>0;
	t11 = t1(1:p1);
	covar.model(id{1}) = t11;
end
if p2>0
	t21 = t1(p1+1:p1+p2);
	covar.c(id{2}) = t21;
end
if p3>0
	t31 = t1(p1+p2+1:p1+p2+p3);
	covar.nugget_t(id{3}) = t31;
end
if p4>0
	t41 = t1(p1+p2+p3+1:p1+p2+p3+p4);
	covar.nugget_l(id{4}) = t41;
end
covara = covar;
