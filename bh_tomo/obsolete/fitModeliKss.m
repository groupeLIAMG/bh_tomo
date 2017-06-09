function [dif]=fitModeliKss(x0libre,covar,id,L,gridx,gridz,x,iktt,afi,lclas,c0,ax1,ax2)
% [dif]=fitModeliKss(x0libre,covar,id,L,gridx,gridz,x,iktt,afi,lclas,c0,ax1,ax2)

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

weight=0;

[g,gt] = modeliKss(x0libre,covar,id,L,...
    gridx,gridz,x,iktt,afi,lclas,c0,ax1,ax2);

q = fliplr((1:length(g)).^2);
q = 1 + q./max(q);
%q=round(q(1:floor(length(q)/afi)));
%q=q(1:floor(length(q)/afi));
g = g(1:length(q));
gt = gt(1:length(q));

if weight == 0
    dif = sum((((gt-g).*q).^2));
else
    bign = 99999;
    
    %  a=(c(1)<0)*bign*abs(c(1));
    %  a=a^2;
    a=0;
    
    %  b=(c(1)>1e8)*bign*abs(c(1)-1e8);
    b=0;
    cc=(c(1)<0)*bign*abs(c(1));
    cc=cc^2;
    d=(c(1)>1e8)*bign*abs(c(1)-1e8);
    e=(pepite_t<0)*bign*abs(pepite_t);
    e=e^2;
    f=(pepite_t>1e4)*bign*abs(pepite_t-1e4);
    
    % mm=(c(4)<0)*bign*abs(c(3));
    mm = 0;
    
    gg = (covar.model(1,2)<0)*bign*abs(covar.model(1,2));
    mxz = max(max(gridx),max(gridz));
    h = (covar.model(1,2)>1*mxz)*bign*abs(covar.model(1,2)-1*mxz);
    if length(covar.model(1,:))>3
        i = (covar.model(1,3)<0)*    bign*abs(covar.model(1,3));
        j = (covar.model(1,3)>1*mxz)*bign*abs(covar.model(1,3)-1*mxz);
        k = (covar.model(1,4)<0)*    bign*abs(covar.model(1,4));
        l = (covar.model(1,4)>180)*  bign*abs(180-covar.model(1,4));
    else
        i=0;j=0;
        k=(covar.model(1,2)<0)*    bign*abs(covar.model(1,2));
        l=(covar.model(1,2)>1*mxz)*bign*abs(1*mxz-covar.model(1,2));
    end
    dif=sum((((gt(:,2)-g(:,2)).*q).^2)+bign*abs(sum((x0libre<0).*x0libre)))+ ...
        a+b+cc+d+e+f+gg+h+i+j+k+l;
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
