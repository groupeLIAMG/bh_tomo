function ha=plotRais2D(t,varargin)
% h=plotRais2D(t)

% Copyright (C) 2005 Bernard Giroux
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

plotRes = 0;
if nargin>=2
	plotRes = varargin{1};
end
dim = size(t.rays{1},2);
if plotRes==0
	h=plot(t.rays{1}(:,1), t.rays{1}(:,dim),'k');
	ha = get(h, 'Parent');
	hold(ha,'on')
	for n=2:length(t.rays)
		plot(ha,t.rays{n}(:,1), t.rays{n}(:,dim),'k')
	end
	hold(ha,'off')
	set(ha,'DataAspectRatio',[1 1 1])
else
	res = t.invData(end).res;
	
	rmin = min(res);
	rmax = max(res);
	%c=jet;
    
    c = [0 0 1;0.8 0.8 0.8;1 0 0];
    c = interp1((-1:1)',c,(-1:0.02:1)');

	m = (size(c,1)-1)/(rmax-rmin);
	b = 1-rmin*m;
	p = m*res(1)+b;
	couleur = interp1(c,p);

	h=plot(t.rays{1}(:,1), t.rays{1}(:,dim), 'Color',couleur);
	ha = get(h, 'Parent');
	hold(ha,'on')
	for n=2:length(t.rays)
		p = m*res(n)+b;
		couleur = interp1(c,p);
		plot(ha,t.rays{n}(:,1), t.rays{n}(:,dim), 'Color',couleur)
	end
	hold(ha,'off')
	set(ha, 'DataAspectRatio',[1 1 1])
	colormap(c)%jet)
	hb=colorbar(ha);
	caxis(ha,[rmin rmax])
	set(get(hb,'Title'),'String','Residuals')
    set(gcf,'renderer','opengl'); 
end