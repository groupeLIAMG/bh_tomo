function [ind,t0] = plotZOP(mog, av, ap, varargin)
% plotZOP(mog, av, ap)
%

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

plotTT=false;
corrSpread=false;
tt = [];
et = [];
offsetMax = mog.data.rstepsz*0.5;
if nargin >=5
    plotTT = varargin{2};
end
if nargin >=6
	corrSpread = varargin{3};
end
if nargin >=7
	offsetMax = varargin{4};
end

if strcmp(mog.data.tunits,'ns')==1
    t0 = corr_t0(length(mog.tt), av, ap);
else
    t0 = zeros(1,length(mog.tt));
end

Tx_z = unique(mog.data.Tx_z(mog.in));
%dTx = min(abs(diff(Tx_z)));
%Tx_z = min(mog.data.Tx_z):dTx:max(mog.data.Tx_z);

if corrSpread
	lrai = sqrt( (mog.data.Tx_x-mog.data.Rx_x).^2 + ...
		(mog.data.Tx_y-mog.data.Rx_y).^2 + ...
		(mog.data.Tx_z-mog.data.Rx_z).^2 );
else
	lrai = ones(1,mog.data.ntrace);
end

ind = false(1,length(mog.data.Tx_z));
traces = zeros(length(Tx_z), mog.data.nptsptrc);
timestp = repmat(mog.data.timestp-min(t0), length(Tx_z), 1);
nos = 1:mog.data.ntrace;
for n=1:length(Tx_z)
    iTx = find(abs(Tx_z(n)-mog.data.Tx_z)<0.001);
    if ~isempty(iTx)
        i = findnear(Tx_z(n), mog.data.Rx_z(iTx));
        no = nos(iTx);
        no = no(i);
        if abs(mog.data.Rx_z(no)-Tx_z(n))<offsetMax && mog.in(no)
            traces(n,:) = lrai(no)*mog.data.rdata(:,no)';
            timestp(n,:) = mog.data.timestp-t0(no);
            ind(no) = true;
        end
    end
end
ind = find(ind);
%[ind' mog.data.Tx_z(ind)' mog.data.Rx_z(ind)']
%indS = abs(mog.data.Tx_z-mog.data.Rx_z)<(mog.data.rstepsz*0.5);


if plotTT
    tt = mog.tt(ind)-t0(ind);
    et = mog.et(ind);
end

%nz = length(Tx_z);
nt=numel(mog.data.timestp);

%prof = Tx_z;
Tx_z = repmat(Tx_z',[1 nt]);

if nargin >=4
    axes(varargin{1})
else
    figure
end

%surf(timestp,Tx_z-0.5*dTx,traces)
surf(timestp,Tx_z,traces)
shading flat
view(0,90)
axis tight
clim = caxis;
cmax = max(abs(clim));
caxis([-cmax cmax])
colormap ramac_cmap(16)
grid on

set(gca,'YDir','normal')

str = get_str_locale();

xlabel([str.s22,' [',mog.data.tunits,']'])
ylabel(str.s120)
title(mog.name,'Interpreter','none','FontSize',14)

if plotTT
    hold on
    z = mog.data.Tx_z(ind);
    plot(tt,z,'yo')
    for n=1:length(et)
        if et(n)~=-1
            plot([tt(n)-et(n) tt(n)+et(n)], [z(n) z(n)],'y')
        end
    end
    hold off
end
