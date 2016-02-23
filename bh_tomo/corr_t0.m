function [t0, fac_dt_av, fac_dt_ap] = corr_t0(ndata, av, ap, varargin)
% t0 = corr_t0(ndata, av, ap)
%
%
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

if nargin>=4
	show = varargin{1};
else
	show = false;
end

fac_dt_av = 1;
fac_dt_ap = 1;
if isempty(av) && isempty(av)
	str = get_str_locale();
	uiwait(warndlg(str.s59))
	t0 = zeros(1,ndata);
	return
end

v_air = 0.2998;
t0av = [];
t0ap = [];
if ~isempty(av)
	if strcmp(av.method, 'fixed_antenna')
		t0av = get_t0_fixed(av, v_air);
	elseif strcmp(av.method, 'walkaway')
		[t0av, fac_dt_av] = get_t0_wa(av, v_air, show);
	end
end
if ~isempty(ap)
	if strcmp(ap.method, 'fixed_antenna')
		t0ap = get_t0_fixed(ap, v_air);
	elseif strcmp(ap.method, 'walkaway')
		[t0ap, fac_dt_ap] = get_t0_wa(ap, v_air, show);
	end
end


if isnan(t0av) || isnan(t0ap)
	str = get_str_locale();
	uiwait(warndlg(str.s59))
	t0 = zeros(1,ndata);
	return
end

if isempty(t0av) && isempty(t0ap)
	t0 = zeros(1,ndata);
elseif isempty(t0av)
	t0 = t0ap*ones(1,ndata);
elseif isempty(t0ap)
	t0 = t0av*ones(1,ndata);
else
	dt0 = t0ap-t0av;
	ddt0 = dt0/(ndata-1);
	t0 = t0av+ddt0*(0:(ndata-1));
end

% /\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\
%
%
function t0 = get_t0_fixed(tir, v)
tt = tir.tt(tir.tt_done);
et = tir.et(tir.tt_done);
ind = tt~=-1;
if all(et == -1)
	tt = mean(tt(ind));
else
	tt = sum(tt(ind).*et(ind))/sum(et(ind));
end
t0 = tt - tir.d_TxRx/v;

% /\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\
%
%
function [t0, fac] = get_t0_wa(tir, v, show)
ind = tir.tt~=-1;
tt = tir.tt(tir.tt_done & ind)';
et = tir.et(tir.tt_done & ind)';
d = tir.d_TxRx(tir.tt_done & ind)';

slown = 1/v;

if all(et == -1)
	b = [d ones(size(d))]\tt;
	t0 = b(2);
	%app_slowness = b(1);
	%fac = true_slowness/app_slowness;
	fac = slown/b(1);

	if show
		figure('Name','Air shot')
		plot(d,tt,'o')
		hold on
		plot([0;d],b(1)*[0;d]+b(2),'k')
		xlabel('Distance')
		ylabel('Time')
		title([tir.name,' - correction factor: ',num2str(fac)])
		text(d(2), b(1), ['t_0 at ',num2str(t0)])
		hold off
	end
else
	W = diag(1./(et.^2));
	x = [d ones(size(d))];
	b = (x'*W*x)\(x'*W*tt);
	t0 = b(2);
	fac = slown/b(1);

	if show
		figure('Name','Air shot')
		subplot(121)
		plot([0;d],b(1)*[0;d]+b(2),'k','LineWidth',1)
		hold on
		errorbar(d,tt,et,'o')
		xlabel('Distance')
		ylabel('Time')
		title([tir.name,' - correction factor: ',num2str(fac)])
		text(d(2), b(1), ['t_0 at ',num2str(t0)])
        ylim=get(gca,'YLim');
        ylim(1)=0;
        set(gca,'YLim',ylim)
		hold off
		subplot(122)
		plot([0;d],slown*[0;d]+b(2)*fac,'g','LineWidth',1)
		hold on
		errorbar(d,tt*fac,et,'o')
		xlabel('Distance')
		ylabel('Time')
		title('After \Delta t corection')
		text(d(2), b(1), ['t_0 at ',num2str(t0*fac)])
        ylim=get(gca,'YLim');
        ylim(1)=0;
        set(gca,'YLim',ylim)
		hold off
	end
end


