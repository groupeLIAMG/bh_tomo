function [data,ind_max] = filtrage_wavelet2(rdata)
% function data = filtrage_wavelet2(rdata)
%
%  rdata: traces a filtrer
%

% Copyright (C) 2005 Bernard Giroux, Abderezzak Bouchedda
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


[nptsptrc,ntrace] = size(rdata);
data = zeros(size(rdata));
N=3;
npts = ceil(nptsptrc/2^N)*2^N;
d=npts-nptsptrc;
rdata = [rdata; rdata(end-d+1:end,:)];
ondelette = 'db2';
ind_max= zeros(ntrace,1);

h=waitbar(0,'Data denoising ...');
inc_wb = round(ntrace/100);
for n=1:ntrace
	%  trace = detrend_rad(rdata(1:npts,n));
	%trace2 = debruite(rdata, ondelette, N);
	trace2 = debruite(rdata(:,n), ondelette, N);
	trace2 = medfilt1(trace2,5);
	data(:,n) = trace2(1:nptsptrc);
    ind1=find((data(:,n))== max((data(:,n))));
    ind_max(n)=ind1(1);
	if rem(n, inc_wb)==0, waitbar(n/ntrace,h); end
end
close(h);


function trace_filtree = debruite(trace, ondelette, N)

swc=swt(trace,N,ondelette);
%th = wthrmngr('sw1ddenoLVL','rigrsure',swc,'sln');
th = wthrmngr('sw1ddenoLVL','sqtwolog',swc,'mln');
%fac = ones(1,N);
%fac(N-1) = 0.25;
%fac(N) = 0.25;
%th = th.*fac;
for n=1:N
	swc(n,:) = wthresh(swc(n,:), 'h', th(n));
end
trace_filtree=iswt(swc,ondelette);
