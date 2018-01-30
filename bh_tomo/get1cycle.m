function [out,fdom,w,taper,idebut] = get1cycle(trace, b, a, dt, f0, seuil_selection, fmin, fmax)

%npt = 10;
split=4;
if isempty(a) || isempty(b)
	fdom = get_fdom(trace, dt, f0, fmin, fmax);
else
	fdom = get_fdom(filter(b,a,trace), dt, f0, fmin, fmax);
end
if fdom == -1
    fdom2 = f0;
else
    fdom2 = fdom;
end
ne = length(trace);

m = 2;
n = 2;
t = 0:dt:n/fdom2;
w = warcone(t,fdom2,0,m,n);

%w = Ricker(dt,fdom2);

xc = xcorr(trace,w)';
xc = xc(ne:end) * (1/max(abs(xc(ne:end))));
%save toto xc
[~,imax] = extr(xc);
i2max = xc(imax) > seuil_selection;
imax = imax(i2max);
if isempty(imax)
	out = trace;
	return
end
idebut = imax(1);
npt = round(length(w)/split);
ifin = round(idebut+length(w)-npt);

taper = zeros(size(trace));
taper(idebut:ifin) = 1;
d = 1/(npt) * (1:npt);

ct = 0.5-0.5*cos(pi*d);
if idebut-npt > 1
	taper(idebut-npt:idebut-1) = ct;
end
taper(ifin+1:ifin+npt) = fliplr(ct);
if numel(trace) ~= numel(taper)
    out = trace;
    return
end
out = trace.*taper;

