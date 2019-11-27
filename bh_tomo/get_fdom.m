function [fdom, Pxx] = get_fdom(trace, dt, f0, varargin)
% fdom = get_fdom(trace, dt, f0, l, u)
%
% Trouve la fréquence dominante d'une trace radar
%
% input
%   trace : trace radar à traiter
%   dt    : pas d'échantillonnade de la trace
%   f0    : fréquence nominale des antennes
%   l     : limite inférieure acceptable (optionnel, par défaut l=0.5*f0)
%   u     : limite supérieure acceptable (optionnel, par défaut u=1.1*f0)
%
%
% output
%   fdom  : la fréquence dominante trouvée.  Si fdom<l ||
%            fdom>u, la fréquence retournée est -1
%
%

lower = 0.5*f0;
upper = 1.1*f0;
if nargin >= 4
	if ~isempty(varargin{1})
		lower = varargin{1};
	end
end
if nargin >= 5
	if ~isempty(varargin{2})
		upper = varargin{2};
	end
end
ordre = 5;
if nargin >= 6
    estimator = varargin{3};
else
    nfft = 2^(1+nextpow2(length(trace)));
    ordre = 2;
    estimator = @(x) pburg(x,ordre,nfft,1/dt);
end

[Pxx,freq] = estimator(trace);
[~,ifreq] = max(Pxx);
fdom = freq(ifreq);
while ( fdom<lower || fdom>upper ) && ordre<5
    ordre = ordre+1;
    [Pxx,freq] = pburg(trace,ordre,nfft,1/dt);
    [~,ifreq] = max(Pxx);
    fdom = freq(ifreq);
end
if fdom == 0 || ( ordre == 4 && ( fdom<lower || fdom>upper ) )
  fdom = -1;
end
