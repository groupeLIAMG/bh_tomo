function [S2, param] = fitS(f, S, type)
% Ajuste une distribution au spectre S
%
% function [S2, param] = fitS(f, S, type)
%
% Entrée
%   f : fréquence
%   S : spectre à ajuster
%   type : 1 -> Maxwell
%          2 -> Rayleigh
%          3 -> Gamma (h=2)
%
% Sortie
%   S2 spectre ajusté
%   param
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

f = f(:);
S = S(:);

aS = sum( S*(f(2)-f(1)) );
S = S/aS;  % on normalise le spectre

f2 = f.^2;

a       = f(find(S==max(S)))^2;         % initial guess
thresh  = 0.995;                        % convergence threshold for the loop
last_cnt= inf;
iter    = 0;

S2 = [];
param = [];
if type == 1
  % Maxwell
  C = sqrt(2/pi)*f2;
  expr  = 'S2 = C*(a^(-1.5)).*exp(-f2/(2*a));';
  expr2 = 'H = S2.*( f2/(2*a^2) - 3/(2*a) );';
elseif type == 2
  % Rayleigh
  expr  = 'S2 = f.*exp( -f2/(2*a) )/a;';
  expr2 = 'H = S2.*(f2/(2*a^2) - 1/a);';
elseif type == 3
  % Gamma (h=2)
  expr  = 'S2 = f.*exp(-f/a)/(a^2);';
  expr2 = 'H = S2.*(f/a^2 - 2/a);';
else
  disp('Non implémenté')
  return
end

while (1)
  iter    = iter + 1;
  eval(expr);
  eval(expr2);
  HT      = H';
  control = inv( HT * H );
  a       = a + control * HT * (S-S2);
  if ( control>(last_cnt * thresh) )
	break;
  else
	last_cnt = control;
  end
end

% summarize results
eval(expr);
S2 = S2*aS;
eval(expr2);
err        = ( S - S2 );
param.a    = a;
param.VAR  = inv( (H') * H );
param.RMS  = sqrt( (err')*err/ (f(2)-f(1))^2 / (length(err)-1) );
param.iter = iter;
param.type = type;

