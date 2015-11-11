function [L, rais, tt, corr_tt]=calcul_rc_ant2(s,gp,Tx,Rx,corr,diam,TxCosDir,RxCosDir, inWater, aniso, varargin)
%
%
%

% Copyright (C) 2008 Bernard Giroux
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

if nargin >=11
    xi = varargin{1};
end
ndata = size(Tx,1);

thetaTx = acos(dot( kron([0 0 -1], ones(ndata,1)), TxCosDir, 2));
signe = cross( kron([0 0 -1], ones(ndata,1)),	TxCosDir, 2);
thetaTx = thetaTx.*sign(signe(:,2));

thetaRx = acos(dot( kron([0 0 -1], ones(ndata,1)), RxCosDir, 2));
signe = cross( kron([0 0 -1], ones(ndata,1)), RxCosDir, 2);
thetaRx = thetaRx.*sign(signe(:,2));

typeC = get(corr,'Type');

if aniso==0
    [tt, rais, L, corr_tt] = ttcr2d_ant(s, gp, Tx, Rx, diam, [thetaTx thetaRx], inWater, typeC);
else
    [tt, rais, L, corr_tt] = ttcr2da_ant(s, xi, gp, Tx, Rx, diam, [thetaTx thetaRx], inWater, typeC);
end
