function S = findS(v1,alpha1,v2,alpha2,epsilon_m,a,m,sigma_nf,epsilon_nf,mu_nf,epsilon_f1,sigma_f1,omega, S0)
% S = findS(v1,alpha1,v2,alpha2,epsilon_m,a,m,sigma_nf,epsilon_nf,mu_nf,epsilon_f1,sigma_f1,omega, S0)
% Estimate nano fluid saturation for the following input parameters
%
% v1:         baseline velocity (m/ns)
% alpha1:     baseline attenuation (Np/m)
% v2:         repeat velocity (m/ns)
% alpha2:     repeat attenuation (Np/m)
% epsilon_m:  relative dielectric permittivity of the rock matrix (-)
% a:          Archie's parameter
% m:          Archie's cementation exponent
% sigma_nf:   electric conductivity of the magnetic nano fluid (S/m)
% epsilon_nf: relative dielectric permittivity of the magnetic nano fluid (-)
% mu_nf:      relative magnetic permeability of the magnetic nano fluid (-)
% epsilon_f1: relative dielectric permittivity of the fluid in place (-)
% sigma_f1:   electric conductivity of the fluid in place (S/m)
% omega:      angular frequency (rad/s)
% S0:         starting value of saturation
%
% Requirement: global optimization toolbox needed to run the simulated
% annealing algorithm
%
%  B. Giroux
%  INRS-ETE
%  2014-06-12

% Copyright (C) 2014 Bernard Giroux
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

global mu0 epsilon0

% some constants
Psi = 1/3;  % depolarization factor
mu0 = 4*pi*1e-7;   % magnetic permeability of vacuum
epsilon0 = 8.8541878176203898505365630317107502606083701665994498081024171524053950954599821142852891607182008932e-12;  % dielectric permittivity of vacuum

% properties of the medium
epsilon1 = 1/mu0 * (1/(v1*1e9)^2 - alpha1^2/omega^2);
sigma1 = 2*alpha1 / (mu0*(v1*1e9));

% average the porosity
phi1 = (a*sigma1/sigma_f1)^(1/m);  % Archie
phi2 = (epsilon1-epsilon_m*epsilon0)/(epsilon_f1*epsilon0-epsilon_m*epsilon0)...
			 *(epsilon_f1*epsilon0/epsilon1)^Psi; % self-similar
phi = 0.5*(phi1+phi2);

% define objective function
f = @(S)func(S,v2,alpha2,sigma_nf,epsilon_nf,mu_nf,epsilon_m,epsilon_f1,...
						 sigma_f1,phi,Psi,a,m,omega);

% call simulated annealing algorithm
S = simulannealbnd(f, S0, 0, 1);
end

function error = func(S,v2,alpha2,sigma_nf,epsilon_nf,mu_nf,epsilon_m,...
    epsilon_f1,sigma_f1,phi,Psi,a,m,omega)

global mu0 epsilon0

epsilon_f2 = real( self_similar(epsilon_f1, epsilon_nf, S, Psi ) );
sigma_f2 = real( self_similar(sigma_f1, sigma_nf, S, Psi ) );

epsilon2 = epsilon0 * self_similar(epsilon_m, epsilon_f2, phi, Psi );
sigma2 = 1/a * sigma_f2 * phi^m;
mu2 = mu0 * self_similar( 1, mu_nf, S*phi, Psi );

alpha2_est = omega*sqrt(mu2*epsilon2*0.5*...
												(sqrt(1+(sigma2/(omega*epsilon2))^2)-1));
v2_est = 1e-9/sqrt(mu2*epsilon2*0.5*(sqrt(1+(sigma2/(omega*epsilon2))^2)+1));

error = abs(v2-v2_est)+abs(alpha2-alpha2_est);
end