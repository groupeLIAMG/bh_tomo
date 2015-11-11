function J = calculJ2(L,e)
%function J = calculJ2(L,e)
%
% L: [L_x L_z]
% e: [s_x; xi; theta]    (column vector)
%

nt = size(L,1);
np = size(L,2)/2;

n = length(e)/3;

J = [];
if np ~= n
    disp('Erreur (calculJ2) - tailles de L et e incompatibles')
    return
end

s = (e(1:n))';
xi = (e((n+1):(2*n)))';
theta = (e((2*n+1):(3*n)))';

co = cos(theta);
si = sin(theta);

%
%    t_j = s_{x',j}\sqrt{\left(l_{x,j}\cos\theta_j+l_{z,j}\sin\theta_j\right)^2 +
%    \xi_{j}^2\left(l_{z,j}\cos\theta_j-l_{x,j}\sin\theta_j\right)^2},

%
%    \mathbf{J}_{s,ij} &{}={}& \frac{\partial t_{ij}(\mathbf{e}_0)}{\partial s_{x',j}} = 
%      \frac{t_{ij}}{s_{0\, x',j}}

tmp = (L(:,1:np) .* kron(ones(nt,1), co) + L(:,(np+1):end) .* kron(ones(nt,1), si)).^2;
tmp = tmp + kron(ones(nt,1), xi.^2) .* ...
    (L(:,1:np) .* kron(ones(nt,1), si) - L(:,(np+1):end) .* kron(ones(nt,1), co)).^2;
Js = sqrt(tmp);


%
%    \mathbf{J}_{\xi,ij} &{}={}& \frac{\partial t_{ij}(\mathbf{e}_0)}{\partial \xi_{j}} = 
%      \frac{s_{0\, x',j}^2\, \xi_{0\, j}\left[l_x\sin \left(\theta_{0\, j} \right) -l_z\cos \left(\theta_{0\, j} \right)\right]^2}{t_{ij}}
  
tmp = (L(:,1:np) .* kron(ones(nt,1), si) - L(:,(np+1):end) .* kron(ones(nt,1), co)).^2;
Jxi = kron(ones(nt,1), s) .* kron(ones(nt,1), xi) .* tmp;

ind = Js~=0;
Jxi(ind) = Jxi(ind)./Js(ind);

%
%    \mathbf{J}_{\theta,ij} &{}={}& \frac{\partial t_{ij}(\mathbf{e}_0)}{\partial \theta_{j}} = 
%      \frac{s_{0\, x',j}^2 \left(\xi_{0\, j} ^2-1\right) \left[l_x^2\sin \left(2 \theta_{0\, j} \right) 
%      -2 l_z l_x\cos \left(2 \theta_{0\, j} \right) -l_z^2\sin \left(2 \theta_{0\, j} \right) \right]}{2 t_{ij}}
tmp = L(:,1:np).^2 - L(:,(np+1):end).^2;
tmp = tmp .* kron(ones(nt,1), sin(2*theta));
tmp = tmp - 2 * L(:,1:np) .* L(:,(np+1):end) .* kron(ones(nt,1), cos(2*theta));
Jtheta = kron(ones(nt,1), s) .* kron(ones(nt,1), (xi.^2-1)) .* tmp;

Jtheta(ind) = Jtheta(ind)./Js(ind);

J = [Js Jxi Jtheta];
