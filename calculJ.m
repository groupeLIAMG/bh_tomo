function J = calculJ(L,e)
%function J = calculJ(L,e)
%
% L: [L_x L_z]
% e: [s_x; xi]    (column vector)
%

nt = size(L,1);
np = size(L,2)/2;

J = L.^2;
Js = J(:,1:np) + J(:,(np+1):end).*kron(ones(nt,1), (e((np+1):end).^2)'); % l_x^2 + l_z^2 * xi^2
Js = sqrt(Js);  % equal to t / s_x

Jxi = J(:,(np+1):end) .* kron(ones(nt,1), (e(1:np))') .* kron(ones(nt,1), (e((np+1):end))');

ind = Js~=0;
Jxi(ind) = Jxi(ind)./Js(ind);
J = [Js Jxi];
