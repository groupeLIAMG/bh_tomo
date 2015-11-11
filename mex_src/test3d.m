
N1 = 5;
N2 = 5;
N3 = 5;
z = 0:N3;
z = z(:);
x = 0:N1;
x = x(:);
y = 0:N2;
y = y(:);

% Vitesse a z=0
a = 1;   % (km/s)
Vmax = 3; % vitesse a 20 km  (Rapport 1:3, tel que Podvin & Lecomte, 1991, fig 13)
b = (Vmax-a)/z(end);

% Source a (0,0,0)
x0 = 0;
y0 = 0;
z0 = 0;

slowness = ones(N1,N2,N3);


Tx=[x0 y0 z0];

Rx=[ones(N3,1)*N1 ones(N3,1)*N2 (1:5)'];
Tx = kron(ones(N3,1), Tx);


g.xmin = 0;
g.ymin = 0;
g.zmin = 0;

g.dx = 1;
g.dy = 1;
g.dz = 1;

g.nx = N1;
g.ny = N2;
g.nz = N3;

g.nsx = 10;
g.nsy = 10;
g.nsz = 10;

[tt, rays, L1] = ttcr3d(slowness(:), g, Rx, Tx);


[L2, gridx, gridy, gridz] = Lsr3d(Rx, Tx, x, y, z);

figure(1)
imagesc(L1); colorbar
figure(2)
imagesc(L2); colorbar

err = (full(sum(L2,2))-full(sum(L1,2)))'

l2 = sqrt(sum( (Tx-Rx).^2, 2));

for n=1:length(rays)
	l1(n,1) = sum(sqrt(sum(diff(rays{n}).^2,2)));
end

(l2 - l1)'


tt2 = L2*slowness(:);


