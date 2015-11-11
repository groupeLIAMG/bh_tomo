N = 20;
z = 0:N;
z = z(:);
x = z;
y = x;

% Vitesse a z=0
a = 1;   % (km/s)
Vmax = 3; % vitesse a 20 km  (Rapport 1:3, tel que Podvin & Lecomte, 1991, fig 13)
b = (Vmax-a)/z(end);

% Source a (0,0,0)
x0 = 0;
y0 = 0;
z0 = 0;

zc = 0.5+z;
V = a+b*zc;  % Vitesse des couches
slowness = zeros(N,N,N);
for i=1:N
	for j=1:N
		for k=1:N
			slowness(i,j,k) = 1/V(k);
		end
	end
end


Tx=[x0 y0 z0];

Rx=[ones(N,1)*20 ones(N,1)*20 (1:20)'];
Tx = kron(ones(N,1), Tx);


g.xmin = 0;
g.ymin = 0;
g.zmin = 0;

g.dx = 1;
g.dy = 1;
g.dz = 1;

g.nx = N;
g.ny = N;
g.nz = N;

g.nsx = 5;
g.nsy = 5;
g.nsz = 5;

[tt, rays, L] = ttcr3d(slowness(:), g, Tx, Rx);

plot3( rays{1}(:,1), rays{1}(:,2), rays{1}(:,3) )
hold on
for n=2:length(rays)
	plot3( rays{n}(:,1), rays{n}(:,2), rays{n}(:,3) )
end
hold off
axis equal
set(gca,'ZDir','reverse')
grid on
