

g.xmin = 0;
g.zmin = 0;
g.dx = 0.2;
g.dz = 0.2;
g.nx = 100;
g.nz = 150;
g.nsnx = 5;
g.nsnz = 5;


x=g.dx*(2:g.nx-1);
z=g.dz*(2:g.nz-1);
Rx = [kron(x(:),ones(numel(z),1)) kron(ones(numel(x),1),z(:))];
Tx = kron([0.5*g.dx*g.nx 0.0 0.5*g.dz*g.nz],ones(size(Rx,1),1));

Rx = [Rx(:,1) zeros(size(Rx,1),1) Rx(:,2)];


slowness = ones(g.nx*g.nz,1);
xi = 2*slowness;
theta = 30*slowness*pi/180;


t = 'iso';

g2 = grid2d_mex('new',g,t,2);

gx = grid2d_mex('get_dx',g2);
nz = grid2d_mex('get_nz',g2);

strcmp(grid2d_mex('get_type',g2), t)

t1 = grid2d_mex('raytrace',g2, slowness, Tx, Rx);
[tt1,r1,L1] = grid2d_mex('raytrace',g2, slowness, Tx(1:10,:), Rx(1:10,:));

grid2d_mex('delete',g2)



t = 'elliptical';
g2 = grid2d_mex('new',g,t,2);
t2 = grid2d_mex('raytrace',g2, slowness, xi, Tx, Rx);
grid2d_mex('delete',g2)

t = 'tilted';
g2 = grid2d_mex('new',g,t,2);
t3 = grid2d_mex('raytrace',g2, slowness, xi, theta, Tx, Rx);
grid2d_mex('delete',g2)




clear g2

figure(1)
imagesc(reshape(t1,g.nz-2,g.nx-2))
axis equal, axis tight

figure(2)
imagesc(reshape(t2,g.nz-2,g.nx-2))
axis equal, axis tight

figure(3)
imagesc(reshape(t3,g.nz-2,g.nx-2))
axis equal, axis tight

