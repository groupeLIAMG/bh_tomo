
clear

grx=-1:0.5:4;
gry=0:0.5:4;
grz=-15:0.5:1;

g = Grid3D(grx,gry,grz);

s = ones(numel(grx)-1,numel(gry)-1,numel(grz)-1);

Tx = [0 2 0; 0 2 0;  0 2 0;  0 2 0;  0 2 0];
Rx = [3 2 0; 3 2 -1; 3 2 -2; 3 2 -3; 3 2 -4];

[tt, rays, L] = g.raytrace(s(:), Tx, Rx);

grz = -16:0.5:1;
g.grz = grz;
s = ones(numel(grx)-1,numel(gry)-1,numel(grz)-1);

[tt2, rays2, L2] = g.raytrace(s(:), Tx, Rx);

tt-tt2



g2 = Grid2D(grx,grz);

s = ones(numel(grx)-1,numel(grz)-1);

Tx = [0 0; 0 0;  0 0;  0 0;  0 0];
Rx = [3 0; 3 -1; 3 -2; 3 -3; 3 -4];

[tt3, rays, L] = g2.raytrace(s(:), Tx, Rx);

grz = -16:0.5:1;
g2.grz = grz;
s = ones(numel(grx)-1,numel(grz)-1);

[tt4, rays2, L2] = g2.raytrace(s(:), Tx, Rx);

tt3-tt4

save /tmp/test

