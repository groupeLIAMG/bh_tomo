

g.xmin = 0;
g.zmin = 0;
g.dx = 0.2;
g.dz = 0.2;
g.nx = 21;
g.nz = 41;
g.nsx = 5;
g.nsz = 5;

s = ones(g.nx*g.nz, 1);

z = 1:0.4:7;
z = z(:);

Tx = [0.2*ones(size(z)) z];
Rx = [3.8*ones(size(z)) z];

Tx = kron(Tx, ones(length(z), 1));
Rx = kron(ones(length(z), 1), Rx);

mex -largeArrayDims ttcr2d.cpp
tic
[tt, rays, L] = ttcr2d(s, g, Tx, Rx);
toc


mex -largeArrayDims CXXFLAGS="\$CFLAGS -fopenmp" LDFLAGS="\$LDFLAGS -fopenmp" ttcr2d.cpp 
tic
[tt2, rays2, L2] = ttcr2d(s, g, Tx, Rx);
toc

sum(tt-tt2)
sum(L-L2)
