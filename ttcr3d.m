% function [tt, rays, L] = ttcr3d(s, g, Tx, Rx)
% TTCR2D - Travel times with curved rays in 3D
%
% Input
%    s: slowness vector ( nSlowness by 1 )
%    g: structured variable defining the grid parameters
%          xmin: origin in X
%          ymin: origin in Y
%          zmin: origin in Z
%          dx: cell size in X
%          dy: cell size in Y
%          dz: cell size in Z
%          nx: number of cells in X
%          ny: number of cells in Y
%          nz: number of cells in Z
%          nsx: number of secondary nodes per cell in X
%          nsy: number of secondary nodes per cell in Y
%          nsz: number of secondary nodes per cell in Z
%    Tx: transmitter coordinates, nTx by 3
%          1st column contains X coordinates, 2nd contains Y coordinates, 3rd contains Z coordinates
%    Rx: receiver coordinates, nRx by 3
%          1st column contains X coordinates, 2nd contains Y coordinates, 3rd contains Z coordinates
%
%    *** nTx must be equal to nRx, i.e. each row define one Tx-Rx pair ***
%    *** nSlowness must equal g.nx*g.ny*g.nz ***
%    *** Indexing of slowness values is done by "vectorizing" a 3D array,
%        i.e. if slowness field s is of size (nx,ny,nz), enter s(:) as
%        first argument
%
%
%
% Output
%    tt: vector of traveltimes, nRx by 1
%    rays: cell object containing the matrices of coordinates of the ray
%          paths, nRx by 1.  Each matrix is nPts by 3
%    L: ray projection matrix, nRx by nSlowness
%
% -----------
%
% Bernard Giroux
% INRS
% 2013-01-23
%
