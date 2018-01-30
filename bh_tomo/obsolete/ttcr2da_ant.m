% function [tt, rays, L] = ttcr2da(s, xi, g, Tx, Rx, diam, theta, inWater, typeC)
% TTCR2DA_ANT - Travel times with curved rays in 2D in elliptic anisotropic media, with antenna length correction
%
% Input
%    s: slowness vector ( nCells by 1 )
%          values are slownesses in X direction
%          nCells is equal to g.nx*g.nz
%    xi: anisotropy ratio vector ( nCells by 1 )
%          values are ratio of slowness in z over slowness in x
%    g: structured variable defining the grid parameters
%          xmin: origin in X
%          zmin: origin in Z
%          dx: cell size in X
%          dz: cell size in Z
%          nx: number of cells in X
%          nz: number of cells in Z
%          nsx: number of secondary nodes per cell in X
%          nsz: number of secondary nodes per cell in Z
%    Tx: transmitter coordinates, nTx by 2
%          1st column contains X coordinates, 2nd contains Z coordinates
%    Rx: receiver coordinates, nRx by 2
%          1st column contains X coordinates, 2nd contains Z coordinates
%    diam: hole diameter, nTx x 2
%          1st column for Tx, 2nd column for Rx
%    theta: antenna angle in radian, nTx x 2
%          1st column for Tx, 2nd column for Rx
%    inWater: antenna in water (true) or not (false), nTx x 2
%          1st column for Tx, 2nd column for Rx
%    typeC: type of antenna for the length correction
%          ( typeC=get(corr,'Type'); )
%
%    *** nTx must be equal to nRx, i.e. each row define one Tx-Rx pair ***
%
% Output
%    tt: vector of traveltimes, nRx by 1
%    rays: cell object containing the matrices of coordinates of the ray
%          paths, nRx by 1.  Each matrix is nPts by 2
%    L: ray projection matrix, nRx by 2*nCells
%          first nCells columns are length component in X direction
%          second nCells columns are length componnent in Z direction
%    tt_corr: time correction
%
% -----------
%
% Bernard Giroux
% INRS-ETE
% 2010-04-21
%