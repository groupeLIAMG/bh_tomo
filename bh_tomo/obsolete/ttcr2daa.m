% function [tt, rays, L] = ttcr2daa(s, xi, theta, g, Tx, Rx)
% TTCR2DAA - Travel times with curved rays in 2D in elliptic anisotropic media
%
% Input
%    s: slowness vector ( nCells by 1 )
%          values are slownesses in X direction
%          nCells is equal to g.nx*g.nz
%    xi: anisotropy ratio vector ( nCells by 1 )
%          values are ratio of slowness in Z over slowness in X
%    theta: angle of rotation of the ellipse of anisotropy ( nCells by 1 ),
%          counter-clockwise from horizontal, units in radian
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
%
% -----------
%
% Bernard Giroux
% INRS-ETE
% 2010-05-03
%
