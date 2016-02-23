% function [tt, rays, L] = ttcr2d(s, g, Tx, Rx)
% TTCR2D - Travel times with curved rays in 2D
%
% Input
%    s: slowness vector ( nSlowness by 1 )
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
%    *** nSlowness must equal g.nx*g.nz ***
%
% Output
%    tt: vector of traveltimes, nRx by 1
%    rays: cell object containing the matrices of coordinates of the ray
%          paths, nRx by 1.  Each matrix is nPts by 2
%    L: ray projection matrix, nRx by nSlowness
%
% -----------
%
% Bernard Giroux
% École Polytechnique de Montréal
% 2008-04-06
%

% Reference paper
%
% @article{gruber:1062,
%  author = {Thomas Gruber and Stewart A. Greenhalgh},
%  collaboration = {},
%  title = {Precision analysis of first-break times in grid models},
%  publisher = {SEG},
%  year = {1998},
%  journal = {Geophysics},
%  volume = {63},
%  number = {3},
%  pages = {1062-1065},
%  url = {http://link.aip.org/link/?GPY/63/1062/1},
%  doi = {10.1190/1.1444384}
% }
