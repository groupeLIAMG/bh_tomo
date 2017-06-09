% function [tt, rays, L, corr_tt] = ttcr2d_ant(s, g, Tx, Rx, diam, theta, inWater, typeC)
% TTCR2D_ANT - Travel times with curved rays in 2D, with antenna length correction
%
% Input
%    s: slownes vector ( nSlowness by 1 )
%    g: structured variable defining the grid parameters
%          xmin: origin in X
%          zmin: origin in Z
%          dx: cell size in X
%          dz: cell size in Z
%          nx: number of cells in X
%          nz: number of cells in Z
%          nsx: number of secondary nodes per cell in X
%          nsz: number of secondary nodes per cell in Z
%          nsgx: number of subgrid cell per main grid cell in X
%          nsgz: number of subgrid cell per main grid cell in Z
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
%    *** nSlowness must equal g.nx*g.nz ***
%
% Output
%    tt: vector of traveltimes, nRx by 1
%    rays: cell object containing the matrices of coordinates of the ray
%          paths, nRx by 1.  Each matrix is nPts by 2
%    L: ray projection matrix, nRx by nSlowness
%    tt_corr: time correction
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
