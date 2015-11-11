% function [tt, rays] = ttcr2dvti_sh(Vp0, Vs0, epsilon, delta, g, Tx, Rx, phase)
% TTCR2DVTI_PSV - P-SV waves travel times with curved rays in 2D in VTI anisotropic media
%
% Input
%    Vp0: velocity vector ( nCells by 1 )
%          values are vertical P-wave velocities
%          nCells is equal to g.nx*g.nz
%    Vs0: velocity vector ( nCells by 1 )
%          values are vertical S-wave velocities
%    epsilon: Thomsen anisotropy parameter vector ( nCells by 1 )
%    delta: Thomsen anisotropy parameter vector ( nCells by 1 )
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
%    phase: boolean equal to true for P-wave or false for SV-wave
%
%    *** nTx must be equal to nRx, i.e. each row define one Tx-Rx pair ***
%
% Output
%    tt: vector of traveltimes, nRx by 1
%    rays: cell object containing the matrices of coordinates of the ray
%          paths, nRx by 1.  Each matrix is nPts by 2
%
% -----------
%
% Bernard Giroux
% INRS-ETE
% 2013-03-06
%
