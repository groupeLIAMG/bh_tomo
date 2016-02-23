% [L,gridx,gridy] = Lsr2da(Tx, Rx, grx, grz)
% LSR2DA - Build projection matrix L for straight rays in 2D elliptical
% anisotripic (VTI) media
%
% Input
%    Tx: transmitter coordinates, nTx by 2
%          1st column contains X coordinates, 2nd contains Z coordinates
%    Rx: receiver coordinates, nRx by 2
%          1st column contains X coordinates, 2nd contains Z coordinates
%    grx: vector of X coordinates of cell boundaries (sorted increasing)
%    grz: vector of Z coordinates of cell boundaries (sorted increasing)
%
%    *** nTx must be equal to nRx, i.e. each row define one Tx-Rx pair ***
%
% Output
%    L: ray projection matrix, nRx by 2x( (length(grx)-1)*(length(grz)-1) )
%       L = [Lx Lz], where Lx contains segments length along X
%                          Lz contains segments length along Z
%    gridx: vector of X coordinates of the center of the cells
%    gridz: vector of Z coordinates of the center of the cells
%
% ------------
%
% Bernard Giroux
% INRS-ETE
% 2009-12-18
%
