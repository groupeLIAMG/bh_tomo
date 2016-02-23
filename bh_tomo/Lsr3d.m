% [L,gridx,gridy,gridz] = Lsr3d(Tx, Rx, grx, gry, grz)
% LSR3D - Build projection matrix L for straight rays in 3D
%
% Input
%    Tx: transmitter coordinates, nTx by 3
%          1st column contains X coordinates, 2nd contains Y coordinates
%          3rd contains Z coordinates
%    Rx: receiver coordinates, nRx by 3
%          1st column contains X coordinates, 2nd contains Y coordinates
%          3rd contains Z coordinates
%    grx: vector of X coordinates of cell boundaries (sorted increasing)
%    gry: vector of Y coordinates of cell boundaries (sorted increasing)
%    grz: vector of Z coordinates of cell boundaries (sorted increasing)
%
%    *** nTx must be equal to nRx, i.e. each row define one Tx-Rx pair ***
%
% Output
%    L: ray projection matrix, 
%          nRx by ( (length(grx)-1)*(length(gry)-1)*(length(grz)-1) )
%    gridx: vector of X coordinates of the center of the cells
%    gridy: vector of Y coordinates of the center of the cells
%    gridz: vector of Z coordinates of the center of the cells
%
% ------------
%
% Bernard Giroux
% INRS
% 2013-01-25
%
