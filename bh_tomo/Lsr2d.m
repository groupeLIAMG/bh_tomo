% [L,gridx,gridy] = Lsr2d(Tx, Rx, grx, grz)
% LSR2D - Build projection matrix L for straight rays in 2D
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
%    L: ray projection matrix, nRx by ( (length(grx)-1)*(length(grz)-1) )
%    gridx: vector of X coordinates of the center of the cells
%    gridz: vector of Z coordinates of the center of the cells
%
% ------------
%
% Bernard Giroux
% École Polytechnique de Montréal
% 2008-05-21
%
