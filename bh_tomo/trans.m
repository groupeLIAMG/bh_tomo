function [cx,rot]=trans(cx,model,im)
% [cx,rot]=trans(cx,model,im)
%
% TRANS is called from COKRI2. It takes as input original coordinates and
%       return the rotated and reduced coordinates following specifications
%       described in model(im,:)
%
% Rotations are all performed anticlockwise with the observer located on
% the positive side of the axis and looking toward the origin. In 3D,
% rotations are performed first along z, then along rotated y and then
% along twice rotated x.
% Author: D. Marcotte
% Version 2.1  97/aug/18

% Copyright (C) 2005 Denis Marcotte
%
% This program is free software; you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation; either version 3 of the License, or
% (at your option) any later version.
%
% This program is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details.
% 
% You should have received a copy of the GNU General Public License
% along with this program.  If not, see <http://www.gnu.org/licenses/>.
% 
% 

  
  % some constants are defined

[n,d]=size(cx);
[m,p]=size(model);

% check for 1-D or isotropic model

if p-1>d,

  % perform rotation counterclockwise

  if d==2,
     ang=model(im,4); cang=cos(ang/180*pi); sang=sin(ang/180*pi);
     rot=[cang,-sang;sang,cang];
  else

     % rotation matrix in 3-D is computed around z, y and x in that order

     angz=model(im,7); cangz=cos(angz/180*pi); sangz=sin(angz/180*pi);
     angy=model(im,6); cangy=cos(angy/180*pi); sangy=sin(angy/180*pi);
     angx=model(im,5); cangx=cos(angx/180*pi); sangx=sin(angx/180*pi);
     rotz=[cangz,-sangz,0;sangz,cangz,0;0 0 1];
     roty=[cangy,0,sangy;0 1 0;-sangy,0,cangy];
     rotx=[1 0 0;0 cangx -sangx;0 sangx cangx];
     rot=rotz*roty*rotx;
  end

  % rotation is performed around z, y and x in that order, the other coordinates are left unchanged.

  dm=min(3,d);
  cx(:,1:dm)=cx(:,1:dm)*rot;
  t=[model(im,2:1+dm),ones(d-dm,1)];
  t=diag(t);
else
  t=eye(d)*model(im,2);
end

% perform contractions or dilatations (reduced h)

 cx=cx/t;

 
