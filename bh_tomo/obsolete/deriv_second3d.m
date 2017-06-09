function D=deriv_second3d(M,N,P)
% D=deriv_second3d(M,N,P)
%
% M size of system in first dimension
% N size of system in second dimension
% P size of system in third dimension
  
% Copyright (C) 2013 Bernard Giroux
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

D = sparse(M*N*P,M*N*P);

% nz = 7(M-2)(N-2)(P-2)
mm=(2:M-1)';
nn=(2:N-1)';
pp=(2:P-1)';

m=kron(mm,ones(length(nn)*length(pp),1));
n=kron(kron(ones(length(mm),1),nn),ones(length(pp),1));
p=kron(ones(length(mm)*length(nn),1),pp);

% derivative along X
im = ((p-1)*N+n-1)*M+m-1;
ic = ((p-1)*N+n-1)*M+m;
ip = ((p-1)*N+n-1)*M+m+1;

D = D + sparse(ic,im, ones(size(ic)), M*N*P, M*N*P);
D = D + sparse(ic,ic, -2+zeros(size(ic)), M*N*P, M*N*P);
D = D + sparse(ic,ip, ones(size(ic)), M*N*P, M*N*P);

% derivative along Y
im = ((p-1)*N+n-2)*M+m;
ip = ((p-1)*N+n  )*M+m;

D = D + sparse(ic,im, ones(size(ic)), M*N*P, M*N*P);
D = D + sparse(ic,ic, -2+zeros(size(ic)), M*N*P, M*N*P);
D = D + sparse(ic,ip, ones(size(ic)), M*N*P, M*N*P);

% derivative along Z
im = ((p-2)*N+n-1)*M+m;
ip = ((p  )*N+n-1)*M+m;

D = D + sparse(ic,im, ones(size(ic)), M*N*P, M*N*P);
D = D + sparse(ic,ic, -2+zeros(size(ic)), M*N*P, M*N*P);
D = D + sparse(ic,ip, ones(size(ic)), M*N*P, M*N*P);

%
% m = 1
%
nn=(1:N)';
pp=(1:P)';

n=kron(nn,ones(length(pp),1));
p=kron(ones(length(nn),1),pp);

% along X
		
ic = ((p-1)*N+n-1)*M+1;
ip = ((p-1)*N+n-1)*M+2;

D = D + sparse(ic,ic, -ones(size(ic)), M*N*P, M*N*P);
D = D + sparse(ic,ip,  ones(size(ic)), M*N*P, M*N*P);

nn=(2:N-1)';
		
n=kron(nn,ones(length(pp),1));
p=kron(ones(length(nn),1),pp);

% along Y

im = ((p-1)*N+n-2)*M+1;
ic = ((p-1)*N+n-1)*M+1;
ip = ((p-1)*N+n  )*M+1;

D = D + sparse(ic,im, ones(size(ic)), M*N*P, M*N*P);
D = D + sparse(ic,ic, -2+zeros(size(ic)), M*N*P, M*N*P);
D = D + sparse(ic,ip, ones(size(ic)), M*N*P, M*N*P);

nn=(1:N)';
pp=(2:P-1)';
		
n=kron(nn,ones(length(pp),1));
p=kron(ones(length(nn),1),pp);

% along Z

im = ((p-2)*N+n-1)*M+1;
ic = ((p-1)*N+n-1)*M+1;
ip = ((p  )*N+n-1)*M+1;

D = D + sparse(ic,im, ones(size(ic)), M*N*P, M*N*P);
D = D + sparse(ic,ic, -2+zeros(size(ic)), M*N*P, M*N*P);
D = D + sparse(ic,ip, ones(size(ic)), M*N*P, M*N*P);

%
% m = M
%
nn=(1:N)';
pp=(1:P)';

n=kron(nn,ones(length(pp),1));
p=kron(ones(length(nn),1),pp);
		
% along X

im = ((p-1)*N+n-1)*M+M-1;
ic = ((p-1)*N+n-1)*M+M;

D = D + sparse(ic,im,  ones(size(ic)), M*N*P, M*N*P);
D = D + sparse(ic,ic, -ones(size(ic)), M*N*P, M*N*P);

nn=(2:N-1)';

n=kron(nn,ones(length(pp),1));
p=kron(ones(length(nn),1),pp);

% along Y

im = ((p-1)*N+n-2)*M+M;
ic = ((p-1)*N+n-1)*M+M;
ip = ((p-1)*N+n  )*M+M;

D = D + sparse(ic,im, ones(size(ic)), M*N*P, M*N*P);
D = D + sparse(ic,ic, -2+zeros(size(ic)), M*N*P, M*N*P);
D = D + sparse(ic,ip, ones(size(ic)), M*N*P, M*N*P);

nn=(1:N)';
pp=(2:P-1)';

n=kron(nn,ones(length(pp),1));
p=kron(ones(length(nn),1),pp);

% along Z

im = ((p-2)*N+n-1)*M+M;
ic = ((p-1)*N+n-1)*M+M;
ip = ((p  )*N+n-1)*M+M;

D = D + sparse(ic,im, ones(size(ic)), M*N*P, M*N*P);
D = D + sparse(ic,ic, -2+zeros(size(ic)), M*N*P, M*N*P);
D = D + sparse(ic,ip, ones(size(ic)), M*N*P, M*N*P);

%
% n = 1
%
mm=(2:M-1)';
pp=(1:P)';

m=kron(mm,ones(length(pp),1));
p=kron(ones(length(mm),1),pp);

% along Y
		
ic = ((p-1)*N  )*M+m;
ip = ((p-1)*N+1)*M+m;

D = D + sparse(ic,ic, -ones(size(ic)), M*N*P, M*N*P);
D = D + sparse(ic,ip,  ones(size(ic)), M*N*P, M*N*P);

% along X
		
im = ((p-1)*N  )*M+m-1;
ic = ((p-1)*N  )*M+m;
ip = ((p-1)*N  )*M+m+1;

D = D + sparse(ic,im, ones(size(ic)), M*N*P, M*N*P);
D = D + sparse(ic,ic, -2+zeros(size(ic)), M*N*P, M*N*P);
D = D + sparse(ic,ip, ones(size(ic)), M*N*P, M*N*P);

mm=(2:M-1)';
pp=(2:P-1)';

m=kron(mm,ones(length(pp),1));
p=kron(ones(length(mm),1),pp);

% along Z

im = ((p-2)*N  )*M+m;
ic = ((p-1)*N  )*M+m;
ip = ((p  )*N  )*M+m;

D = D + sparse(ic,im, ones(size(ic)), M*N*P, M*N*P);
D = D + sparse(ic,ic, -2+zeros(size(ic)), M*N*P, M*N*P);
D = D + sparse(ic,ip, ones(size(ic)), M*N*P, M*N*P);


%
% n = N
%
pp=(1:P)';

m=kron(mm,ones(length(pp),1));
p=kron(ones(length(mm),1),pp);

% along Y
		
im = ((p-1)*N+N-2)*M+m;
ic = ((p-1)*N+N-1)*M+m;

D = D + sparse(ic,im,  ones(size(ic)), M*N*P, M*N*P);
D = D + sparse(ic,ic, -ones(size(ic)), M*N*P, M*N*P);

% along X

im = ((p-1)*N+N-1)*M+m-1;
ic = ((p-1)*N+N-1)*M+m;
ip = ((p-1)*N+N-1)*M+m+1;

D = D + sparse(ic,im, ones(size(ic)), M*N*P, M*N*P);
D = D + sparse(ic,ic, -2+zeros(size(ic)), M*N*P, M*N*P);
D = D + sparse(ic,ip, ones(size(ic)), M*N*P, M*N*P);

pp=(2:P-1)';

m=kron(mm,ones(length(pp),1));
p=kron(ones(length(mm),1),pp);

% along Z
		
im = ((p-2)*N+N-1)*M+m;
ic = ((p-1)*N+N-1)*M+m;
ip = ((p  )*N+N-1)*M+m;

D = D + sparse(ic,im, ones(size(ic)), M*N*P, M*N*P);
D = D + sparse(ic,ic, -2+zeros(size(ic)), M*N*P, M*N*P);
D = D + sparse(ic,ip, ones(size(ic)), M*N*P, M*N*P);


%
% p = 1
%
nn=(2:N-1)';

m=kron(mm,ones(length(nn),1));
n=kron(ones(length(mm),1),nn);

% along Z
		
ic = (  n-1)*M+m;
ip = (N+n-1)*M+m;

D = D + sparse(ic,ic, -ones(size(ic)), M*N*P, M*N*P);
D = D + sparse(ic,ip,  ones(size(ic)), M*N*P, M*N*P);

% along X
		
im = (n-1)*M+m-1;
ic = (n-1)*M+m;
ip = (n-1)*M+m+1;

D = D + sparse(ic,im, ones(size(ic)), M*N*P, M*N*P);
D = D + sparse(ic,ic, -2+zeros(size(ic)), M*N*P, M*N*P);
D = D + sparse(ic,ip, ones(size(ic)), M*N*P, M*N*P);

% along Y
		
im = (n-2)*M+m;
ic = (n-1)*M+m;
ip = (n  )*M+m;

D = D + sparse(ic,im, ones(size(ic)), M*N*P, M*N*P);
D = D + sparse(ic,ic, -2+zeros(size(ic)), M*N*P, M*N*P);
D = D + sparse(ic,ip, ones(size(ic)), M*N*P, M*N*P);


%
% p = P
%

% along Z
		
im = ((P-2)*N+n-1)*M+m;
ic = ((P-1)*N+n-1)*M+m;

D = D + sparse(ic,im,  ones(size(ic)), M*N*P, M*N*P);
D = D + sparse(ic,ic, -ones(size(ic)), M*N*P, M*N*P);

% along X

im = ((P-1)*N+n-1)*M+m-1;
ic = ((P-1)*N+n-1)*M+m;
ip = ((P-1)*N+n-1)*M+m+1;

D = D + sparse(ic,im, ones(size(ic)), M*N*P, M*N*P);
D = D + sparse(ic,ic, -2+zeros(size(ic)), M*N*P, M*N*P);
D = D + sparse(ic,ip, ones(size(ic)), M*N*P, M*N*P);

% along Y
		
im = ((P-1)*N+n-2)*M+m;
ic = ((P-1)*N+n-1)*M+m;
ip = ((P-1)*N+n  )*M+m;

D = D + sparse(ic,im, ones(size(ic)), M*N*P, M*N*P);
D = D + sparse(ic,ic, -2+zeros(size(ic)), M*N*P, M*N*P);
D = D + sparse(ic,ip, ones(size(ic)), M*N*P, M*N*P);
