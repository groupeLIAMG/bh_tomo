function D=deriv_second3d(M,N,P)
% D=deriv_second3d(M,N,O)
%
% M size of system in first dimension
% N size of system in second dimension
% N size of system in third dimension
  
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

D = zeros(M*N*P);

% nz = 7(M-2)(N-2)(P-2)
for m=2:M-1
	for n=2:N-1
		for p=2:P-1
			
			% derivative along X
			im = ((p-1)*N+n-1)*M+m-1;
			ic = ((p-1)*N+n-1)*M+m;
			ip = ((p-1)*N+n-1)*M+m+1;
			
			D(ic,im) = D(ic,im) + 1;
			D(ic,ic) = D(ic,ic) - 2;
			D(ic,ip) = D(ic,ip) + 1;
			
			% derivative along Y
			im = ((p-1)*N+n-2)*M+m;
			ip = ((p-1)*N+n  )*M+m;
			
			D(ic,im) = D(ic,im) + 1;
			D(ic,ic) = D(ic,ic) - 2;
			D(ic,ip) = D(ic,ip) + 1;
			
			% derivative along Z
			im = ((p-2)*N+n-1)*M+m;
			ip = ((p  )*N+n-1)*M+m;
			
			D(ic,im) = D(ic,im) + 1;
			D(ic,ic) = D(ic,ic) - 2;
			D(ic,ip) = D(ic,ip) + 1;
			
		end
	end
end


%
% m = 1
%
for n=1:N
	for p=1:P
		
		% along X
		
		ic = ((p-1)*N+n-1)*M+1;
		ip = ((p-1)*N+n-1)*M+2;
		
		D(ic,ic) = D(ic,ic) - 1;
		D(ic,ip) = D(ic,ip) + 1;
	end
end
for n=2:N-1
	for p=1:P
		
 		% along Y

		im = ((p-1)*N+n-2)*M+1;
		ic = ((p-1)*N+n-1)*M+1;
		ip = ((p-1)*N+n  )*M+1;

		D(ic,im) = D(ic,im) + 1;
		D(ic,ic) = D(ic,ic) - 2;
		D(ic,ip) = D(ic,ip) + 1;
	end
end
for n=1:N
	for p=2:P-1
		
 		% along Z

		im = ((p-2)*N+n-1)*M+1;
		ic = ((p-1)*N+n-1)*M+1;
		ip = ((p  )*N+n-1)*M+1;

		D(ic,im) = D(ic,im) + 1;
		D(ic,ic) = D(ic,ic) - 2;
		D(ic,ip) = D(ic,ip) + 1;
	end
end

%
% m = M
%
for n=1:N
	for p=1:P
		
		% along X
		
		im = ((p-1)*N+n-1)*M+M-1;
		ic = ((p-1)*N+n-1)*M+M;

		D(ic,im) = D(ic,im) + 1;
		D(ic,ic) = D(ic,ic) - 1;
	end
end
for n=2:N-1
	for p=1:P
		
 		% along Y

		im = ((p-1)*N+n-2)*M+M;
		ic = ((p-1)*N+n-1)*M+M;
		ip = ((p-1)*N+n  )*M+M;

		D(ic,im) = D(ic,im) + 1;
		D(ic,ic) = D(ic,ic) - 2;
		D(ic,ip) = D(ic,ip) + 1;
	end
end
for n=1:N
	for p=2:P-1
		
 		% along Z

		im = ((p-2)*N+n-1)*M+M;
		ic = ((p-1)*N+n-1)*M+M;
		ip = ((p  )*N+n-1)*M+M;

		D(ic,im) = D(ic,im) + 1;
		D(ic,ic) = D(ic,ic) - 2;
		D(ic,ip) = D(ic,ip) + 1;
	end
end


%
% n = 1
%
for m=2:M-1
	for p=1:P
		
		% along Y
		
		ic = ((p-1)*N  )*M+m;
		ip = ((p-1)*N+1)*M+m;
		D(ic,ic) = D(ic,ic) - 1;
		D(ic,ip) = D(ic,ip) + 1;
	end
end
for m=2:M-1
	for p=1:P

		% along X
		
		im = ((p-1)*N  )*M+m-1;
		ic = ((p-1)*N  )*M+m;
		ip = ((p-1)*N  )*M+m+1;

		D(ic,im) = D(ic,im) + 1;
		D(ic,ic) = D(ic,ic) - 2;
		D(ic,ip) = D(ic,ip) + 1;
	end
end
for m=2:M-1
	for p=2:P-1

		% along Z
		
		im = ((p-2)*N  )*M+m;
		ic = ((p-1)*N  )*M+m;
		ip = ((p  )*N  )*M+m;

		D(ic,im) = D(ic,im) + 1;
		D(ic,ic) = D(ic,ic) - 2;
		D(ic,ip) = D(ic,ip) + 1;
	end
end


%
% n = N
%
for m=2:M-1
	for p=1:P
		
		% along Y
		
		im = ((p-1)*N+N-2)*M+m;
		ic = ((p-1)*N+N-1)*M+m;
		D(ic,im) = D(ic,im) + 1;
		D(ic,ic) = D(ic,ic) - 1;
	end
end
for m=2:M-1
	for p=1:P

		% along X
		
		im = ((p-1)*N+N-1)*M+m-1;
		ic = ((p-1)*N+N-1)*M+m;
		ip = ((p-1)*N+N-1)*M+m+1;

		D(ic,im) = D(ic,im) + 1;
		D(ic,ic) = D(ic,ic) - 2;
		D(ic,ip) = D(ic,ip) + 1;
	end
end
for m=2:M-1
	for p=2:P-1

		% along Z
		
		im = ((p-2)*N+N-1)*M+m;
		ic = ((p-1)*N+N-1)*M+m;
		ip = ((p  )*N+N-1)*M+m;

		D(ic,im) = D(ic,im) + 1;
		D(ic,ic) = D(ic,ic) - 2;
		D(ic,ip) = D(ic,ip) + 1;
	end
end


%
% p = 1
%
for m=2:M-1
	for n=2:N-1

		% along Z
		
		ic = (  n-1)*M+m;
		ip = (N+n-1)*M+m;
		D(ic,ic) = D(ic,ic) - 1;
		D(ic,ip) = D(ic,ip) + 1;
	end
end
for m=2:M-1
	for n=2:N-1

		% along X
		
		im = (n-1)*M+m-1;
		ic = (n-1)*M+m;
		ip = (n-1)*M+m+1;

		D(ic,im) = D(ic,im) + 1;
		D(ic,ic) = D(ic,ic) - 2;
		D(ic,ip) = D(ic,ip) + 1;
	end
end
for m=2:M-1
	for n=2:N-1

		% along Y
		
		im = (n-2)*M+m;
		ic = (n-1)*M+m;
		ip = (n  )*M+m;

		D(ic,im) = D(ic,im) + 1;
		D(ic,ic) = D(ic,ic) - 2;
		D(ic,ip) = D(ic,ip) + 1;
	end
end


%
% p = P
%
for m=2:M-1
	for n=2:N-1

		% along Z
		
		im = ((P-2)*N+n-1)*M+m;
		ic = ((P-1)*N+n-1)*M+m;
		D(ic,im) = D(ic,im) + 1;
		D(ic,ic) = D(ic,ic) - 1;
	end
end
for m=2:M-1
	for n=2:N-1

		% along X
		
		im = ((P-1)*N+n-1)*M+m-1;
		ic = ((P-1)*N+n-1)*M+m;
		ip = ((P-1)*N+n-1)*M+m+1;

		D(ic,im) = D(ic,im) + 1;
		D(ic,ic) = D(ic,ic) - 2;
		D(ic,ip) = D(ic,ip) + 1;
	end
end
for m=2:M-1
	for n=2:N-1

		% along Y
		
		im = ((P-1)*N+n-2)*M+m;
		ic = ((P-1)*N+n-1)*M+m;
		ip = ((P-1)*N+n  )*M+m;

		D(ic,im) = D(ic,im) + 1;
		D(ic,ic) = D(ic,ic) - 2;
		D(ic,ip) = D(ic,ip) + 1;
	end
end
