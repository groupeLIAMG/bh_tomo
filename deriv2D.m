function [Cx,Cz]=deriv2D(pasx,pasz,flag)


nxx=numel(pasx); % nombre d'elements dans la direction x
nzz=numel(pasz); % nombre d'elements dans la direction z
N=nxx*nzz; % nombre total d'elements
%nx=(nxx-1)*nzz;
nz=nxx*(nzz-1);
Cx=sparse(N,N);





if flag ==1 % first derivative
   
   pasz=(pasz(1:end-1)+pasz(1:end-1))/2;
   pasx=(pasx(1:end-1)+pasx(1:end-1))/2; 
   pasx=1./pasx(:); % 
   pasz=1./pasz(:); 
   
   % in x direction 
   for i=1:nzz
   
       Cx1 = sparse(1+(i-1)*nxx:nxx-1+(i-1)*nxx,1+nxx*(i-1):nxx-1+nxx*(i-1),-pasx,N,N);
       Cx2 = sparse(1+(i-1)*nxx:nxx-1+(i-1)*nxx,2+nxx*(i-1):nxx+nxx*(i-1),pasx,N,N);
       Cx=Cx+Cx1+Cx2;

   end

   % in z direction
   pasz=pasz';

   pasz=repmat(pasz,nxx,1);
   pasz=pasz(:); % mise sous forme de vecteur 
% 
%    
   Cz=sparse(N,N); 
   Cz1 = sparse(1:(nzz-1)*nxx,1:(nzz-1)*nxx,-pasz,N,N);
   Cz2 = sparse(1:(nzz-1)*nxx,1+nxx:N,pasz,N,N);
   Cz=Cz+Cz1+Cz2;

elseif flag==2 % second derivative
    
   ax=zeros(nxx,1); bx=ax; az=zeros(nzz,1);bz=az;
   ax(1)=0;
   ax(end)=3/((pasx(end-1)+2*pasx(end))*pasx(end));
   bx(1)=3/((pasx(2)+2*pasx(1))*pasx(1));
   bx(end)=0;
   
   az(1)=0;
   az(end)=3./((pasz(end-1)+2*pasz(end))*pasz(end));
   
   bz(1)=3./((pasz(2)+2*pasz(1))*pasz(1));
   bz(end)=0;
   
   
   ax(2:end-1)=4./((pasx(3:end)+2*pasx(2:end-1)+pasx(1:end-2)).*(pasx(2:end-1)+pasx(1:end-2)));
   bx(2:end-1)=4./((pasx(3:end)+2*pasx(2:end-1)+pasx(1:end-2)).*(pasx(2:end-1)+pasx(3:end)));
   az(2:end-1)=4./((pasz(3:end)+2*pasz(2:end-1)+pasz(1:end-2)).*(pasz(2:end-1)+pasz(1:end-2)));
   bz(2:end-1)=4./((pasz(3:end)+2*pasz(2:end-1)+pasz(1:end-2)).*(pasz(2:end-1)+pasz(3:end)));

   ax=kron(ones(nzz,1),ax);
   bx=kron(ones(nzz,1),bx);
   az=kron(az,ones(nxx,1));
   bz=kron(bz,ones(nxx,1));
   
   
   cx=-(ax+bx);
   cz=-(az+bz);
   % reordering to use spdiags function
   ax=[ax(2:end);0];
   az=[az(nxx+1:end);zeros(nxx,1)];
   bx=[0;bx(1:end-1)]; bx(2)=0; % to avoid null space matrix
   bz=[zeros(nxx,1);bz(1:end-nxx)]; bz(2)=0; % to avoid null space matrix
   
   Cx = spdiags([ax cx bx],-1:1, N, N);
   Cz = spdiags([az cz bz],[-nxx 0 nxx], N, N);
   
end




    