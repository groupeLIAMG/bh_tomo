Tx = [1 1; 1 1; 1 1; 1.0 1; 1 1; 1 5; 1 5; 1 5]; 
Rx = [5 1; 5 3; 5 5; 1.6 5; 1 5; 5 1; 5 3; 5 5];
Tx2 = [Tx(:,1) zeros(8,1) Tx(:,2)]
Rx2 = [Rx(:,1) zeros(8,1) Rx(:,2)]

%Tx = [1 1; 1 5];
%Rx = [5 3; 5 3];
%Tx2 = [1 0 1; 1 0 5];
%Rx2 = [5 0 3; 5 0 3];

grx = -0.5:1:5.5;
grz = -0.5:1:6.5;

[L,gridx,gridz]=Lsr2d(Tx,Rx,grx,grz);

[L2,gridx2,gridz2]=rais_droits(Tx2,Rx2,grx,grz,1,3);

L2 = sparse(L2);

d = L-L2;

imagesc(d)
colorbar

pause


load ~/Projets/BHRS/toto

grx = grille.grx';
grz = grille.grz';

Tx = P1(:,[1 3]);
Rx = P2(:,[1 3]);

P1(:,2)=0;
P2(:,2)=0;

tic
[L,gridx,gridz]=Lsr2d(Tx,Rx,grx,grz);
toc

tic
[L2,gridx2,gridz2]=rais_droits(P1,P2,grx,grz,1,3,'y');
toc
L2 = sparse(L2);

d = L-L2;

imagesc(d)
colorbar
set(gca,'DataAspectRatio',[1 1 1])

max(d(:))

