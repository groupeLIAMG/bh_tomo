function [L,gridx,gridy]=rais_droits(P1,P2,grx,gry,Dir1,Dir2,varargin)
% [L,gridx,gridy]=rais_droits(P1,P2,grx,gry,Dir1,Dir2)

% Copyright (C) 2005 Erwan Gloaguen, Bernard Giroux
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

%format long
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%5  
grx=grx(:);
gry=gry(:);


%**************************************************

if (Dir1==1 && Dir2==2)
	Dir3=3;
elseif (Dir1==2 && Dir2==3)
	Dir3=1;
else
	Dir3=2;
end

if nargin>=7
	choix=varargin{1};
else
	choix='n';
end
if nargin>=8
	barre=varargin{1};
else
	barre=false;
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%5
%Regarde et transforme le fichier
% position de telle sorte que la
% coordonnée x la plus petite soit
% dans le forage P1
indo=find(P1(:,1)>P2(:,1));

if ~isempty(indo)
	temp=P1;
	P1(indo,1:3)=P2(indo,1:3);
	P2(indo,1:3)=temp(indo,1:3);
	clear temp indo
end

str = get_str_locale();
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if strcmp(choix,'n')==1
	h=figure('numbertitle','off','name',str.s64);
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
while strcmp(choix,'n')==1
	set(h,'Visible','on')
	clf(h)
	hold on
	plot([P1(:,Dir1) P2(:,Dir1)]',[P1(:,Dir2),P2(:,Dir2)]','r');
	aa=[min(grx)*ones(length(gry),1) max(grx)*ones(length(gry),1);grx grx];
	bb=[gry gry; min(gry)*ones(length(grx),1) max(gry)*ones(length(grx),1)];
	plot(aa',bb','Color',[0.5 0.5 0.5])
	% set(gca,'DataAspectRatio',[1 1 1]);
	set(gca,'Ydir','reverse');
	axis([min(grx) max(grx) min(gry) max(gry)])
	hold off
	choix = 'y';
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%55
Xm=grx;
Ym=gry;

%Calcul des vecteurs position des centres des blocs
gridx = 0.5*(grx(1:end-1)+grx(2:end));
gridy = 0.5*(gry(1:end-1)+gry(2:end));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%5
%Définition des longueurs des vecteurs
l1=length(Xm)-1;
l2=length(Ym)-1;

nbrais=length(P1(:,1));
%minz=NaN.*zeros(l1*l2,1);
%maxz=NaN.*zeros(l1*l2,1);
L=zeros(nbrais,l1*l2);
dist = zeros(nbrais,1);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Début de l'itération par tir

if barre, h = waitbar(0,str.s259); end
for j=1:nbrais
	
	if barre, waitbar(j/nbrais,h), end
	
	Xs=P1(j,Dir1);
	Ys=P1(j,Dir2);
	Zs=P1(j,Dir3);
	Xr=P2(j,Dir1);
	Yr=P2(j,Dir2);
	Zr=P2(j,Dir3);
	%*******************************************************
	%Calcul des distances et des paramètres des pentes
	%et ordonnées des rais
	
	dist(j)=sqrt((Xs-Xr).^2+(Ys-Yr).^2+(Zs-Zr).^2);
	
	if Xs ~= Xr % le rai n'est pas vertical
		a=(Yr-Ys)/(Xr-Xs);
		b=Ys-a*Xs;
	
		%Si émetteur au même niveau que récepteur
		if Ys==Yr
			is=find(Xs>=Xm, 1, 'last' );
			js=find(Ys>Ym, 1, 'last' );
			ir=find(Xr>Xm, 1, 'last' );
			%jr=find(Yr>Ym, 1, 'last' );
			%disp('Ys = Yr')
			
%			DX=[(Ym(js+1:jr)-b)./a Ym(js+1:jr) ;Xm(is+1:ir) a.*Xm(is+1:ir)+b];
			DX=[Xm(is+1:ir) a.*Xm(is+1:ir)+b];
			DX=sort(DX);
			DX=[Xs Ys;DX;Xr Yr];
			%Calcul des longueurs de rais dans les cellules
			ind=1:length(DX)-1;
			ind2=2:length(DX);
%			dii=((DX(ind2,2)-DX(ind,2)).^2+(DX(ind2,1)-DX(ind,1)).^2).^0.5;
			dii=DX(ind2,1)-DX(ind,1);
			%Trouvers les rais interceptant noeuds de la grille
			ind0=find(dii>1e-10);
			dii=dii(ind0);
			ind0=[ind0;length(DX(:,1))];
			DX=DX(ind0,:);
			%Calcul des indices dans la grille
			if isempty(js)
				js=1;
			end
			indx=zeros(length(DX(:,1))-1,1);
			indy=ones(length(DX(:,1))-1,1)*js;
			indx(1)=is;
			for ii=2:length(DX(:,1))-1
				indx(ii)=isempty(find(Xm==DX(ii,1), 1));
				indx(ii)=indx(ii-1)+abs(indx(ii)-1);
			end
		else
			%Définition de la droite y=ax+b
			
			is=find(Xs>=Xm, 1, 'last' );
			js=find(Ys>=Ym, 1, 'last' );
			ir=find(Xr>Xm, 1, 'last' );
			jr=find(Yr>Ym, 1, 'last' );
			
			if Ys < Yr
				%disp('Ys < Yr')
				
				DX=[(Ym(js+1:jr)-b)./a Ym(js+1:jr) ;Xm(is+1:ir) a.*Xm(is+1:ir)+b];
				DX=sort(DX);
				DX=[Xs Ys;DX;Xr Yr];
				%Calcul des longueurs de rais dans les cellules
				ind=1:length(DX)-1;
				ind2=2:length(DX);
				dii=((DX(ind2,2)-DX(ind,2)).^2+(DX(ind2,1)-DX(ind,1)).^2).^0.5;
				%Trouvers les rais interceptant noeuds de la grille
				ind0=find(dii>1e-10);
				dii=dii(ind0);
				ind0=[ind0;length(DX(:,1))];
				DX=DX(ind0,:);
				%Calcul des indices dans la grille
				indx=zeros(length(DX(:,1))-1,1);
				indy=zeros(length(DX(:,1))-1,1);
				indx(1)=is;indy(1)=js;
				for ii=2:length(DX(:,1))-1
					indx(ii)=isempty(find(Xm==DX(ii,1), 1));
					indy(ii)=isempty(find(Ym==DX(ii,2), 1));
					indx(ii)=indx(ii-1)+abs(indx(ii)-1);
					indy(ii)=indy(ii-1)+abs(indy(ii)-1);
				end
				
			end
			
			%Si la source est au dessus du récepteur
			if Ys > Yr
				%disp('Ys > Yr')
				jr=find(Ys>Ym, 1, 'last' );
				js=find(Yr>=Ym, 1, 'last' );
				
				DX=[(Ym(js+1:jr)-b)./a Ym(js+1:jr) ;Xm(is+1:ir) a.*Xm(is+1:ir)+b];
				
				DX=[Xs Ys;DX;Xr Yr];
				DX=sort(DX);
				DX(:,2)=flipud(DX(:,2));
				%Calcul des longueurs de rais dans les cellules
				ind=1:length(DX)-1;
				ind2=2:length(DX);
				dii=((DX(ind2,2)-DX(ind,2)).^2+(DX(ind2,1)-DX(ind,1)).^2).^0.5;
				%Trouver les rais interceptant noeuds de la grille
				ind0=find(dii>1e-10);
				dii=dii(ind0);
				ind0=[ind0;length(DX(:,1))];
				DX=DX(ind0,:);
				%Calcul des indices dans la grille
				indx=zeros(length(DX(:,1))-1,1);
				indy=zeros(length(DX(:,1))-1,1);
				indx(1)=is;indy(1)=jr;
				for ii=2:length(DX(:,1))-1
					indx(ii)=isempty(find(Xm==DX(ii,1), 1));
					indy(ii)=isempty(find(Ym==DX(ii,2), 1));
					indx(ii)=indx(ii-1)+abs(indx(ii)-1);
					indy(ii)=indy(ii-1)-abs(indy(ii)-1);
				end
				%	indy=flipud(indy);
			end
		end
		
	else  % on a un rai vertical (ajoutÃ© BG)
		ymin = min([Ys Yr]);
		ymax = max([Ys Yr]);
		indy = find( Ym > ymin & Ym < ymax );
		Y = [ymin; Ym(indy); ymax];
		dii = diff(Y);
		indy = [indy(1)-1; indy];
		indx = zeros(size(indy)) + find(Xs>=Xm, 1, 'last' );
	end
		
	for k=1:length(indx)
%		a(k)=l2*(indx(k)-1)+indy(k);
		L(j,l2*(indx(k)-1)+indy(k))=dii(k);
	end
end
if barre, close(h); end


%vérification si les distances sont bien calculées
%verif=sum(L,2)-dist;
%trc=find(verif>1e-10);
%if isempty(trc)==0
%  disp('Y''a une couille!')
%end
