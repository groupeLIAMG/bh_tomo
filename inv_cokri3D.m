function tomo = inv_cokri3D(param, data, covar, grille, L, t_handle, g_handles)
% tomo = inv_cokri3D(param, data, covar, grille, L, t_handle, g_handles)
%
%

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

global kss ks0 Css
str = get_str_locale();
if ~isempty(t_handle)
	set(t_handle, 'String',[str.s225,' - ',str.s227])
else
	disp([str.s225,' - ',str.s227])
end
tomo.rays = {};
tomo.L = [];
tomo.invData = [];
if size(data,2)>=9
	tomo.no_trace = data(:,9);
end

if any(sum(data(:,[1 2 3])==data(:,[4 5 6]),2)==3)
    uiwait(errordlg('Coincident Tx & Rx positions, aborting inversion'))
    return
end

if isempty(L)
	[L,gridx,gridy,gridz] = Lsr3d(data(:,[1 2 3]),data(:,[4 5 6]),grille.grx,...
																grille.gry,grille.grz);
else
	%Calcul des vecteurs position des centres des blocs
	gridx = 0.5*(grille.grx(1:end-1)+grille.grx(2:end));
	gridy = 0.5*(grille.gry(1:end-1)+grille.gry(2:end));
	gridz = 0.5*(grille.grz(1:end-1)+grille.grz(2:end));
end

gp.xmin = grille.grx(1);
gp.ymin = grille.gry(1);
gp.zmin = grille.grz(1);
gp.dx = grille.grx(2)-grille.grx(1);
gp.dy = grille.gry(2)-grille.gry(1);
gp.dz = grille.grz(2)-grille.grz(1);
gp.nx = length(grille.grx)-1;
gp.ny = length(grille.gry)-1;
gp.nz = length(grille.grz)-1;

g3d = grid3d(gp);

cont = [];
if param.tomo_amp==0
	if ~isempty( grille.cont.slowness )
		cont = grille.cont.slowness.data;
	end
else
	if ~isempty( grille.cont.attenuation )
		cont = grille.cont.attenuation.data;
	end
end


np=length(L(1,:));
nt=length(L(:,1));
ind=1:np;

%Calcul de la matrice de covariance
if param.ival==2
	x=grille3d(min(gridx),gridx(2)-gridx(1),length(gridx),...
						 min(gridy),gridy(2)-gridy(1),length(gridy),...
						 min(gridz),gridz(2)-gridz(1),length(gridz));
	
	%	np=calcul_kss_global(np,x,ind,covar,gridx,gridz,m,n,t_handle);
	calcul_kss_global3D(x,ind,covar,t_handle);
	
	if ~isempty( cont ) && param.use_cont
		indc = zeros(size(cont,1),1);
		for i=1:size(cont,1)
			ind11=findnear(cont(i,1),x(:,1));
			ind22=findnear(cont(i,2),x(ind11,2));
			ind33=findnear(cont(i,3),x(ind11(ind22),3));
			indc(i)=min(ind11(ind22(ind33)));
		end
		clear ind11 ind22 ind33
		
		%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
		
		%Terme de gauche
		ks0=kss(indc,:);
		%Covariance des donnees de contrainte
		
		% contraintes dures
		Css = kss(indc,indc);
		
		% contraintes souples
		% la variance est contenue ds la 5e colonne
		if size(cont,2)==5
			Css = Css+diag(cont(:,5));
		end
	end
end

c0=[];
if length(data(1,:))>7 && covar.use_c0==1
	if data(:,8)~=0
		c0=data(:,8).^2;
	end
end


for noIter=1:param.nbreitrd+param.nbreitrc
	
	if ~isempty(t_handle)
		set(t_handle, 'String',[str.s225,' - ',str.s228,' ',num2str(noIter)])
		drawnow
	else
		disp([str.s225,' - ',str.s228,' ',num2str(noIter)]); drawnow
	end
	
	if noIter == 1
		l_moy = mean(data(:,7)./sum(L,2));
	else
		l_moy = mean(tomo.s);
	end
	mta = sum(L*l_moy,2);
	
	%vecteur des delta temps (t-moy(t))
	dt = data(:,7) - mta;
	
	mode = 1;
	if param.mode==2 && noIter==(param.nbreitrd+param.nbreitrc)
		mode = 2;
	end
	
	if ~isempty( cont ) && param.use_cont
		if mode==1 && ~isempty(t_handle)
			set(t_handle,'String',[str.s225,', it ',num2str(noIter),' - ', ...
										str.s232]);
			drawnow
		elseif mode==2 && ~isempty(t_handle)
			set(t_handle,'String',[str.s225,', it ',num2str(noIter),' - ', ...
										str.s233]);
			drawnow
		end
		scont = cont(:,4)-l_moy;
		[S_sim,ss] = simu_contrainte3d(L, dt, covar, indc, nt, ...
																	 gridx, gridy, gridz, param.nbresimu, ...
																	 scont, mode, c0);
	else
		if mode==1 && ~isempty(t_handle)
			set(t_handle,'String',[str.s225,', it ',num2str(noIter),' - ', ...
										str.s234]);
			drawnow
		elseif mode==2 && ~isempty(t_handle)
			set(t_handle,'String',[str.s225,', it ',num2str(noIter),' - ', ...
										str.s234b]);
			drawnow
		end
		[S_sim,ss] = cokri_simu3d(L, dt, covar, nt, gridx, gridy, gridz, ...
															param.nbresimu, mode, c0);
	end
	
	if param.tomo_amp==1
		ss(ss<-l_moy) = -l_moy;   % on fixe l'attenuation a zero
	end
	if mode==2
		[diff1,~,diff1_min]=choixsimu(L(:,ind),S_sim(ind,:),dt,c0);
		
		if param.use_cont==false || isempty( cont )
			sgr = deformationGraduelle(S_sim, L, c0, dt);
		else
			sgr = S_sim(:,diff1_min);
		end
		
		if  ~isempty(g_handles)
% 			set(g_handles{3},'DataAspectRatio',[1 1 1],'YDir','normal')
% 			title(str.s230,'FontSize',14)
% 			xlabel(str.s119)
%             ylabel(str.s119)
% 			zlabel(str.s120)
% 			if ~isempty(g_handles{1}), caxis(g_handles{3},g_handles{1}), end
% 			colorbar('peer',g_handles{3})
            ssplot = reshape((sgr+l_moy),gp.nx,gp.ny,gp.nz);
            [nx,ny,nz] = size(ssplot);
            ix = round(nx/2);
            iy = round(ny/2);
            iz = round(nz/2);
			if param.tomo_amp==0
                slice(1./ssplot,iy,ix,iz,'Parent',g_handles{4});
% 				imagesc(gridx,gridz,1./reshape((sgr+l_moy),length(gridz),...
% 																			 length(gridx)),'Parent',g_handles{4})
			else
				sgr(sgr<-l_moy) = -l_moy;
                ssplot = reshape((sgr+l_moy),gp.nx,gp.ny,gp.nz);
                slice(ssplot,iy,ix,iz,'Parent',g_handles{4});
% 				imagesc(gridx,gridz,reshape((sgr+l_moy),length(gridz),...
% 																		length(gridx)),'Parent',g_handles{4})
			end
			set(g_handles{4},'DataAspectRatio',[1 1 1],'YDir','normal')
			title(str.s231,'FontSize',12,'Parent',g_handles{4})
			xlabel(str.s119,'Parent',g_handles{4})
            ylabel(str.s119,'Parent',g_handles{4})
			zlabel(str.s120,'Parent',g_handles{4})
			if ~isempty(g_handles{1}), caxis(g_handles{4},g_handles{1}), end
			colorbar('peer', g_handles{4})
			eval(['colormap(',g_handles{2},')'])
			drawnow
		end
		tomo.s=ss+l_moy;
		tomo.simu=S_sim;%+l_moy;
		tomo.sgr=sgr;
		tomo.x=gridx;
		tomo.y=gridy;
		tomo.z=gridz;
		tomo.lmoy=l_moy;
		tomo.diff1=diff1;
		tomo.diff1_min=diff1_min;
		%tomo.dt=dt;
	else
		if  ~isempty(g_handles)
            ssplot = reshape((ss+l_moy),gp.nx,gp.ny,gp.nz);
            [nx,ny,nz] = size(ssplot);
            ix = round(nx/2);
            iy = round(ny/2);
            iz = round(nz/2);
			if param.tomo_amp==0
                slice(1./ssplot,iy,ix,iz,'Parent',g_handles{3});
            else
                slice(ssplot,iy,ix,iz,'Parent',g_handles{3});
			end
			set(g_handles{3},'DataAspectRatio',[1 1 1],'YDir','normal')
			title(str.s230,'FontSize',12,'Parent',g_handles{3})
			xlabel(str.s119,'Parent',g_handles{3})
			ylabel(str.s119,'Parent',g_handles{3})
            zlabel(str.s120,'Parent',g_handles{3})
			if ~isempty(g_handles{1}), caxis(g_handles{3},g_handles{1}), end
			colorbar('peer', g_handles{3})
			eval(['colormap(',g_handles{2},')'])
			drawnow
		end
		
		tomo.s=ss+l_moy;
		tomo.x=gridx;
		tomo.y=gridy;
		tomo.z=gridz;
		%tomo.dt=dt;
	end
	tomo.L = sparse(L);
	
	if param.tomo_amp==0 && noIter>=param.nbreitrd && param.nbreitrc > 0
		
		% update les rais
		s = tomo.s;
		% RADAR : si on a des points ou la vitesse > 0.2998 m/ns
		if param.radar == true
			s(s<3.3356) = 3.3356;
		end
		%        if param.use_cont == true
		%            s = appliqueCont(s, cont, gp);
		%        end
		if ~isempty(t_handle)
			set(t_handle, 'String',[str.s225,', it ',num2str(noIter),' - ',...
										str.s229])
			drawnow
		end
        [~, tomo.rays, L] = g3d.raytrace(s, data(:,[1 2 3]), data(:,[4 5 6]));
        tt = L*s;  % for some reason, L*s more accurate
	end
	
	if param.saveInvData == 1
		if noIter <= param.nbreitrd+param.nbreitrc || param.tomo_amp==1 || ...
					param.nbreitrc==0
			tt = tomo.L*tomo.s;
		end
		tomo.invData(noIter).res = data(:,7)-tt;
		tomo.invData(noIter).s = tomo.s;
	end
end
tomo.L = sparse(L);

if param.saveInvData == 1
	tomo.invData(end).date = datestr(now);
end
