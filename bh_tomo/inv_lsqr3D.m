function  tomo =  inv_lsqr3D(param,  data, grille,  L, t_handle,  g_handles)
% tomo = inv_lsqr3D(param, data, grille, L, t_handle, g_handles)
%
% input
%   param : struct variable with the following members
%       use_cont:  use constraints (1) ot not (0)
%       tomo_amp:  invert traveltimes (0) or amplitudes (1)
%       nbreitrd:  number of straight ray iterations to perform
%       nbreitrc:  number of curved ray iterations to perform
%       alpha:     Lagrangian of second derivative smoothing
%       saveInvData: save intermediate inversion results (not yet
%                     implemented, set to 0)
%       dv_max:    max velocity variation per iteration (0<dv_max<1)
%       radar:     1 if we invert gpr data, 0 for seismic/acoustic data
%       nbreiter:  LSQR -> max number of iterations
%       tol:       LSQR -> tolerance on objective function
%       gradmin:   LSQR -> min gradient on objective function
%       correlmin: LSQR -> min correlation from one iteration to the next
%
%   data: matrix nDataPts x 9, with the following columns
%       Tx_x Tx_y Tx_z Rx_x Rx_y Rx_z travel_time trav_time_variance ray_number
%
%   grille: struct variable of the grid, with the following members
%       grx:     X coordinates of grid cell edges
%       gry:     Y coordinates of grid cell edges
%       grz:     Z coordinates of grid cell edges
%       cont: struct variable for the constraints, with the following  members
%           slowness:    matrix nConstrPts x 5 with the following columns
%                        x_coord y_coord z_coord slowness variance_of_slowness
%           attenuation: matrix nConstrPts x 5 with the following columns
%                        x_coord y_coord z_coord attenuation variance_of_attenuation
%
%   L: precomputed ray matrix (set to [] if not readily available)
%
%   t_handle: text handle used in GUI, set to [] if using command line
%   g_handle: graphics handle used in GUI, set to [] if using command line
%
%
% output
%   tomo : struct variable with the following members
%       s:    slowness vector
%       x:    x coord of slowness pts
%       y:    y coord of slowness pts
%       z:    z coord of slowness pts
%       res : residuals
%       L:    ray matrix
%       rais: array of cells containing ray paths
%
%
%

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

str = get_str_locale();
if ~isempty(t_handle)
	set(t_handle, 'String',[str.s226,' - ',str.s227])
else
	disp([str.s226,' - ',str.s227])
end

tomo.rays = {};
tomo.L = [];
corr_tt_Lant = [];
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
x=grille3d(min(gridx),gridx(2)-gridx(1),length(gridx),...
    min(gridy),gridy(2)-gridy(1),length(gridy),...
    min(gridz),gridz(2)-gridz(1),length(gridz));

gp.xmin = grille.grx(1);
gp.ymin = grille.gry(1);
gp.zmin = grille.grz(1);
gp.dx = grille.grx(2)-grille.grx(1);
gp.dy = grille.gry(2)-grille.gry(1);
gp.dz = grille.grz(2)-grille.grz(1);
gp.nx = length(grille.grx)-1;
gp.ny = length(grille.gry)-1;
gp.nz = length(grille.grz)-1;

nThreads = feature('numcores');
g3d = grid3d(gp, nThreads);

cont = [];
if param.use_cont == true
    if param.tomo_amp==0
        if ~isempty( grille.cont.slowness )
            cont = grille.cont.slowness.data;
        end
    else
        if ~isempty( grille.cont.attenuation )
            cont = grille.cont.attenuation.data;
        end
	end
end

for noIter=1:param.nbreitrc+param.nbreitrd
	if ~isempty(t_handle)
		set(t_handle, 'String',[str.s226,' - ',str.s228,' ', ...
			num2str(noIter)])
	else
		disp([str.s226,' - ',str.s228,' ', num2str(noIter)])
	end
	np=length(L(1,:));
	%nt=length(L(:,1));

	%Ne pas tenir compte des pixels non visites
	% 	if strcmp(param.pixnonvisite,'oui')==1
	% 		mask=sum(L);
	% 		for i=1:np
	% 			mask(i)=isempty(find(mask(i)==0));
	% 		end
	% 		ind=find(mask==1);
	% 		ind2=find(mask==0);
	% 		mask(ind2)=nan;
	% 	else
	ind=1:np;
	ind2=[];
	% 	end

	%Calcul de la moyenne
	if noIter == 1
		l_moy=mean(data(:,7)./sum(L,2));
	else
		l_moy = mean(tomo.s);
	end
	mta=sum(L*l_moy,2);

	if noIter==1
		s_o = l_moy*ones(np,1);
	end

	%vecteur des delta temps (t-moy(t))
	if isempty(corr_tt_Lant)
		dt = data(:,7) - mta;
	else
		dt = data(:,7) - mta + corr_tt_Lant;
	end

	%%%%%%%%%%%%%%%%%%%%%%%%%%%
	%longueur du vecteur solution
	ss=zeros(np,1);
    p=length(gridz);
    n=length(gridy);
    m=length(gridx);
	D=deriv_second3d(m,n,p);  
	D=D(ind,ind);
	ld=[L(:,ind);D*param.alpha];
	dls=[dt;zeros(length(ind),1)];
	[xxx,res_m,var_res]=lsqr_erwan2(ld,dls,param.nbreiter,param.tol,...
		param.gradmin,param.correlmin,cont,x,l_moy,ind);

	if max(abs(s_o./(xxx+l_moy) - 1))>param.dv_max
		fac = min(abs( (s_o/(param.dv_max+1)-l_moy)./xxx ));
		xxx = fac*xxx;
		s_o = xxx+l_moy;
	end

	ss(ind)=xxx;
	ss(ind2)=nan;

	if ~isempty(g_handles)
        ssplot = reshape((ss+l_moy),gp.nx,gp.ny,gp.nz);
        [nx,ny,nz] = size(ssplot);
        ix = gridx(round(nx/2));
        iy = gridy(round(ny/2));
        iz = gridz(round(nz/2));
		if param.tomo_amp==0
            slice(gridx,gridy,gridz,1./ssplot,ix,iy,iz,'Parent',g_handles{3});
        else
            slice(gridx,gridy,gridz,ssplot,ix,iy,iz,'Parent',g_handles{3});
		end
		set(g_handles{3},'DataAspectRatio',[1 1 1])
		title('LSQR','Parent',g_handles{3})
        xlabel(['X ',str.s119],'Parent',g_handles{3})
        ylabel(['Y ',str.s119],'Parent',g_handles{3})
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
	tomo.res=res_m;
	tomo.var_res=var_res;
	tomo.corr_tt_Lant = corr_tt_Lant;

	if param.tomo_amp==0 && noIter>=param.nbreitrd && ...
			param.nbreitrc > 0
		% update les rais
		s = tomo.s;
		% RADAR : si on a des points ou la vitesse > 0.2998 m/ns
		if param.radar == true
			s(s<3.3356) = 3.3356;
        end

        %        if param.use_cont == true
        %            s = appliqueCont(s, cont, gp);
        %        end
        [~, tomo.rays, L] = g3d.raytrace(s, data(:,[1 2 3]), data(:,[4 5 6]));
        tt = L*s;  % for some reason, L*s more accurate
        
	end
	if param.saveInvData == 1
		if noIter < param.nbreitrd || param.tomo_amp==1 || ...
				param.nbreitrc==0
			tt = L*tomo.s;
		end
		tomo.invData(noIter).res = data(:,7)-tt;
		tomo.invData(noIter).s = tomo.s;
	end
end


if param.saveInvData == 1
    tomo.invData(end).date = datestr(now);
end

