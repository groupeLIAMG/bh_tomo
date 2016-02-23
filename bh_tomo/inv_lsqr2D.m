function  tomo =  inv_lsqr2D(param,  data, grille,  L, t_handle,  g_handles)
% tomo = inv_lsqr2D(param, data, grille, L, t_handle, g_handles)
%
% input
%   param : struct variable with the following members
%       use_cont:  use constraints (1) ot not (0)
%       tomo_amp:  invert traveltimes (0) or amplitudes (1)
%       nbreitrd:  number of straight ray iterations to perform
%       nbreitrc:  number of curved ray iterations to perform
%       nbreitra:  number of iterations with antenna correction
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
%         *** Important ***
%         Tx and Rx should be in the (x,z) plane.  If not, swap one pair
%         of coordinate axes.
%
%   grille: struct variable of the grid, with the following members
%       grx:     X coordinates of grid cell edges
%       grz:     Z coordinates of grid cell edges
%       cont: struct variable for the constraints, with the following  members
%           slowness:    matrix nConstrPts x 4 with the following columns
%                        z_coord x_coord slowness variance_of_slowness
%           attenuation: matrix nConstrPts x 4 with the following columns
%                        z_coord x_coord attenuation variance_of_attenuation
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
%       z:    z coord of slowness pts
%       res : residuals
%       L:    ray matrix
%       rais: array of cells containing ray paths
%
%  to view the results
%     imagesc(tomo.x, tomo.z, reshape(tomo.s,length(tomo.z),length(tomo.x)))
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

%Dir1 = 1;
%Dir2 = 3;
tomo.rays = {};
tomo.L = [];
corr_tt_Lant = [];
if size(data,2)>=9
	tomo.no_trace = data(:,9);
end

if any(sum(data(:,[1 3])==data(:,[4 6]),2)==2)
    uiwait(errordlg('Coincident Tx & Rx positions, aborting inversion'))
    return
end

if isempty(L)
	%	[L,gridx,gridz]=rais_droits(data(:,1:3),data(:,4:6),grille.grx,grille.grz,Dir1,Dir2,'y');
	[L,gridx,gridz] = Lsr2d(data(:,[1 3]),data(:,[4 6]),grille.grx,...
													grille.grz);
else
	%Calcul des vecteurs position des centres des blocs
	gridx = 0.5*(grille.grx(1:end-1)+grille.grx(2:end));
	gridz = 0.5*(grille.grz(1:end-1)+grille.grz(2:end));
end
n=length(gridz);
m=length(gridx);
x=grille2d(min(gridz),gridz(2)-gridz(1),n,min(gridx),gridx(2)-gridx(1),m);

gp.xmin = grille.grx(1);
gp.zmin = grille.grz(1);
gp.dx = grille.grx(2)-grille.grx(1);
gp.dz = grille.grz(2)-grille.grz(1);
gp.nx = length(grille.grx)-1;
gp.nz = length(grille.grz)-1;
gp.nsx = 11;
gp.nsz = 11;
gp.nsgx = 6;
gp.nsgz = 6;


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

for noIter=1:param.nbreitrc+param.nbreitrd+param.nbreitra
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

	D=deriv_second(m,n);
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
		if param.tomo_amp==0
			imagesc(gridx,gridz,1./reshape((ss+l_moy),n,m),'Parent',g_handles{3})
		else
			imagesc(gridx,gridz,reshape((ss+l_moy),n,m),'Parent',g_handles{3})
		end
		set(g_handles{3},'DataAspectRatio',[1 1 1],'YDir','normal')
		title('LSQR')
		xlabel(str.s119)
		ylabel(str.s120)
		if ~isempty(g_handles{1})
            caxis(g_handles{3}, g_handles{1})
        end
        colorbar('peer', g_handles{3})
		eval(['colormap(',g_handles{2},')'])
		drawnow
	end

	tomo.s=ss+l_moy;
	tomo.x=gridx;
	tomo.z=gridz;
	tomo.res=res_m;
	tomo.var_res=var_res;
	tomo.corr_tt_Lant = corr_tt_Lant;

	if param.tomo_amp==0 && noIter>=param.nbreitrd && ...
			param.nbreitrc+param.nbreitra > 0
		% update les rais
		s = tomo.s;
		% RADAR : si on a des points ou la vitesse > 0.2998 m/ns
		if param.radar == true
			s(s<3.3356) = 3.3356;
		end

%        if param.use_cont == true
%            s = appliqueCont(s, cont, gp);
%        end
		if param.nbreitra>0 && noIter>=param.nbreitrd+param.nbreitrc && ...
				noIter<param.nbreitrd+param.nbreitrc+param.nbreitra
			% on fait la derniere iteration sans faire la correction, pour pouvoir
			%    calculer l'angle du rai par rapport à l'antenne pour le
			%    diagramme de rayonnement plus facilement
            set(t_handle, 'String',[str.s226,', it ',num2str(noIter),' - ',str.s278])
            drawnow
			inWater = [data(:,3)<grille.Tx_Z_eau(3) data(:,3)<grille.Rx_Z_eau(3)];
			[L,tomo.rays,tt,corr_tt_Lant]=calcul_rc_ant2(s,gp, data(:,[1 3]), ...
				data(:,[4 6]), param.cla.corr,param.cla.diam,data(:,10:12),...
				data(:,13:15), inWater, 0);
        else
            set(t_handle, 'String',[str.s226,', it ',num2str(noIter),' - ',str.s229])
            drawnow
            [tt, tomo.rays, L] = ttcr2d(s, gp, data(:,[1 3]), data(:,[4 6]));
            corr_tt_Lant = [];
		end
		tomo.L=sparse(L);

	end
	if param.saveInvData == 1
		if noIter < param.nbreitrd || param.tomo_amp==1 || ...
				param.nbreitrc+param.nbreitra==0
			tt = L*tomo.s;
		end
		tomo.invData(noIter).res = data(:,7)-tt;
		tomo.invData(noIter).s = tomo.s;
	end
end

if param.saveInvData == 1 && param.nbreitra>0 && param.tomo_amp==0
	set(t_handle, 'String',str.s277)
	drawnow
	% on calcul les résidus avec la correction
	s = tomo.s;
	% RADAR : si on a des points ou la vitesse > 0.2998 m/ns
	if param.radar == true
		s(s<3.3356) = 3.3356;
	end
%	if param.use_cont == true
%		s = appliqueCont(s, cont, gp);
%	end
	inWater = [data(:,3)<grille.Tx_Z_eau(3) data(:,3)<grille.Rx_Z_eau(3)];
	[L,rais,tt]=calcul_rc_ant2(s,gp, data(:,[1 3]), ...
		data(:,[4 6]), param.cla.corr,param.cla.diam,data(:,10:12),...
		data(:,13:15), inWater, 0);
	tomo.invData(end).res = data(:,7)-tt;
end
	
if param.saveInvData == 1
    tomo.invData(end).date = datestr(now);
end

