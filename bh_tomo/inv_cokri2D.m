function tomo = inv_cokri2D(param, data, covar, grille, L, t_handle, g_handles)
% tomo = inv_cokri2D(param, data, covar, grille, L, t_handle, g_handles)
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
%   covar: struct variable with following members
%       model : Each row of this matrix describes a different elementary structure.
%               The first column is a code for the model type, the d following
%               columns give the ranges along the different coordinates and the
%               subsequent columns give rotation angles (a maximum of three).
%               For more details on how to specify rotations, type help trans.
%               The codes for the current models are:
%                  1: nugget effect
%                  2: exponential model
%                  3: gaussian model
%                  4: spherical model
%                  5: linear model
%               Note: a linear model is specified by arbitrary ranges and a sill
%                     such that sill/range gives the desired slope in the direction
%                     condidered.
%       c : The (rp x p) coefficient matrix of the coregionalization model.
%           Position (i,j) in each submatrix of size p x p give the sill of the
%           elementary component for each cross-variogram (variogram) between
%           variable i and variable j.
%       use_c0 : use experimental covariance (1) or not (0)
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
corr_tt_Lant = [];

if any(sum(data(:,[1 3])==data(:,[4 6]),2)==2)
    uiwait(errordlg('Coincident Tx & Rx positions, aborting inversion'))
    return
end

if isempty(L)
    [L,gridx,gridz] = Lsr2d(data(:,[1 3]),data(:,[4 6]),grille.grx,...
        grille.grz);
else
    %Calcul des vecteurs position des centres des blocs
    gridx = 0.5*(grille.grx(1:end-1)+grille.grx(2:end));
    gridz = 0.5*(grille.grz(1:end-1)+grille.grz(2:end));
end

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
if param.tomo_amp==0
    if ~isempty( grille.cont.slowness )
        cont = grille.cont.slowness.data;
    end
else
    if ~isempty( grille.cont.attenuation )
        cont = grille.cont.attenuation.data;
    end
end


n=length(gridz);
m=length(gridx);
np=length(L(1,:));
nt=length(L(:,1));
ind=1:np;

%Calcul de la matrice de covariance
if param.ival==2
    x=grille2d(min(gridz),gridz(2)-gridz(1),n,min(gridx),gridx(2)-gridx(1),m);
    %	np=calcul_kss_global(np,x,ind,covar,gridx,gridz,m,n,t_handle);
    calcul_kss_global(np,x,ind,covar,gridx,gridz,m,n,t_handle);
    if ~isempty( cont ) && param.use_cont
        indc = zeros(size(cont,1),1);
        for i=1:size(cont,1)
            ind11=findnear(cont(i,1),x(:,1));
            ind22=findnear(cont(i,2),x(ind11,2));
            indc(i)=min(ind11(ind22));
        end
        clear ind11 ind22
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        %Terme de gauche
        ks0=kss(indc,:);
        %Covariance des donnees de contrainte
        
        % contraintes dures
        Css = kss(indc,indc);
        
        % contraintes souples
        % la variance est contenue ds la 4e colonne
        if size(cont,2)==4
            Css = Css+diag(cont(:,4));
        end
    end
end

c0=[];
if length(data(1,:))>7 && covar.use_c0==1
    if data(:,8)~=0
        c0=data(:,8).^2;
    end
end


for noIter=1:param.nbreitrd+param.nbreitrc+param.nbreitra
    
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
    if isempty(corr_tt_Lant)
        dt = data(:,7) - mta;
    else
        dt = data(:,7) + corr_tt_Lant - mta;
    end
    mode = 1;
    if param.mode==2 && noIter==(param.nbreitrd+param.nbreitrc+param.nbreitra)
        mode = 2;
    end
    
    if ~isempty( cont ) && param.use_cont
        if mode==1 && ~isempty(t_handle)
            set(t_handle,'String',[str.s225,', it ',num2str(noIter),' - ',str.s232]); drawnow
        elseif mode==2 && ~isempty(t_handle)
            set(t_handle,'String',[str.s225,', it ',num2str(noIter),' - ',str.s233]); drawnow
        end
        scont = cont(:,3)-l_moy;
        [S_sim,ss] = simu_contrainte2(L, dt, covar, indc, nt, ...
            gridx, gridz, param.nbresimu, scont, mode, c0);
    else
        if mode==1 && ~isempty(t_handle)
            set(t_handle,'String',[str.s225,', it ',num2str(noIter),' - ',str.s234]); drawnow
        elseif mode==2 && ~isempty(t_handle)
            set(t_handle,'String',[str.s225,', it ',num2str(noIter),' - ',str.s234b]); drawnow
        end
        [S_sim,ss] = cokri_simu(L, dt, covar, nt, gridx, gridz, ...
            param.nbresimu, mode, c0);
    end
    
    if param.tomo_amp==1
        ss(ss<-l_moy) = -l_moy;   % on fixe l'attŽnuation a zŽro
    end
    if mode==2
        [diff1,~,diff1_min]=choixsimu(L(:,ind),S_sim(ind,:),dt,c0);
        
        if param.use_cont==false || isempty( cont )
            sgr = deformationGraduelle(S_sim, L, c0, dt);
        else
            sgr = S_sim(:,diff1_min);
        end
        
        if  ~isempty(g_handles)
            set(g_handles{3},'DataAspectRatio',[1 1 1],'YDir','normal')
            title(g_handles{3},str.s230,'FontSize',14)
            xlabel(g_handles{3},str.s119)
            ylabel(g_handles{3},str.s120)
            if ~isempty(g_handles{1}), caxis(g_handles{3}, g_handles{1}), end
            colorbar('peer',g_handles{3})
            if param.tomo_amp==0
                imagesc(gridx,gridz,1./reshape((sgr+l_moy),length(gridz),length(gridx)),'Parent',g_handles{4})
            else
                sgr(sgr<-l_moy) = -l_moy;
                imagesc(gridx,gridz,reshape((sgr+l_moy),length(gridz),length(gridx)),'Parent',g_handles{4})
            end
            set(g_handles{4},'DataAspectRatio',[1 1 1],'YDir','normal')
            title(g_handles{4},str.s231,'FontSize',14)
            xlabel(g_handles{4},str.s119)
            ylabel(g_handles{4},str.s120)
            if ~isempty(g_handles{1}), caxis(g_handles{4}, g_handles{1}), end
            colorbar('peer', g_handles{4})
            eval(['colormap(',g_handles{2},')'])
            drawnow
        end
        tomo.s=ss+l_moy;
        tomo.simu=S_sim;%+l_moy;
        tomo.sgr=sgr;
        tomo.x=gridx;
        tomo.z=gridz;
        tomo.lmoy=l_moy;
        tomo.diff1=diff1;
        tomo.diff1_min=diff1_min;
        %tomo.dt=dt;
        tomo.corr_tt_Lant = corr_tt_Lant;
    else
        if  ~isempty(g_handles)
            if param.tomo_amp==0
                imagesc(gridx,gridz,1./reshape((ss+l_moy),n,m),'Parent',g_handles{3})
            else
                imagesc(gridx,gridz,reshape((ss+l_moy),n,m),'Parent',g_handles{3})
            end
            set(g_handles{3},'DataAspectRatio',[1 1 1],'YDir','normal')
            title(g_handles{3},str.s230,'FontSize',14)
            xlabel(g_handles{3},str.s119)
            ylabel(g_handles{3},str.s120)
            if ~isempty(g_handles{1}), caxis(g_handles{3},g_handles{1}), end
            colorbar('peer', g_handles{3})
            eval(['colormap(',g_handles{2},')'])
            drawnow
        end
        
        tomo.s=ss+l_moy;
        tomo.x=gridx;
        tomo.z=gridz;
        %tomo.dt=dt;
        tomo.corr_tt_Lant = corr_tt_Lant;
    end
    tomo.L = sparse(L);
    
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
            % on fait la derniere mise à jour des rais sans faire la correction, pour pouvoir
            %    calculer l'angle du rai par rapport à l'antenne pour le
            %    diagramme de rayonnement plus facilement
            if ~isempty(t_handle)
                set(t_handle, 'String',[str.s225,', it ',num2str(noIter),' - ',str.s278])
                drawnow
            end
            inWater = [data(:,3)<grille.Tx_Z_eau(3) data(:,3)<grille.Rx_Z_eau(3)];
            [L,tomo.rays,tt,corr_tt_Lant]=calcul_rc_ant2(s,gp, data(:,[1 3]), ...
                data(:,[4 6]), param.cla.corr,param.cla.diam,data(:,10:12),...
                data(:,13:15), inWater, 0);
        else
            if ~isempty(t_handle)
                set(t_handle, 'String',[str.s225,', it ',num2str(noIter),' - ',str.s229])
                drawnow
            end
            [tt, tomo.rays, L] = ttcr2d(s, gp, data(:,[1 3]), data(:,[4 6]));
            corr_tt_Lant = [];
        end
        
    end
    
    if param.saveInvData == 1
        if noIter <= param.nbreitrd+param.nbreitrc || param.tomo_amp==1 || ...
                param.nbreitrc+param.nbreitra==0
            tt = tomo.L*tomo.s;
        end
        tomo.invData(noIter).res = data(:,7)-tt;
        tomo.invData(noIter).s = tomo.s;
    end
end
tomo.L = sparse(L);

if param.saveInvData == 1 && param.nbreitra>0 && param.tomo_amp==0
    if ~isempty(t_handle)
        set(t_handle, 'String',str.s277)
        drawnow
	end

	% on calcul les résidus avec la correction (s est déjà défini)
	if param.radar == true
		s(s<3.3356) = 3.3356;
	end
    %	if param.use_cont == true
    %		s = appliqueCont(s, cont, gp);
    %	end
    inWater = [data(:,3)<grille.Tx_Z_eau(3) data(:,3)<grille.Rx_Z_eau(3)];
    [~,~,tt]=calcul_rc_ant2(s,gp, data(:,[1 3]), ...
        data(:,[4 6]), param.cla.corr,param.cla.diam,data(:,10:12),...
        data(:,13:15), inWater, 0);
    tomo.invData(end).res = data(:,7)-tt;
elseif param.tomo_amp==1
    tomo.invData(end).res = data(:,7)-tt;
end

if param.saveInvData == 1
    tomo.invData(end).date = datestr(now);
end
