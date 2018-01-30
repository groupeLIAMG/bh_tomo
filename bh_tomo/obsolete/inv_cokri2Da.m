function tomo = inv_cokri2Da(param, data, covar, grille, L, t_handle, g_handles)
% tomo = inv_cokri2Da(param, data, covar, grille, L, t_handle, g_handles)
%
%

% Copyright (C) 2010 Bernard Giroux
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
    [L,gridx,gridz] = Lsr2da(data(:,[1 3]),data(:,[4 6]),grille.grx,...
        grille.grz);
else
    %Calcul des vecteurs position des centres des blocs
    gridx = 0.5*(grille.grx(1:end-1)+grille.grx(2:end));
    gridz = 0.5*(grille.grz(1:end-1)+grille.grz(2:end));
end

np=size(L,2)/2;

Lx = L(:,1:np);
Lz = L(:,(1+np):end);


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


%cont = grille.cont.slowness;
cont = [];
if param.tomo_amp==0
    if ~isempty( grille.cont.slowness )
        cont = grille.cont.slowness;
    end
%else
%    if ~isempty( grille.cont.attenuation )
%        cont = grille.cont.attenuation.data;
%    end
end

% 	if isfield( grille.cont, 'variance_s')
% 		variance_cont = grille.cont.variance_s;
% 	end


n=length(gridz);
m=length(gridx);
nCells=m*n;
np=size(L,2);
nt=size(L,1);
ind=1:np;



%Calcul de la matrice de covariance
if param.ival==2
    x=grille2d(min(gridz),gridz(2)-gridz(1),n,min(gridx),gridx(2)-gridx(1),m);
    
    indc = calcul_kss_aniso(x, covar, cont, param.use_cont, nCells, t_handle);
end

c0=[];
if length(data(1,:))>7 && covar.use_c0==1
    if data(:,8)~=0
        c0=data(:,8).^2;
    end
end



nbreitJ = 3;

rayLength = sqrt(Lx.^2 + Lz.^2);
rayLength = sum(rayLength,2);
s0 = mean(data(:,7)./rayLength);
t0 = s0*rayLength;
s0 = s0 + zeros(size(Lx,2),1);
xi0 = ones(size(Lz,2),1);
e0 = [s0; xi0];


ordreNorme = 2;
error = covar.nugget_t;
%facReducC = 0.1;


modeJ = 1;

for noIter=1:param.nbreitrd+param.nbreitrc+param.nbreitra
    if ~isempty(t_handle)
        set(t_handle, 'String',[str.s225,' - ',str.s228,' ',num2str(noIter)])
        drawnow
		else
			disp([str.s225,' - ',str.s228,' ',num2str(noIter)]); drawnow
    end

    mode = 1;
    if param.mode==2 && noIter==(param.nbreitrd+param.nbreitrc+param.nbreitra)
        mode = 2;
    end

    res = norm(data(:,7) - t0, ordreNorme);
    resJ = zeros(nbreitJ,1)-1;
    resJ(1) = res;

    for itJ=1:nbreitJ
			if mode==2 && itJ==nbreitJ
				modeJ = 2;
			end
        %disp(itJ)

        if isempty(corr_tt_Lant)
            dt = data(:,7) - t0;
        else
            dt = data(:,7) + corr_tt_Lant - t0;
        end
    
        J = calculJ(L, e0);
        
        % center w/respect to J
%        moy = mean(dt./sum(J,2));
%        mt = sum(J*moy,2);
%        dt = dt - mt;
        
        if ~isempty( cont.data ) && param.use_cont   %%YH cont->cont.data
            if mode==1 && ~isempty(t_handle)
                set(t_handle,'String',[str.s225,', it ',num2str(noIter),' - ',str.s232]); drawnow
            elseif mode==2 && ~isempty(t_handle)
                set(t_handle,'String',[str.s225,', it ',num2str(noIter),' - ',str.s233]); drawnow
            end
            scont = cont.data(:,3);
            if ~isempty(cont.data_xi)
                scont = [scont; cont.data_xi(:,3)];
            end
            scont = scont-e0(indc);
            [e_sim,de] = simu_contrainte2a(J,dt,covar,indc,nt,...
                gridx,gridz,param.nbresimu,scont,modeJ,c0);
        else
            if mode==1 &&  ~isempty(t_handle)
                set(t_handle,'String',[str.s225,', it ',num2str(noIter),' - ',str.s234]); drawnow
            elseif mode==2 &&  ~isempty(t_handle)
                set(t_handle,'String',[str.s225,', it ',num2str(noIter),' - ',str.s234b]); drawnow
            end
            [e_sim,de] = cokri_simua(J,dt,covar,nt,gridx,gridz,param.nbresimu,modeJ,c0);
        end
    
%        ds = ds+moy;
        
        s = s0 + de(1:nCells);
        xi = xi0 + de((nCells+1):end);
        t = Lx.^2 + Lz.^2.*kron(ones(nt,1),xi'.^2);
        t = kron(ones(nt,1),s') .* sqrt(t);
        t = sum(t,2);
        res = norm(data(:,7)-t, ordreNorme);
				
				if modeJ==2
					e_sim = [e_sim(1:nCells,:)+kron(ones(1,param.nbresimu),s0);
									 e_sim((nCells+1):end,:)+kron(ones(1,param.nbresimu),xi0)];
				end
        t0 = t;
        s0 = s;
        xi0 = xi;
        e0 = [s0; xi0];
        resJ(itJ) = res;
        if res < error
            break
        end
    end

    
    
    
    
    
    if mode==2

        [diff1,~,diff1_min] = choixsimua(L,e_sim,dt,c0);
        
%        if param.use_cont==false || isempty( cont )
%            sgr = deformationGraduelle(e_sim, L, c0, dt);
%        else
            sgr = e_sim(:,diff1_min);
%        end
        if ~isempty(g_handles)
            imagesc(gridx,gridz,1./reshape((ss+l_moy),length(gridz),length(gridx)),'Parent',g_handles{3})
            set(g_handles{3},'DataAspectRatio',[1 1 1],'YDir','normal')
            title(g_handles{3},str.s230,'FontSize',14)
            xlabel(g_handles{3},str.s119)
            ylabel(g_handles{3},str.s120)
            if ~isempty(g_handles{1}), caxis(g_handles{3},g_handles{1}), end
            colorbar('peer',g_handles{3})
            imagesc(gridx,gridz,1./reshape((sgr+l_moy),length(gridz),length(gridx)),'Parent',g_handles{4})
            set(g_handles{4},'DataAspectRatio',[1 1 1],'YDir','normal')
            title(g_handles{4},'\xi','FontSize',14)
            xlabel(g_handles{4},str.s119)
            ylabel(g_handles{4},str.s120)
            if ~isempty(g_handles{1}), caxis(g_handles{4},g_handles{1}), end
            colorbar('peer', g_handles{4})
            eval(['colormap(',g_handles{2},')'])
            drawnow
        end
        
        tomo.s    = s0;
        tomo.xi   = xi0;
        tomo.simu = e_sim;
        tomo.sgr  = sgr;
        tomo.x    = gridx;
        tomo.z    = gridz;
        %tomo.lmoy = l_moy;
        tomo.diff1 = diff1;
        tomo.diff1_min = diff1_min;
        %tomo.dt=dt;
        tomo.corr_tt_Lant = corr_tt_Lant;
    else
        tomo.s  = s0;
        tomo.xi = xi0;
        tomo.x  = gridx;
        tomo.z  = gridz;
        %tomo.dt=dt;
        tomo.corr_tt_Lant = corr_tt_Lant;
        sz = s0 .* xi0;

        if ~isempty(g_handles)
            imagesc(gridx,gridz,1./reshape(sz,n,m),'Parent',g_handles{3})
            set(g_handles{3},'DataAspectRatio',[1 1 1],'YDir','normal')
            title(g_handles{3},str.s285b,'FontSize',14)
            xlabel(g_handles{3},str.s119)
            ylabel(g_handles{3},str.s120)
            if ~isempty(g_handles{1}), caxis(g_handles{3},g_handles{1}), end
            colorbar('peer', g_handles{3})
        
            imagesc(gridx,gridz,reshape(xi0,n,m),'Parent',g_handles{4})
            set(g_handles{4},'DataAspectRatio',[1 1 1],'YDir','normal')
            title(g_handles{4},'\xi','FontSize',14)
            xlabel(g_handles{4},str.s119)
            ylabel(g_handles{4},str.s120)
            %if ~isempty(g_handles{1}), caxis(g_handles{4},g_handles{1}), end
            colorbar('peer', g_handles{4})
            eval(['colormap(',g_handles{2},')'])
            drawnow
        end        
    end
    
    if noIter>=param.nbreitrd && param.nbreitrc+param.nbreitra > 0
        
        % update les rais
        s = tomo.s;
        % RADAR : si on a des points ou la vitesse > 0.2998 m/ns
				if param.radar == true
					s(s<3.3356) = 3.3356;
        end
				
        if param.nbreitra>0 && noIter>=param.nbreitrd+param.nbreitrc && ...
                noIter<param.nbreitrd+param.nbreitrc+param.nbreitra
            % on fait la derniere iteration sans faire la correction, pour pouvoir
            %    calculer l'angle du rai par rapport a l'antenne pour le
            %    diagramme de rayonnement plus facilement
            if ~isempty(t_handle)
                set(t_handle, 'String',[str.s225,', it ',num2str(noIter),' - ',str.s278])
                drawnow
            end
            inWater = [data(:,3)<grille.Tx_Z_eau(3) data(:,3)<grille.Rx_Z_eau(3)];
            [L,tomo.rays,tt,corr_tt_Lant] = calcul_rc_ant2(s,gp, data(:,[1 3]), ...
                data(:,[4 6]), param.cla.corr,param.cla.diam,data(:,10:12),...
                data(:,13:15), inWater, 1, tomo.xi);
        else
            if ~isempty(t_handle)
                set(t_handle, 'String',[str.s225,', it ',num2str(noIter),' - ',str.s229])
                drawnow
            end
            [tt, tomo.rays, L] = ttcr2da(tomo.s, tomo.xi, gp, data(:,[1 3]), data(:,[4 6]));

             corr_tt_Lant = [];
        end
        Lx = L(:,1:nCells);
        Lz = L(:,(1+nCells):np);
        tomo.L = sparse(L);
        
    end
    
    if param.saveInvData == 1
        if noIter <= param.nbreitrd+param.nbreitrc || ...
                param.nbreitrc+param.nbreitra==0
            tt = Lx.^2 + Lz.^2.*kron(ones(nt,1),(tomo.xi'.^2));
            tt = sqrt(tt) .* kron(ones(nt,1),(tomo.s'));
            tt = sum(tt,2);
        end
        tomo.invData(noIter).res = data(:,7)-tt;
        tomo.invData(noIter).resJ = resJ;
        tomo.invData(noIter).s  = tomo.s;
        tomo.invData(noIter).xi = tomo.xi;
    end
end

if param.saveInvData == 1 && param.nbreitra>0 && param.tomo_amp==0
    if ~isempty(t_handle)
        set(t_handle, 'String',str.s277)
        drawnow
    end
    % on calcul les residus avec la correction (s est deja defini)
		if param.radar == true
			s(s<3.3356) = 3.3356;
		end

    inWater = [data(:,3)<grille.Tx_Z_eau(3) data(:,3)<grille.Rx_Z_eau(3)];
    [~,tomo.rays,tt]=calcul_rc_ant2(s,gp, data(:,[1 3]), ...
        data(:,[4 6]), param.cla.corr,param.cla.diam,data(:,10:12),...
        data(:,13:15), inWater, 1, tomo.xi);
    tomo.invData(end).res = data(:,7)-tt;
end

if param.saveInvData == 1
    tomo.invData(end).date = datestr(now);
end

