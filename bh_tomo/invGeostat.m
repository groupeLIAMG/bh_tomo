function tomo = invGeostat(param,data,idata,grid,cm,L,t_handle,g_handles,gv)
% tomo = invGeostat(param, data, idata, grid, cm, L, t_handle, g_handles)
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
%   grid: grid object, with the following members
%       grx:     X coordinates of grid cell edges
%       grz:     Z coordinates of grid cell edges
%       cont: struct variable for the constraints, with the following  members
%           slowness:    matrix nConstrPts x 4 with the following columns
%                        z_coord x_coord slowness variance_of_slowness
%           attenuation: matrix nConstrPts x 4 with the following columns
%                        z_coord x_coord attenuation variance_of_attenuation
%
%   cm: Covariance object with following members
%
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

if ~isempty(t_handle)
    t_handle.String = 'Geostatistical Inversion - Starting ...';
    drawnow
else
    disp('Geostatistical Inversion - Starting ...');
end

if ~isempty(g_handles)
    hideUnvisited = g_handles{6};
    showTxRx = g_handles{7};
end

tomo.rays = {};
tomo.L = [];
tomo.invData = [];
if size(data,2)>=9
    tomo.no_trace = data(:,9);
end

if any(sum(data(:,1:3)==data(:,4:6),2)==3)
    uiwait(errordlg('Coincident Tx & Rx positions, aborting inversion'))
    return
end

if isempty(L)
    L = grid.getForwardStraightRays(idata);
end
tomo.x = 0.5*(grid.grx(1:end-1)+grid.grx(2:end));
tomo.z = 0.5*(grid.grz(1:end-1)+grid.grz(2:end));
if ~isempty(grid.gry)
    tomo.y = 0.5*(grid.gry(1:end-1)+grid.gry(2:end));
else
    tomo.y = [];
end

cont = [];
if param.tomoAtt==0
    if ~isempty( grid.cont.slowness )
        cont = grid.cont.slowness.data;
    end
else
    if ~isempty( grid.cont.attenuation )
        cont = grid.cont.attenuation.data;
    end
end

x = grid.getCellCenter();
Cm = cm.compute(x,x);

if ~isempty(cont) && param.useCont==1
    indc = grid.getContIndices(cont,x);

    Cm0 = Cm(indc,:);
    Cmc = Cm(indc,indc);
    Cmc = Cmc+diag(cont(:,end));
end

c0=[];
if length(data(1,:))>7 && cm.use_c0==1
    if data(:,8)~=0
        c0=data(:,8).^2;
    end
end

for noIter=1:param.numItStraight + param.numItCurved

    if ~isempty(t_handle)
        t_handle.String = ['Geostatistical Inversion - Solving System, Iteration ',num2str(noIter)];
        drawnow
    else
        disp(['Geostatistical Inversion -  Solving System, Iteration ',num2str(noIter)])
    end

    if noIter == 1
        l_moy = mean(data(:,7)./sum(L,2));
    else
        l_moy = mean(tomo.s);
    end
    mta = sum(L*l_moy,2);
    dt = data(:,7) - mta;
    doSim = 0;
    if param.doSim==1 && noIter==(param.numItStraight+param.numItCurved)
        doSim = 1;
    end

    if isempty(c0)
        Cd = L*Cm*L' + cm.nugget_d*eye(size(L,1));
    else
        Cd = L*Cm*L' + cm.nugget_d*diag(c0);
    end
    Cdm = L*Cm;

    if doSim==0
        if ~isempty(cont) && param.useCont==1
            scont = cont(:,end-1)-l_moy;
            C = [Cmc Cdm(:,indc)';Cdm(:,indc) Cd];
            C=C+eye(size(C))*1e-6;
            % dual cokriging (see Gloaguen et al 2005)
            Gamma = (C\[scont;dt])';
            m = (Gamma*[Cm0;Cdm])';
        else
            Gamma = (Cd\dt)';
            m = (Gamma*Cdm)';
        end

        if param.tomoAtt==1
            % negative attenuation set to zero
            m(m<-l_moy) = -l_moy;
        end
        tomo.s = m+l_moy;

        if ~isempty(g_handles)
            if hideUnvisited == 1
                rd = full(sum(L));
                mask = zeros(size(tomo.s));
                mask(rd == 0) = nan;
            else
                mask = zeros(size(tomo.s));
            end

            if param.tomoAtt==0
                gv.plotTomo(mask+1./tomo.s,'Cokriging','Distance [m]','Elevation [m]',g_handles{3})
            else
                gv.plotTomo(mask+tomo.s,'Cokriging','Distance [m]','Elevation [m]',g_handles{3})
            end
            if ~isempty(g_handles{1}), caxis(g_handles{3},g_handles{1}), end
            colorbar('peer', g_handles{3})
            eval(['colormap(',g_handles{2},')'])
            if showTxRx == 1 && ~isa(grid, 'Grid3D')
                hold(g_handles{3}, 'on')
                plot(g_handles{3}, data(:,1), data(:,3), 'r+')
                plot(g_handles{3}, data(:,4), data(:,6), 'ro')
                hold(g_handles{3}, 'off')
            end
            drawnow
        end
    else
        if ~isempty(cont) && param.useCont==1
            scont = cont(:,end-1)-l_moy;
            C = [Cmc Cdm(:,indc)';Cdm(:,indc) Cd];
            C=C+eye(size(C))*1e-6;
            % primal cokriging (see Gloaguen et al 2005)
            Lambda = (C\[Cm0;Cdm])';
            m = (Lambda*[scont;dt]);
        else
            Lambda = (Cd\Cdm)';
            m = (Lambda*dt);
        end

        mSim = zeros(length(m),param.nSim);
        G = grid.preFFTMA(cm.covar);
        for ns=1:param.nSim
            ms = grid.FFTMA(G);
            ms = ms(:);

            % data computed with simulated model
            ds = L*ms + randn(length(dt),1)*cm.nugget_d;
            % cokriging estimator of the simulated data
            if ~isempty(cont) && param.useCont==1
                mss = Lambda*[ms(indc);ds];
            else
                mss = Lambda*ds;
            end
            % conditional simulated model  (eq 34 of Giroux et al 2007)
            mSim(:,ns) = m + (ms-mss);

        end

        [~,~,diff1_min]=choixsimu(L,mSim,dt,c0);
        if ~isempty(cont) && param.useCont==1
            ms = mSim(:,diff1_min);
        else
            ms = deformationGraduelle(mSim,L,c0,dt,t_handle);
        end

        tomo.s = m+l_moy;
        tomo.simu = mSim+l_moy;
        tomo.sgr = ms+l_moy;
        tomo.lmoy=l_moy;
        tomo.diff1_min=diff1_min;

        if ~isempty(g_handles)
            if hideUnvisited == 1
                rd = full(sum(L));
                mask = zeros(size(tomo.s));
                mask(rd == 0) = nan;
            else
                mask = zeros(size(tomo.s));
            end

            if param.tomoAtt==0
                gv.plotTomo(mask+1./tomo.s,'Cokriging','Distance [m]','Elevation [m]',g_handles{3})
            else
                gv.plotTomo(mask+tomo.s,'Cokriging','Distance [m]','Elevation [m]',g_handles{3})
            end
            if ~isempty(g_handles{1}), caxis(g_handles{3},g_handles{1}), end
            colorbar('peer', g_handles{3})
            eval(['colormap(',g_handles{2},')'])

            if param.tomoAtt==0
                gv.plotTomo(mask+1./(tomo.sgr),'Simulation','Distance [m]',[],g_handles{4},2)
            else
                gv.plotTomo(mask+tomo.sgr,'Simulation','Distance [m]',[],g_handles{4},2)
            end
            if ~isempty(g_handles{1}), caxis(g_handles{4},g_handles{1}), end
            colorbar('peer', g_handles{4})
            if showTxRx == 1 && ~isa(grid, 'Grid3D')
                hold(g_handles{3}, 'on')
                hold(g_handles{4}, 'on')
                plot(g_handles{3}, data(:,1), data(:,3), 'r+')
                plot(g_handles{3}, data(:,4), data(:,6), 'ro')
                plot(g_handles{4}, data(:,1), data(:,3), 'r+')
                plot(g_handles{4}, data(:,4), data(:,6), 'ro')
                hold(g_handles{3}, 'off')
                hold(g_handles{4}, 'off')
            end
            drawnow
        end
    end

    if param.tomoAtt==0 && noIter>=param.numItStraight && ...
            param.numItCurved > 0
        if any(tomo.s<0)
            warndlg('Negative Slownesses: Change Inversion Parameters')
            tomo = [];
        end
        % update Rays
        if ~isempty(t_handle)
            t_handle.String = ['Geostatistical Inversion - Raytracing, Iteration ',num2str(noIter)];
            drawnow
        else
            disp(['Geostatistical Inversion -  Raytracing, Iteration ',num2str(noIter)])
        end
        [~,tomo.rays,L] = grid.raytrace(tomo.s,data(:,1:3),data(:,4:6));

    end
    if param.saveInvData == 1
        tt = L*tomo.s;
        tomo.invData(noIter).res = data(:,7)-tt;
        tomo.invData(noIter).s = tomo.s;
    end
end
tomo.L = L;

if param.saveInvData == 1
    tomo.invData(end).date = datestr(now);
end
if ~isempty(t_handle)
    t_handle.String = '';
end

end
