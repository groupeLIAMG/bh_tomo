function tomo = invLSQR(param,data,idata,grid,L,th,gh,gv)
% tomo = invLSQR(param,data,idata,grid,L,th,gh,gv)
%
% input
%   param : struct variable with the following members
%       use_cont:  use constraints (1) or not (0)
%       tomo_amp:  invert traveltimes (0) or amplitudes (1)
%       nbreitrd:  number of straight ray iterations to perform
%       nbreitrc:  number of curved ray iterations to perform
%       alpha:     Lagrangian of second derivative smoothing
%       saveInvData: save intermediate inversion results (1) or not (0)
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
%   idata: 
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
%   L: precomputed ray matrix (set to [] if not readily available)
%
%   th: text handle used in GUI, set to [] if using command line
%   gh: graphics handle used in GUI, set to [] if using command line
%   gv: GridViewer
%
%
% output
%   tomo : struct variable with the following members
%       s:    slowness vector
%       x:    x coord of slowness pts
%       z:    z coord of slowness pts
%       res : residuals
%       L:    ray matrix
%       rays: array of cells containing ray paths
%
%  to view the results
%     imagesc(tomo.x, tomo.z, reshape(tomo.s,length(tomo.z),length(tomo.x)))
%
%

if ~isempty(th)
    th.String = 'LSQR Inversion - Starting ...';
else
    disp('LSQR Inversion - Starting ...');
end

if ~isempty(gh)
    hideUnvisited = gh{6};
    showTxRx = gh{7};
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
    if ~isempty(th)
        th.String = 'LSQR Inversion - Straight rays ...';
        drawnow
    else
        disp('LSQR Inversion - Straight rays ...');
    end
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

[Dx,Dy,Dz] = grid.derivative(param.order);

for noIter=1:param.numItStraight + param.numItCurved
    
    if ~isempty(th)
        th.String = ['LSQR Inversion - Solving System, Iteration ',num2str(noIter)];
        drawnow
    else
        disp(['LSQR Inversion - Solving System, Iteration ',num2str(noIter)])
    end

    if param.delta == true && noIter == 1
        l_moy = 0;  % mean slowness is initially zero for delta data
    elseif noIter == 1
        l_moy = mean(data(:,7)./sum(L,2));
    else
        l_moy = mean(tomo.s);
    end
    
    mta = sum(L*l_moy,2);
    dt = data(:,7) - mta;

    if noIter==1
        s_o = l_moy * ones(size(L,2),1);
    end

    A = [L; Dx*param.alphax; Dy*param.alphay; Dz*param.alphaz];
    b = [dt;zeros(size(Dx,1),1);zeros(size(Dy,1),1);zeros(size(Dz,1),1)];
    if ~isempty(cont) && param.useCont==1
        indc = grid.getContIndices(cont);
        npx=length(cont(:,1));
        K=zeros(npx,size(L,2));
        for ndc=1:npx
            K(ndc,indc(ndc))=1;
        end
        A = [A; param.wCont*sparse(K)]; %#ok<AGROW>
        b = [b; param.wCont*cont(:,end-1)-l_moy]; %#ok<AGROW>
    end
    [x,~,~,~,resvec,lsvec] = lsqr(A,b,param.tol,param.nbreiter);
    
    if param.delta == false
        if max(abs(s_o./(x+l_moy) - 1))>param.dv_max
            fac = min(abs( (s_o/(param.dv_max+1)-l_moy)./x ));
            x = fac*x;
            s_o = x+l_moy;
        end
    else
        if max(abs(x))>param.dv_max
            x = x * param.dv_max/max(abs(x));
            s_o = x;
        end
    end
    tomo.s = x+l_moy;
    tomo.res = resvec;
    tomo.var_res = lsvec;
    
    if ~isempty(gh)
        if hideUnvisited == 1
            rd = full(sum(L));
            mask = zeros(size(tomo.s));
            mask(rd == 0) = nan;
        else
            mask = zeros(size(tomo.s));
        end
        
        if  param.tomoAtt==0 && param.delta==true
            % we want to see change in percent
            ch = 100*tomo.s ./ param.delta_prev_m;
            gv.plotTomo(mask+ch,'LSQR','Distance [m]','Elevation [m]',gh{3})
            cb_title = '% change in slowness';
        elseif param.tomoAtt==0 && param.delta == false
            gv.plotTomo(mask+1./tomo.s,'LSQR','Distance [m]','Elevation [m]',gh{3})
            cb_title = 'Velocity';
        elseif param.delta==true
            % we want to see change in percent
            ch = 100*tomo.s ./ param.delta_prev_m;
            gv.plotTomo(mask+ch,'LSQR','Distance [m]','Elevation [m]',gh{3})
            cb_title = '% change in attenuation';
        else
            gv.plotTomo(mask+tomo.s,'LSQR','Distance [m]','Elevation [m]',gh{3})
            cb_title = 'Attenuation';
        end
        if ~isempty(gh{1}), caxis(gh{3},gh{1}), end
        hcb = colorbar('peer', gh{3});
        hcb.Label.String = cb_title;
        hcb.Label.FontSize = 12;
        colormap( gh{3}, gh{2})
        if showTxRx == 1 && ~isa(grid, 'Grid3D')
            hold(gh{3}, 'on')
            plot(gh{3}, data(:,1), data(:,3), 'r+')
            plot(gh{3}, data(:,4), data(:,6), 'ro')
            hold(gh{3}, 'off')
        end
        drawnow
    end
    
    if param.tomoAtt==0 && noIter>=param.numItStraight && ...
            param.numItCurved > 0
        if any(tomo.s<0) && param.delta==false
            warndlg('Negative Slownesses: Change Inversion Parameters')
            tomo = [];
        end
        % update Rays
        if ~isempty(th)
            th.String = ['LSQR Inversion - Raytracing, Iteration ',num2str(noIter)];
            drawnow
        else
            disp(['LSQR Inversion - Raytracing, Iteration ',num2str(noIter)])
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
if ~isempty(th)
    th.String = '';
end

end

