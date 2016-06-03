function tomo = invGeostat(param,data,idata,grid,cm,L,t_handle,g_handles,gv)

if ~isempty(t_handle)
    t_handle.String = 'Geostatistical Inversion - Starting ...';
    drawnow
else
    disp('Geostatistical Inversion - Starting ...');
end

tomo.rays = {};
tomo.L = [];
tomo.invData = [];
if size(data,2)>=9
    tomo.no_trace = data(:,9);
end

if any(sum(data(:,1:3)==data(:,4:6),2)==2)
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
            if param.tomoAtt==0
                gv.plotTomo(1./tomo.s,'Cokriging','Distance [m]','Elevation [m]',g_handles{3})
            else
                gv.plotTomo(tomo.s,'Cokriging','Distance [m]','Elevation [m]',g_handles{3})
            end
            if ~isempty(g_handles{1}), caxis(g_handles{3},g_handles{1}), end
            colorbar('peer', g_handles{3})
            eval(['colormap(',g_handles{2},')'])
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
            if param.tomoAtt==0
                gv.plotTomo(1./tomo.s,'Cokriging','Distance [m]','Elevation [m]',g_handles{3})
            else
                gv.plotTomo(tomo.s,'Cokriging','Distance [m]','Elevation [m]',g_handles{3})
            end
            if ~isempty(g_handles{1}), caxis(g_handles{3},g_handles{1}), end
            colorbar('peer', g_handles{3})
            eval(['colormap(',g_handles{2},')'])
            
            if param.tomoAtt==0
                gv.plotTomo(1./(tomo.sgr),'Simulation','Distance [m]',[],g_handles{4},2)
            else
                gv.plotTomo(tomo.sgr,'Simulation','Distance [m]',[],g_handles{4},2)
            end
            if ~isempty(g_handles{1}), caxis(g_handles{4},g_handles{1}), end
            colorbar('peer', g_handles{4})
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