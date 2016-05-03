function tomo = invLSQR(param,data,idata,grid,L,t_handle,g_handles,gv)

if ~isempty(t_handle)
    t_handle.String = 'LSQR Inversion - Starting ...';
else
    disp('LSQR Inversion - Starting ...');
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

for noIter=1:param.numItStraight + param.numItCurved
    
    if ~isempty(t_handle)
        t_handle.String = ['LSQR Inversion - Beginning of iteration ',num2str(noIter)];
        drawnow
    else
        disp(['LSQR Inversion -  Beginning of iteration ',num2str(noIter)])
    end

    if noIter == 1
        l_moy = mean(data(:,7)./sum(L,2));
    else
        l_moy = mean(tomo.s);
    end
    mta = sum(L*l_moy,2);
    dt = data(:,7) - mta;

    if noIter==1
        s_o = l_moy * ones(size(L,2),1);
    end

    [Dx,Dy,Dz] = grid.derivative(2);
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
    [x,flag,relres,iter,resvec,lsvec] = lsqr(A,b,param.tol,param.nbreiter);
    
    if max(abs(s_o./(x+l_moy) - 1))>param.dv_max
        fac = min(abs( (s_o/(param.dv_max+1)-l_moy)./x ));
        x = fac*x;
        s_o = x+l_moy;
    end
    tomo.s = x+l_moy;
    tomo.res = resvec;
    tomo.var_res = lsvec;
    
    if ~isempty(g_handles)
        if param.tomoAtt==0
            gv.plotTomo(1./tomo.s,'LSQR','Distance [m]','Elevation [m]',g_handles{3})
        else
            gv.plotTomo(tomo.s,'LSQR','Distance [m]','Elevation [m]',g_handles{3})
        end
        if ~isempty(g_handles{1}), caxis(g_handles{3},g_handles{1}), end
        colorbar('peer', g_handles{3})
        eval(['colormap(',g_handles{2},')'])
        drawnow
    end
    
    if param.tomoAtt==0 && noIter>=param.numItStraight && ...
            param.numItCurved > 0
        % update Rays
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

