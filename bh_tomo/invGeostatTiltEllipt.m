function tomo = invGeostatTiltEllipt(param,data,idata,grid,cm,L,th,gh,gv )

if ~isempty(th)
    th.String = 'Geostatistical Inversion - Starting ...';
else
    disp('Geostatistical Inversion - Starting ...');
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
    L = grid.getForwardStraightRays(idata,[],[],[],true);
end
tomo.x = 0.5*(grid.grx(1:end-1)+grid.grx(2:end));
tomo.z = 0.5*(grid.grz(1:end-1)+grid.grz(2:end));
if ~isempty(grid.gry)
    tomo.y = 0.5*(grid.gry(1:end-1)+grid.gry(2:end));
else
    tomo.y = [];
end


nCells=size(L,2)/2;
np=size(L,2);
nt=size(L,1);

Lx = L(:,1:nCells);
Lz = L(:,(1+nCells):end);

cont = [];
if param.tomoAtt==0
    if ~isempty( grid.cont.slowness )
        cont = grid.cont.slowness;
    end
end

x = grid.getCellCenter();
Cm = cm.compute(x,x);

if ~isempty(cont) && param.useCont==1
    indc = grid.getContIndices(cont,x);
    Cm0 = Cm(indc,:);
    Cmc = Cm(indc,indc);
    if size(cont.data,2)==4 && size(cont.data_xi,2)==4 && size(cont.data_theta,2)==4
        Cmc = Cmc + diag([cont.data(:,4); cont.data_xi(:,4); cont.data_theta(:,4)]);
    end
end

c0=[];
if length(data(1,:))>7 && cm.use_c0==1
    if data(:,8)~=0
        c0=data(:,8).^2;
    end
end

rayLength = sqrt(Lx.^2 + Lz.^2);
rayLength = sum(rayLength,2);
s0 = mean(data(:,7)./rayLength);
t0 = s0*rayLength;
s0 = s0 + zeros(size(Lx,2),1);
xi0 = ones(size(Lz,2),1)+0.001;     % add 1/1000 so that J_th ~= 0
te0 = zeros(size(Lz,2),1)+0.0044;   % add a quarter of a degree so that J_th ~= 0
e0 = [s0; xi0; te0];

nbreitJ = 3;
normOrder = 2;
error = cm.nugget_d;
modeJ = 1;

for noIter=1:param.numItStraight + param.numItCurved

    if ~isempty(th)
        th.String = ['Geostatistical Inversion - Solving System, Iteration ',num2str(noIter)];
        drawnow
    else
        disp(['Geostatistical Inversion -  Solving System, Iteration ',num2str(noIter)])
    end


    doSim = 0;
    if param.doSim==1 && noIter==(param.numItStraight+param.numItCurved)
        doSim = 1;
    end


    res = norm(data(:,7) - t0, normOrder);
    resJ = zeros(nbreitJ,1)-1;
    resJ(1) = res;

    if doSim==0
        for itJ=1:nbreitJ

            dt = data(:,7) - t0;

            J = calculJ2(L, e0);

            if isempty(c0)
                Cd = J*Cm*J' + cm.nugget_d*eye(size(L,1));
            else
                Cd = J*Cm*J' + cm.nugget_d*diag(c0);
            end
            Cdm = J*Cm;

            % center w/respect to J
            %        moy = mean(dt./sum(J,2));
            %        mt = sum(J*moy,2);
            %        dt = dt - mt;

            if ~isempty( cont.data ) && param.useCont==1
                scont = [];
                if ~isempty(cont.data)
                    scont = [scont; cont.data(:,3)]; %#ok<AGROW>
                end
                if ~isempty(cont.data_xi)
                    scont = [scont; cont.data_xi(:,3)]; %#ok<AGROW>
                end
                if ~isempty(cont.data_theta)
                    scont = [scont; cont.data_theta(:,3)]; %#ok<AGROW>
                end
                scont = scont-e0(indc);
                C = [Cmc Cdm(:,indc)';Cdm(:,indc) Cd];
                C = C+eye(size(C))*1e-6;

                Gamma = (C\[scont;dt])';
                de = (Gamma*[Cm0;Cdm])';
            else
                Gamma = (Cd\dt)';
                de = (Gamma*Cdm)';
            end

            s = s0 + de(1:nCells);
            xi = xi0 + de((nCells+1):(2*nCells));
            te = te0 + de((2*nCells+1):end);
            co = kron(ones(nt,1),cos(te)');
            si = kron(ones(nt,1),sin(te)');
            t = (Lx.*co + Lz.*si).^2 + (Lz.*co - Lx.*si).^2.*kron(ones(nt,1),xi'.^2);
            t = kron(ones(nt,1),s') .* sqrt(t);
            t = sum(t,2);
            res = norm(data(:,7)-t, normOrder);

            t0 = t;
            s0 = s;
            xi0 = xi;
            te0 = te;
            e0 = [s0; xi0; te0];
            resJ(itJ) = res;
            if res < error
                break
            end
        end

        tomo.s  = s0;
        tomo.xi = xi0;
        tomo.th = te0;

    else

        for itJ=1:nbreitJ
            if doSim==1 && itJ==nbreitJ
                modeJ = 2;
            end

            dt = data(:,7) - t0;

            J = calculJ2(L, e0);

            if isempty(c0)
                Cd = J*Cm*J' + cm.nugget_d*eye(size(L,1));
            else
                Cd = J*Cm*J' + cm.nugget_d*diag(c0);
            end
            Cdm = J*Cm;

            % center w/respect to J
            %        moy = mean(dt./sum(J,2));
            %        mt = sum(J*moy,2);
            %        dt = dt - mt;

            if ~isempty( cont.data ) && param.useCont==1
                scont = [];
                if ~isempty(cont.data)
                    scont = [scont; cont.data(:,3)]; %#ok<AGROW>
                end
                if ~isempty(cont.data_xi)
                    scont = [scont; cont.data_xi(:,3)]; %#ok<AGROW>
                end
                if ~isempty(cont.data_theta)
                    scont = [scont; cont.data_theta(:,3)]; %#ok<AGROW>
                end
                scont = scont-e0(indc);
                C = [Cmc Cdm(:,indc)';Cdm(:,indc) Cd];
                C = C+eye(size(C))*1e-6;

                Gamma = (C\[scont;dt])';
                de = (Gamma*[Cm0;Cdm])';
            else
                Gamma = (Cd\dt)';
                de = (Gamma*Cdm)';
            end

            s = s0 + de(1:nCells);
            xi = xi0 + de((nCells+1):(2*nCells));
            te = te0 + de((2*nCells+1):end);
            co = kron(ones(nt,1),cos(te)');
            si = kron(ones(nt,1),sin(te)');
            t = (Lx.*co + Lz.*si).^2 + (Lz.*co - Lx.*si).^2.*kron(ones(nt,1),xi'.^2);
            t = kron(ones(nt,1),s') .* sqrt(t);
            t = sum(t,2);
            res = norm(data(:,7)-t, normOrder);

            t0 = t;
            s0 = s;
            xi0 = xi;
            te0 = te;
            e0 = [s0; xi0; te0];
            resJ(itJ) = res;
            if res < error
                break
            end

            if modeJ==2

                e_sim = zeros(length(e0),param.nSim);
                Gs = grid.preFFTMA(cm.covar);
                Gxi = grid.preFFTMA(cm.covar_xi);
                Gth = grid.preFFTMA(cm.covar_tilt);

                for ns=1:param.nSim
                    ms = grid.FFTMA(Gs);
                    xis = grid.FFTMA(Gxi);
                    ths = grid.FFTMA(Gth);
                    ms = ms(:);
                    xis = xis(:);
                    ths = ths(:);

                    % data computed with simulated model
                    co = kron(ones(nt,1),cos(ths)');
                    si = kron(ones(nt,1),sin(ths)');
                    t = (Lx.*co + Lz.*si).^2 + (Lz.*co - Lx.*si).^2.*kron(ones(nt,1),xis'.^2);
                    t = kron(ones(nt,1),ms') .* sqrt(t);
                    t = sum(t,2);

                    ds = t + randn(length(dt),1)*cm.nugget_d;
                    % cokriging estimator of the simulated data
                    es = [ms; xis; ths];
                    if ~isempty(cont) && param.useCont==1
                        ess = Lambda*[es(indc);ds];
                    else
                        ess = Lambda*ds;
                    end
                    % conditional simulated model  (eq 34 of Giroux et al 2007)
                    e_sim(:,ns) = de + (es-ess);
                end

                e_sim = [e_sim(1:nCells,:)+kron(ones(1,param.nSim),s0);
                    e_sim((nCells+1):2*nCells,:)+kron(ones(1,param.nSim),xi0);
                    e_sim((2*nCells+1):end,:)+kron(ones(1,param.nSim),te0);];
            end

        end

        [~,~,diff1_min] = choixsimuaa(L,e_sim,dt,c0);
        if ~isempty(cont) && param.useCont==1
            sgr = e_sim(:,diff1_min);
        else
            sgr = deformationGraduelle(e_sim,L,c0,dt,th);
        end

        tomo.s  = s0;
        tomo.xi = xi0;
        tomo.th = te0;
        tomo.simu = e_sim;
        tomo.sgr  = sgr;
        tomo.diff1_min = diff1_min;
    end

    if ~isempty(gh)
        if hideUnvisited == 1
            rd = full(sum(sqrt(Lx.^2+Lz.^2)));
            mask = zeros(size(tomo.s));
            mask(rd == 0) = nan;
        else
            mask = zeros(size(tomo.s));
        end

        gv.plotTomo(mask+1./(s0.*xi0),'V_z','Distance [m]','Elevation [m]',gh{3})
        if ~isempty(gh{1}), caxis(gh{3},gh{1}), end
        colorbar('peer', gh{3})

        gv.plotTomo(mask+xi0,'\xi','Distance [m]','',gh{4})
        colorbar('peer', gh{4})

        gv.plotTomo(mask+te0*180/pi,'\theta','Distance [m]','',gh{5})
        colorbar('peer', gh{5})

        colormap( gh{3}, gh{2})
        if showTxRx == 1 && ~isa(grid, 'Grid3D')
            hold(gh{3}, 'on')
            hold(gh{4}, 'on')
            hold(gh{5}, 'on')
            plot(gh{3}, data(:,1), data(:,3), 'r+')
            plot(gh{3}, data(:,4), data(:,6), 'ro')
            plot(gh{4}, data(:,1), data(:,3), 'r+')
            plot(gh{4}, data(:,4), data(:,6), 'ro')
            plot(gh{5}, data(:,1), data(:,3), 'r+')
            plot(gh{5}, data(:,4), data(:,6), 'ro')
            hold(gh{3}, 'off')
            hold(gh{4}, 'off')
            hold(gh{5}, 'off')
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
            th.String = ['Geostatistical Inversion - Raytracing, Iteration ',num2str(noIter)];
            drawnow
        else
            disp(['Geostatistical Inversion -  Raytracing, Iteration ',num2str(noIter)])
        end
        [~,tomo.rays,L] = grid.raytrace(tomo.s,tomo.xi,tomo.th,data(:,1:3),data(:,4:6));
        Lx = L(:,1:nCells);
        Lz = L(:,(1+nCells):np);
    end

    if param.saveInvData == 1
        co = kron(ones(nt,1),cos(tomo.th)');
            si = kron(ones(nt,1),sin(tomo.th)');
            tt = (Lx.*co + Lz.*si).^2 + (Lz.*co - Lx.*si).^2.*kron(ones(nt,1),tomo.xi'.^2);
            tt = kron(ones(nt,1),tomo.s') .* sqrt(tt);
            tt = sum(tt,2);
            tomo.invData(noIter).res = data(:,7)-tt;
            tomo.invData(noIter).resJ = resJ;
            tomo.invData(noIter).s  = tomo.s;
            tomo.invData(noIter).xi = tomo.xi;
            tomo.invData(noIter).th = tomo.th;
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
