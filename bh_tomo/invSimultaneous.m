function tomo = invSimultaneous(model,p,hmessage,gh,gv)

if ~isempty(hmessage)
    hmessage.String = 'Simultaneous Inversion - Starting ...';
    drawnow
else
    disp('Simultaneous Inversion - Starting ...');
end

tomo.rays = {};
tomo.L = [];
tomo.invData = [];

[t0_obs,ind0] = Model.getModelData(model,p.db_file,p.typeData,p.mog_no0);
[t1_obs,ind1] = Model.getModelData(model,p.db_file,p.typeData,p.mog_no1);
Tx0=model.grid.Tx(ind0,:);
Rx0=model.grid.Rx(ind0,:);
Tx1=model.grid.Tx(ind1,:);
Rx1=model.grid.Rx(ind1,:);

if p.tomoAtt == 1
    no=t0_obs(:,3);
    ind = zeros(size(no));
    for n=1:length(no)
        nn=find(no(n) == model.inv_res(p.L_tomo_no).tomo.no_trace0);
        if isempty(nn)
            errordlg('Inconsistency between baseline traveltime and amplitude datasets', ...
                'Inconsistent data','modal')
            return
        end
        ind(n)=nn;
    end
    L0 = model.inv_res(p.L_tomo_no).tomo.L0(ind,:);
    
    no=t1_obs(:,3);
    ind = zeros(size(no));
    for n=1:length(no)
        nn=find(no(n) == model.inv_res(p.L_tomo_no).tomo.no_trace);
        if isempty(nn)
            errordlg('Inconsistency between repeat traveltime and amplitude datasets', ...
                'Inconsistent data','modal')
            return
        end
        ind(n)=nn;
    end
    L1 = model.inv_res(p.L_tomo_no).tomo.L(ind,:);
else
    L0 = model.grid.getForwardStraightRays(ind0);
    L1 = model.grid.getForwardStraightRays(ind1);
end

tomo.x = 0.5*(model.grid.grx(1:end-1)+model.grid.grx(2:end));
tomo.z = 0.5*(model.grid.grz(1:end-1)+model.grid.grz(2:end));
if ~isempty(model.grid.gry)
    tomo.y = 0.5*(model.grid.gry(1:end-1)+model.grid.gry(2:end));
else
    tomo.y = [];
end


[Dx,Dy,Dz] = model.grid.derivative(p.Dorder);
if isempty(Dy)
    Dy=0;
end

nc = model.grid.getNumberOfCells();

wt=ones(nc,1);
if ~isempty(p.ind_reservoir)
    wt(p.ind_reservoir(:)) = p.weight_reservoir;
end
WT = sparse(diag(wt));

Wm = WT'*(p.lambda*speye(nc) + p.mux*(Dx'*Dx) + p.muy*(Dy'*Dy) + p.muz*(Dz'*Dz))*WT;

WmWm = Wm'*Wm;
WdiffWdiff=WmWm;

WmWm = p.beta*WmWm;
WdiffWdiff = p.alpha*WdiffWdiff;
clear Wm

Wd0 = sparse(diag(1./(t0_obs(:,2).^2)));
Wd1 = sparse(diag(1./(t1_obs(:,2).^2)));
Wd0Wd0 = Wd0'*Wd0;
Wd1Wd1 = Wd1'*Wd1;

clear Wd0 Wd1

s_ref = model.inv_res(p.ref_inv_no).tomo.s;
s0 = s_ref;
s1 = s_ref;

tt0=L0*s0;
tt1=L1*s1;

ds_reff = 0;

for noIter=1:p.max_it

    if p.tomoAtt==0
        if ~isempty(hmessage)
            hmessage.String = ['Traveltime inversion - Solving System, Iteration ',num2str(noIter)];
            drawnow
        else
            disp(['Traveltime inversion - Solving System, Iteration ',num2str(noIter)]); drawnow
        end
    else
        if ~isempty(hmessage)
            hmessage.String = ['Amplitude inversion - Solving System, Iteration ',num2str(noIter)];
            drawnow
        else
            disp(['Amplitude inversion - Solving System, Iteration ',num2str(noIter)]); drawnow
        end
    end

    A = L0'*(Wd0Wd0)*L0 + WmWm + WdiffWdiff;
    b = -(L0'*(Wd0Wd0)*(tt0 - t0_obs(:,1)) + WmWm*(s0-s_ref)) + ...
        WdiffWdiff*(s0-s1 - ds_reff);

    delta_s0 = lsqr(A,b,[],500);
    fac = p.damping*max(abs(delta_s0))/max(abs(s0));
    s0 = s0+fac*delta_s0;

    A = L1'*(Wd1Wd1)*L1 + WmWm + WdiffWdiff;
    b = -(L1'*(Wd1Wd1)*(tt1 - t1_obs(:,1)) + WmWm*(s1-s_ref)) + ...
        WdiffWdiff*(s0-s1 - ds_reff);

    delta_s1 = lsqr(A,b,[],500);
    fac = p.damping*max(abs(delta_s1))/max(abs(s1));
    s1 = s1+fac*delta_s1;

    if p.tomoAtt==0

        if ~isempty(gh)
            gv.plotTomo(1./s0,['Velocity 0 - iteration ',num2str(noIter)],...
                'Distance [m]','Elevation [m]',gh{3})
            if ~isempty(gh{1}), caxis(gh{3},gh{1}), end
            colorbar('peer', gh{3})
            
            gv.plotTomo(1./s1,['Velocity 1 - iteration ',num2str(noIter)],...
                'Distance [m]','Elevation [m]',gh{4})
            if ~isempty(gh{1}), caxis(gh{4},gh{1}), end
            colorbar('peer', gh{4})
            
            colormap( gh{3}, gh{2})
            colormap( gh{4}, gh{2})
            drawnow
        end
    
        if ~isempty(hmessage)
            hmessage.String = ['Traveltime inversion - Raytracing, Iteration ',num2str(noIter)];
            drawnow
        else
            disp(['Traveltime inversion - Raytracing, Iteration ',num2str(noIter)]); drawnow
        end
        
        [tt0,~,L0] = model.grid.raytrace(s0,Tx0,Rx0);
        if noIter==p.max_it
            [tt1,tomo.rays,L1] = model.grid.raytrace(s1,Tx1,Rx1);
        else
            [tt1,~,L1] = model.grid.raytrace(s1,Tx1,Rx1);
        end
        
    else
        
        if ~isempty(gh)
            gv.plotTomo(s0,['Attenuation 0 - iteration ',num2str(noIter)],...
                'Distance [m]','Elevation [m]',gh{3})
            if ~isempty(gh{1}), caxis(gh{3},gh{1}), end
            colorbar('peer', gh{3})
            
            gv.plotTomo(s1,['Attenuation 1 - iteration ',num2str(noIter)],...
                'Distance [m]','Elevation [m]',gh{4})
            if ~isempty(gh{1}), caxis(gh{4},gh{1}), end
            colorbar('peer', gh{4})
            
            colormap( gh{3}, gh{2})
            colormap( gh{4}, gh{2})
            drawnow
        end
        
        tt0 = L0*s0;
        tt1 = L1*s1;
    end
    
    if p.saveInvData==1
        tomo.invData(noIter).res0 = t0_obs(:,1)-tt0;
        tomo.invData(noIter).s0 = s0;
        tomo.invData(noIter).res = t1_obs(:,1)-tt1;
        tomo.invData(noIter).s = s1;
    end
end

load(p.db_file,'mogs')

tomo.s0 = s0;
tomo.s = s1;
tomo.L0 = L0;
tomo.L = L1;
tomo.no_trace0 = mogs(p.mog_no0).no_traces;
tomo.no_trace = mogs(p.mog_no1).no_traces;
tomo.date = mogs(p.mog_no1).date;

if p.saveInvData == 1
    tomo.invData(end).date = datestr(now);
end

end
