function tomo = invDifference(model,p,hmessage,gh,gv)

if ~isempty(hmessage)
    hmessage.String = 'Difference Inversion - Starting ...';
else
    disp('Difference Inversion - Starting ...');
end

tomo.rays = {};
tomo.L = [];
tomo.invData = [];

ref_mog_no = model.inv_res(p.ref_inv_no).param.selectedMogs;

[t0_obs,ind0] = Model.getModelData(model,p.db_file,p.typeData,ref_mog_no);
[t1_obs,ind1] = Model.getModelData(model,p.db_file,p.typeData,p.mog_no);

TxRx0 = [model.grid.Tx(ind0,:) model.grid.Rx(ind0,:)];
TxRx1 = [model.grid.Tx(ind1,:) model.grid.Rx(ind1,:)];

[i0,i1] = findCommonTxRx(TxRx0, TxRx1);

if isempty(i0)
    errordlg('Tx and Rx coordinates inconsistent between traveltime datasets', ...
        'Inconsistent data','modal')
    return
end

s0 = model.inv_res(p.ref_inv_no).tomo.s;
L0 = model.inv_res(p.ref_inv_no).tomo.L;
s = s0;
if p.tomoAtt==1
    no=t1_obs(i1,3);
    ind = zeros(size(no));
    for n=1:length(no)
        nn=find(no(n) == model.inv_res(p.L_tomo_no).tomo.no_trace);
        if isempty(nn)
            errordlg('Inconsistency between traveltime and amplitude datasets', ...
                'Inconsistent data','modal')
            return
        end
        ind(n)=nn;
    end
    L = model.inv_res(p.L_tomo_no).tomo.L(ind,:);
else
    L = L0;
end
ind = model.inv_res(p.ref_inv_no).tomo.no_trace;

Tx = model.grid.Tx(ind,:);
Rx = model.grid.Rx(ind,:);

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
if any(t1_obs(i1,2)>0)
    Wd = sparse(diag(1./(t1_obs(i1,2).^2)));
else
    Wd = speye(length(t1_obs(i1,2)));
end

WmWm = p.beta*(Wm'*Wm);
WdWd = Wd'*Wd;

clear Wd Wm wt Wt

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

    DeltaD = ((L*s - L0*s0) - (t1_obs(i1,1)-t0_obs(i0,1)));
    
    b = -(L'*(WdWd)*DeltaD + WmWm*(s-s0));
    A = L'*(WdWd)*L + WmWm;

    delta_s = lsqr(A,b,[],500);

    fac = p.damping*max(abs(delta_s))/max(abs(s));
    s = s+fac*delta_s;

    if p.tomoAtt==0

        if ~isempty(gh)
            gv.plotTomo(1./s,['Repeat Survey: Velocity - Iteration ',num2str(noIter)],...
                'Distance [m]','Elevation [m]',gh{4})
            if ~isempty(gh{1}), caxis(gh{4},gh{1}), end
            colorbar('peer', gh{4})
            
            eval(['colormap(',gh{2},')'])
            drawnow
        end

        if ~isempty(hmessage)
            hmessage.String = ['Traveltime inversion - Raytracing, Iteration ',num2str(noIter)];
            drawnow
        else
            disp(['Traveltime inversion - Raytracing, Iteration ',num2str(noIter)]); drawnow
        end
        if noIter==p.max_it
            [tt,tomo.rays,L] = model.grid.raytrace(s,Tx,Rx);
        else
            [tt,~,L] = model.grid.raytrace(s,Tx,Rx);
        end
            
    else
            
        if ~isempty(gh)
            gv.plotTomo(s,['Repeat Survey: Attenuation - Iteration ',num2str(noIter)],...
                'Distance [m]','Elevation [m]',gh{4})
            if ~isempty(gh{1}), caxis(gh{4},gh{1}), end
            colorbar('peer', gh{4})
            
            eval(['colormap(',gh{2},')'])
            drawnow
        end
        
        tt = L*s;
    end
            
    if p.saveInvData==1
        tomo.invData(noIter).res = t1_obs(i1,1)-tt;
        tomo.invData(noIter).s = s;
    end
end

load(p.db_file,'mogs')

tomo.L = L;
tomo.no_trace = mogs(p.mog_no).no_traces;
tomo.s = s;
tomo.s0 = s0;
tomo.date = mogs(p.mog_no).date;

if p.saveInvData == 1
    tomo.invData(end).date = datestr(now);
end

end         