function tomo = inv_difference2D(panel, p, t_handle, g_handles)
%
%
%

str = get_str_locale();
if ~isempty(t_handle)
    set(t_handle, 'String',['Difference inversion - ',str.s227])
else
    disp(['Difference inversion - ',str.s227])
end

tomo.rays = {};
tomo.L = [];
tomo.invData = [];

ref_mog_no = panel.inv_res(p.ref_inv_no).param.selected_mogs;

[t1_obs,ind1] = getPanneauData(panel, p.db_file, p.type_data, ref_mog_no);
[t2_obs,ind2] = getPanneauData(panel, p.db_file, p.type_data, p.mog_no);

TxRx1 = [panel.grid.Tx(ind1,:) panel.grid.Rx(ind1,:)];
TxRx2 = [panel.grid.Tx(ind2,:) panel.grid.Rx(ind2,:)];

[i1,i2] = findCommonTxRx(TxRx1, TxRx2);

if isempty(i1)
    errordlg('Tx and Rx coordinates inconsistent between traveltime datasets', ...
        'Inconsistent data','modal')
    return
end

s1 = panel.inv_res(p.ref_inv_no).tomo.s;
L1 = panel.inv_res(p.ref_inv_no).tomo.L;
s = s1;
if p.tomo_amp==1
    no=t2_obs(i2,3);
    ind = zeros(size(no));
    for n=1:length(no)
        nn=find(no(n) == panel.inv_res(p.L_tomo_no).tomo.no_trace);
        if isempty(nn)
            errordlg('Inconsistency between traveltime and amplitude datasets', ...
                'Inconsistent data','modal')
            return
        end
        ind(n)=nn;
    end
    L = panel.inv_res(p.L_tomo_no).tomo.L(ind,:);
else
    L = L1;
end
ind = panel.inv_res(p.ref_inv_no).tomo.no_trace;

x = panel.inv_res(p.ref_inv_no).tomo.x;
z = panel.inv_res(p.ref_inv_no).tomo.z;
Tx = [panel.grid.Tx(ind,1) panel.grid.Tx(ind,3)];
Rx = [panel.grid.Rx(ind,1) panel.grid.Rx(ind,3)];

g.xmin = panel.grid.grx(1);
g.zmin = panel.grid.grz(1);
g.dx = panel.grid.grx(2)-panel.grid.grx(1);
g.dz = panel.grid.grz(2)-panel.grid.grz(1);
g.nx = length(panel.grid.grx)-1;
g.nz = length(panel.grid.grz)-1;
g.nsx = 9;
g.nsz = 9;

M=length(z);
N=length(x);

[Gx,Gz]=deriv2D(ones(size(x)),ones(size(z)),1);
Gx = sparse(Gx);
Gz = sparse(Gz);

wt=ones(M*N,1);
if ~isempty(p.ind_reservoir)
    wt(p.ind_reservoir(:)) = p.weight_reservoir;
end
WT = sparse(diag(wt));


Wm = WT'*(p.lambda*speye(length(s)) + p.mu*(Gx'*Gx) + p.eta*(Gz'*Gz))*WT;
if any(t2_obs(i2,2)>0)
    Wd = sparse(diag(1./(t2_obs(i2,2).^2)));
else
    Wd = speye(length(t2_obs(i2,2)));
end

WmWm = p.beta*(Wm'*Wm);
WdWd = Wd'*Wd;

clear Wd Wm wt Wt

for noIter=1:p.max_it
    
    if p.tomo_amp==0
        if ~isempty(t_handle)
            set(t_handle, 'String',['Traveltime inversion - ',str.s228,' ',num2str(noIter)])
            drawnow
        else
            disp(['Traveltime inversion - ',str.s228,' ',num2str(noIter)]); drawnow
        end
    else
        if ~isempty(t_handle)
            set(t_handle, 'String',['Amplitude inversion - ',str.s228,' ',num2str(noIter)])
            drawnow
        else
            disp(['Amplitude inversion - ',str.s228,' ',num2str(noIter)]); drawnow
        end
    end
    
    DeltaD = ((L*s - L1*s1) - (t2_obs(i2,1)-t1_obs(i1,1)));
    
    b = -(L'*(WdWd)*DeltaD + WmWm*(s-s1));
    A = L'*(WdWd)*L + WmWm;

    delta_s = lsqr(A,b,[],500);

    fac = p.damping*max(abs(delta_s))/max(abs(s));
    s = s+fac*delta_s;
    
    if p.tomo_amp==0

        if ~isempty(g_handles)
            imagesc(x,z,1./reshape(s,M,N),'Parent',g_handles{4})
            set(g_handles{4},'DataAspectRatio',[1 1 1],'YDir','normal')
            title(g_handles{4},['Repeat survey: velocity - iteration ',num2str(noIter)],'FontSize',14)
            xlabel(g_handles{4},str.s119)
            ylabel(g_handles{4},str.s120)
            if ~isempty(g_handles{1}), caxis(g_handles{4}, g_handles{1}), end
            colorbar('peer',g_handles{4})
            eval(['colormap(',g_handles{2},')'])
            drawnow
        end
    
        if noIter==p.max_it
            [tt,tomo.rays,L] = ttcr2d(s,g,Tx,Rx);
        else
            [tt,~,L] = ttcr2d(s,g,Tx,Rx);
        end
        
    else
        
        if ~isempty(g_handles)
            imagesc(x,z,reshape(s,M,N),'Parent',g_handles{4})
            set(g_handles{4},'DataAspectRatio',[1 1 1],'YDir','normal')
            title(g_handles{4},['Repeat survey: attenuation - iteration ',num2str(noIter)],'FontSize',14)
            xlabel(g_handles{4},str.s119)
            ylabel(g_handles{4},str.s120)
            if ~isempty(g_handles{5}), caxis(g_handles{4}, g_handles{5}), end
            colorbar('peer',g_handles{4})
            eval(['colormap(',g_handles{2},')'])
        end
        
        tt = L*s;
    end
    
    if p.saveInvData==1
        tomo.invData(noIter).res = t2_obs(i2,1)-tt;
        tomo.invData(noIter).s = s;
    end
end

load(p.db_file,'mogs')

tomo.L = L;
tomo.no_trace = mogs(p.mog_no).no_traces;
tomo.s = s;
tomo.x = x;
tomo.z = z;
tomo.corr_tt_Lant = [];
tomo.date = mogs(p.mog_no).date;

if p.saveInvData == 1
    tomo.invData(end).date = datestr(now);
end

end