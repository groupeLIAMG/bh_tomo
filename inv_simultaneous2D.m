function tomo = inv_simultaneous2D(panel, p, t_handle, g_handles)
%
%

str = get_str_locale();
if ~isempty(t_handle)
    set(t_handle, 'String',['Simultaneous inversion - ',str.s227])
else
    disp(['Simultaneous inversion - ',str.s227])
end

tomo.rays = {};
tomo.L = [];
tomo.invData = [];

[t0_obs,ind0] = getPanneauData(panel, p.db_file, p.type_data, p.mog_no0);
[t1_obs,ind1] = getPanneauData(panel, p.db_file, p.type_data, p.mog_no1);
Tx0=panel.grid.Tx(ind0,[1 3]);
Rx0=panel.grid.Rx(ind0,[1 3]);
Tx1=panel.grid.Tx(ind1,[1 3]);
Rx1=panel.grid.Rx(ind1,[1 3]);

if p.tomo_amp==1
    no=t0_obs(:,3);
    ind = zeros(size(no));
    for n=1:length(no)
        nn=find(no(n) == panel.inv_res(p.L_tomo_no).tomo.no_trace0);
        if isempty(nn)
            errordlg('Inconsistency between baseline traveltime and amplitude datasets', ...
                'Inconsistent data','modal')
            return
        end
        ind(n)=nn;
    end
    L0 = panel.inv_res(p.L_tomo_no).tomo.L0(ind,:);
    
    no=t1_obs(:,3);
    ind = zeros(size(no));
    for n=1:length(no)
        nn=find(no(n) == panel.inv_res(p.L_tomo_no).tomo.no_trace);
        if isempty(nn)
            errordlg('Inconsistency between repeat traveltime and amplitude datasets', ...
                'Inconsistent data','modal')
            return
        end
        ind(n)=nn;
    end
    L1 = panel.inv_res(p.L_tomo_no).tomo.L(ind,:);
    
else
    L0=Lsr2d(Tx0,Rx0,panel.grid.grx,panel.grid.grz);
    L1=Lsr2d(Tx1,Rx1,panel.grid.grx,panel.grid.grz);
end


x = panel.inv_res(p.mog_no0).tomo.x;
z = panel.inv_res(p.mog_no0).tomo.z;

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

Wm = WT'*(p.lambda*speye(M*N) + p.mu*(Gx'*Gx) + p.eta*(Gz'*Gz))*WT;
WmWm = Wm'*Wm;
WdiffWdiff=WmWm;

WmWm = p.beta*WmWm;
WdiffWdiff = p.alpha*WdiffWdiff;

Wd0 = sparse(diag(1./(t0_obs(:,2).^2)));
Wd1 = sparse(diag(1./(t1_obs(:,2).^2)));
Wd0Wd0 = Wd0'*Wd0;
Wd1Wd1 = Wd1'*Wd1;

clear Wm Wd0 Wd1

s_ref = panel.inv_res(p.ref_inv_no).tomo.s;
s0 = s_ref;
s1 = s_ref;

tt0=L0*s0;
tt1=L1*s1;

ds_reff = 0;

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

    if p.tomo_amp==0

        if ~isempty(g_handles)

            imagesc(x,z,reshape(1./s0,length(z),length(x)),'Parent',g_handles{3})
            set(g_handles{3},'DataAspectRatio',[1 1 1],'YDir','normal')
            title(g_handles{3},['Velocity 0 - iteration ',num2str(noIter)],'FontSize',14)
            xlabel(g_handles{3},str.s119)
            ylabel(g_handles{3},str.s120)
            if ~isempty(g_handles{1}), caxis(g_handles{3}, g_handles{1}), end
            colorbar('peer',g_handles{3})
            eval(['colormap(',g_handles{2},')'])

            imagesc(x,z,reshape(1./s1,length(z),length(x)),'Parent',g_handles{4})
            set(g_handles{4},'DataAspectRatio',[1 1 1],'YDir','normal')
            title(g_handles{4},['Velocity 1 - iteration ',num2str(noIter)],'FontSize',14)
            xlabel(g_handles{4},str.s119)
            ylabel(g_handles{4},str.s120)
            if ~isempty(g_handles{1}), caxis(g_handles{4}, g_handles{1}), end
            colorbar('peer',g_handles{4})
            eval(['colormap(',g_handles{2},')'])
    
            drawnow
        end

        [tt0,~,L0] = ttcr2d(s0,g,Tx0,Rx0);
        if noIter==p.max_it
            [tt1,tomo.rays,L1] = ttcr2d(s1,g,Tx1,Rx1);
        else
            [tt1,~,L1] = ttcr2d(s1,g,Tx1,Rx1);
        end
        
    else
        
        if ~isempty(g_handles)

            imagesc(x,z,reshape(s0,length(z),length(x)),'Parent',g_handles{3})
            set(g_handles{3},'DataAspectRatio',[1 1 1],'YDir','normal')
            title(g_handles{3},['Attenuation 0 - iteration ',num2str(noIter)],'FontSize',14)
            xlabel(g_handles{3},str.s119)
            ylabel(g_handles{3},str.s120)
            if ~isempty(g_handles{1}), caxis(g_handles{3}, g_handles{1}), end
            colorbar('peer',g_handles{3})
            eval(['colormap(',g_handles{2},')'])

            imagesc(x,z,reshape(s1,length(z),length(x)),'Parent',g_handles{4})
            set(g_handles{4},'DataAspectRatio',[1 1 1],'YDir','normal')
            title(g_handles{4},['Attenuation 1 - iteration ',num2str(noIter)],'FontSize',14)
            xlabel(g_handles{4},str.s119)
            ylabel(g_handles{4},str.s120)
            if ~isempty(g_handles{1}), caxis(g_handles{4}, g_handles{1}), end
            colorbar('peer',g_handles{4})
            eval(['colormap(',g_handles{2},')'])
    
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
tomo.s = s;
tomo.x = x;
tomo.z = z;
tomo.corr_tt_Lant = [];
tomo.date = mogs(p.mog_no1).date;

if p.saveInvData == 1
    tomo.invData(end).date = datestr(now);
end

end




