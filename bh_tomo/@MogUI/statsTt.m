function statsTt(obj,varargin)

no = obj.handles.listMogs.Value;
if no>0 && no<=length(obj.mogs)
    mog = obj.mogs(no);
    done = mog.tt_done & mog.in;
    if sum(done)==0
        errordlg('Data not processed!!!')
        return
    end
    ind = mog.tt~=-1 & mog.in;
    [tt,t0] = mog.getCorrectedTravelTimes(obj.air);
    obj.handles.editDtCorr.String = num2str( mog.fac_dt );
    et = mog.et(ind);
    tt = tt(ind);
    t0 = t0(ind);
    
    hyp = sqrt( (mog.data.Tx_x(ind)-mog.data.Rx_x(ind)).^2 + ...
        (mog.data.Tx_y(ind)-mog.data.Rx_y(ind)).^2 + ...
        (mog.data.Tx_z(ind)-mog.data.Rx_z(ind)).^2 );
    dz = mog.data.Rx_z(ind)-mog.data.Tx_z(ind);
    theta = 180/pi*asin(dz./hyp);
    h_stat = figure;
    set(h_stat, 'Position',[256 71 512 620]);
    
    subplot(321)
    plot(hyp, tt, 'o')
    xlabel('Straight ray length [m]')
    ylabel(['Time',' [',mog.data.tunits,']'])
    
    subplot(322)
    plot(theta, hyp./tt, 'o')
    xlabel('Angle w/r to horizontal [deg]')
    ylabel(['Apparent velocity', '[',mog.data.cunits,'/',mog.data.tunits,']'])
    title('Velocity before correction')
    
    vapp = hyp./(tt-t0(ind));
    n = 1:length(ind);
    n=n(ind);
    ind2 = vapp<0;
    if ~isempty(n(ind2))
        disp(['velocity apparente negative: trace', num2str(n(ind2))])
    end
    
    subplot(323)
    plot(theta, hyp./(tt-t0(ind)), 'o')
    xlabel('Angle w/r to horizontal [deg]')
    ylabel(['Apparent velocity', '[',mog.data.cunits,'/',mog.data.tunits,']'])
    title('Velocity after correction (air)')
    
    subplot(324)
    plot(t0)
    xlabel('Shot number')
    ylabel(['Time',' [',mog.data.tunits,']'])
    title('t_0 drift in air')
    
    subplot(325)
    plot(hyp,et,'o')
    xlabel('Straight ray length [m]')
    ylabel('Standard deviation')
    
    subplot(326)
    plot(theta,et,'o')
    xlabel('Angle w/r to horizontal [deg]')
    ylabel('Standard deviation')
    
    suptitle(mog.name,'Interpreter','none')
    
    h_stat2 = figure;
    set(h_stat2, 'Position',[356 71 512 620]);
    drawnow
    
    vapp = hyp./tt;
    %lapp = tt./hyp;
    Tx = [mog.data.Tx_x(ind)' mog.data.Tx_y(ind)' mog.data.Tx_z(ind)'];
    Rx = [mog.data.Rx_x(ind)' mog.data.Rx_y(ind)' mog.data.Rx_z(ind)'];
    
    vmin = min(vapp);
    vmax = max(vapp);
    c=cmr;
    
    [~, a]=lsplane([Tx; Rx]);
    el = (pi-a(3))*180/pi;
    az = atan( cos(a(2))/cos(a(1)) )*180/pi;
    
    
    m = (size(c,1)-1)/(vmax-vmin);
    b = 1-vmin*m;
    p = m*vapp(1)+b;
    couleur = interp1(c,p);
    
    plot3([Tx(1,1) Rx(1,1)], [Tx(1,2) Rx(1,2)], [Tx(1,3) Rx(1,3)], 'Color',couleur)
    hold on
    for n=2:length(vapp)
        p = m*vapp(n)+b;
        couleur = interp1(c,p);
        plot3([Tx(n,1) Rx(n,1)], [Tx(n,2) Rx(n,2)], [Tx(n,3) Rx(n,3)], 'Color',couleur)
    end
    hold off
    set(gca, 'DataAspectRatio',[1 1 1],'Units','normalized')
    axis tight
    view(az,el)
    caxis([vmin vmax])
    colormap(cmr)
    colorbar()
    title([mog.name,' - apparent velocity'],'Interpreter','none')
    tlabh = get(gca,'Title');
    set(tlabh,'FontSize',12,'Position',get(tlabh,'Position') + [0 0 .4])
    
end
end
