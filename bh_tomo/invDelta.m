function tomo = invDelta(model,p,hmessage,gh,gv)

if ~isempty(hmessage)
    hmessage.String = [char(916),'t / ',char(916),char(964),' Inversion - Starting ...'];
else
    disp([char(916),'t / ',char(916),char(964),' Inversion - Starting ...']);
end

tomo.rays = {};
tomo.L = [];
tomo.invData = [];

[data,~,delta] = Model.getModelData(model,p.db_file,p.typeData,p.mog_no);

if delta == false
    errordlg('Data should be of type delta')
    return
end

s0 = model.inv_res(p.ref_inv_no).tomo.s;
L = model.inv_res(p.ref_inv_no).tomo.L;

[i0, i1] = findCommonTraces(data(:,3), model.inv_res(p.ref_inv_no).tomo.no_trace);

dt = data(i0,1);
L = L(i1,:);
    
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

if any(data(:,2)>0)
    Wd = sparse(diag(1./(data(i0,2).^2)));
else
    Wd = speye(length(data(i0,2)));
end

Wm = p.mux*(Dx'*Dx) + p.muy*(Dy'*Dy) + p.muz*(Dz'*Dz);

WmWm = p.beta*(Wm'*Wm);
WdWd = Wd'*Wd;
A = L'*(WdWd)*L + WmWm;

clear Wd Wm

s = zeros(size(s0));
for noIter=1:p.max_it

    if ~isempty(hmessage)
        hmessage.String = [char(916),'t / ',char(916),char(964),' Inversion - Solving System, Iteration ',num2str(noIter)];
        drawnow
    else
        disp([char(916),'t / ',char(916),char(964),' Inversion - Solving System, Iteration ',num2str(noIter)])
    end

    b = -(L'*(WdWd)*dt + WmWm*s);

    x = lsqr(A,b,[],500);
    
    s = s + x*p.damping;
    
    if p.tomoAtt==0

        if ~isempty(gh)
            gv.plotTomo(1./(s0+s),['Repeat Survey: Velocity - Iteration ',num2str(noIter)],...
                'Distance [m]','Elevation [m]',gh{4})
            if ~isempty(gh{1}), caxis(gh{4},gh{1}), end
            colorbar('peer', gh{4})
            
            colormap( gh{4}, gh{2})
            drawnow
        end

    else
        
        if ~isempty(gh)
            gv.plotTomo(s0+s,['Repeat Survey: Attenuation - Iteration ',num2str(noIter)],...
                'Distance [m]','Elevation [m]',gh{4})
            if ~isempty(gh{1}), caxis(gh{4},gh{1}), end
            colorbar('peer', gh{4})
            
            colormap( gh{4}, gh{2})
            drawnow
        end
    end
    
    if p.saveInvData == 1
        tomo.invData(noIter).res = dt-L*s;
        tomo.invData(noIter).s = s0+s;
    end

end

load(p.db_file,'mogs')

tomo.L = L;
tomo.no_trace = mogs(p.mog_no).no_traces;
tomo.s = s0 + s;
tomo.s0 = s0;
tomo.date = mogs(p.mog_no).date;

if p.saveInvData == 1
    tomo.invData(end).date = datestr(now);
end

end