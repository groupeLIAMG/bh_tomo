function verify_struct_fields(panels,mogs,boreholes)
%created 2012-11-14  YH
msg = {};
%panels
nb_panels = length(panels);
for i = 1:nb_panels
    if ~isfield(panels(i),'name')
        msg{end+1} = 'panels.name';
    end
    
    if ~isfield(panels(i),'mogs')
        msg{end+1} = 'panels.mogs';
    end
    
    if ~isfield(panels(i),'boreholes')
        msg{end+1} = 'panels.boreholes';
    end
    
    if ~isfield(panels(i),'type')
        msg{end+1} = 'panels.type';
    end
    
    if ~isfield(panels(i),'grid')
        msg{end+1} = 'panels.grid';
    end
    
    if ~isfield(panels(i),'tt_covar')
        msg{end+1} = 'panels.tt_covar';
    end
    
    if ~isfield(panels(i),'inv_res')
        msg{end+1} = 'panels.inv_res';
    end
    
    if ~isfield(panels(i),'amp_covar')
        msg{end+1} = 'panels.amp_covar';
    end
    
    %panels(i).grid
    if ~isfield(panels(i).grid,'grx')
        msg{end+1} = ['panels(' num2str(i) ').grid.grx'];
    end
    
    if ~isfield(panels(i).grid,'grz')
        msg{end+1} = ['panels(' num2str(i) ').grid.grz'];
    end
    
    if ~isfield(panels(i).grid,'cont')
        msg{end+1} = ['panels(' num2str(i) ').grid.cont'];
    end
    
    if ~isfield(panels(i).grid,'Tx')
        msg{end+1} = ['panels(' num2str(i) ').grid.Tx'];
    end
    
    if ~isfield(panels(i).grid,'Rx')
        msg{end+1} = ['panels(' num2str(i) ').grid.Rx'];
    end
    
    if ~isfield(panels(i).grid,'TxCosDir')
        msg{end+1} = ['panels(' num2str(i) ').grid.TxCosDir'];
    end
    
    if ~isfield(panels(i).grid,'RxCosDir')
        msg{end+1} = ['panels(' num2str(i) ').grid.RxCosDir'];
    end
    
    if ~isfield(panels(i).grid,'x0')
        msg{end+1} = ['panels(' num2str(i) ').grid.x0'];
    end
    
    if ~isfield(panels(i).grid,'borehole_x0')
        msg{end+1} = ['panels(' num2str(i) ').grid.forage_x0'];
    end
    
    if ~isfield(panels(i).grid,'bord')
        msg{end+1} = ['panels(' num2str(i) ').grid.bord'];
    end
    
    if ~isfield(panels(i).grid,'flip')
        msg{end+1} = ['panels(' num2str(i) ').grid.flip'];
    end
    
    if ~isfield(panels(i).grid,'Tx_Z_water')
        msg{end+1} = ['panels(' num2str(i) ').grid.Tx_Z_water'];
    end
    
    if ~isfield(panels(i).grid,'Rx_Z_water')
        msg{end+1} = ['panels(' num2str(i) ').grid.Rx_Z_water'];
    end
    
    if ~isfield(panels(i).grid,'in')
        msg{end+1} = ['panels(' num2str(i) ').grid.in'];
    end
    
    %panels(i).inv_res
    if isfield(panels(i),'inv_res')
        %if isempty(panels(i).inv_res)
        
        nb_inv = length(panels(i).inv_res);
        for j = 1:nb_inv
            if ~isfield(panels(i).inv_res(j),'name')
                msg{end+1} = ['panels(' num2str(i) ').inv_res(' num2str(j) ').name'];
            end
            if ~isfield(panels(i).inv_res(j),'tomo')
                msg{end+1} = ['panels(' num2str(i) ').inv_res(' num2str(j) ').tomo'];
            end
            if ~isfield(panels(i).inv_res(j),'param')
                msg{end+1} = ['panels(' num2str(i) ').inv_res(' num2str(j) ').param'];
            end
            if ~isfield(panels(i).inv_res(j),'covar')
                msg{end+1} = ['panels(' num2str(i) ').inv_res(' num2str(j) ').covar'];
            end
        end
        
        %panels(i).inv_res(i).tomo
        if isfield(panels(i).inv_res,'tomo')
            nb_inv = length(panels(i).inv_res);
            for j = 1:nb_inv
                if ~isfield(panels(i).inv_res(j).tomo,'rays')
                    msg{end+1} = ['panels(' num2str(i) ').inv_res(' num2str(j) ').tomo.rays'];
                end
                if ~isfield(panels(i).inv_res(j).tomo,'L')
                    msg{end+1} = ['panels(' num2str(i) ').inv_res(' num2str(j) ').tomo.L'];
                end
                if ~isfield(panels(i).inv_res(j).tomo,'invData')
                    msg{end+1} = ['panels(' num2str(i) ').inv_res(' num2str(j) ').tomo.invData'];
                end
                if ~isfield(panels(i).inv_res(j).tomo,'no_trace')
                    msg{end+1} = ['panels(' num2str(i) ').inv_res(' num2str(j) ').tomo.no_trace'];
                end
                if ~isfield(panels(i).inv_res(j).tomo,'s')
                    msg{end+1} = ['panels(' num2str(i) ').inv_res(' num2str(j) ').tomo.s'];
                end
                if ~isfield(panels(i).inv_res(j).tomo,'x')
                    msg{end+1} = ['panels(' num2str(i) ').inv_res(' num2str(j) ').tomo.x'];
                end
                if ~isfield(panels(i).inv_res(j).tomo,'z')
                    msg{end+1} = ['panels(' num2str(i) ').inv_res(' num2str(j) ').tomo.z'];
                end
                if ~isfield(panels(i).inv_res(j).tomo,'corr_tt_Lant')
                    msg{end+1} = ['panels(' num2str(i) ').inv_res(' num2str(j) ').tomo.corr_tt_Lant'];
                end
                if ~isfield(panels(i).inv_res(j).tomo,'date')
                    msg{end+1} = ['panels(' num2str(i) ').inv_res(' num2str(j) ').tomo.date'];
                end
            end
        end
        
        % %panels(i).inv_res(i).param
        if isfield(panels(i).inv_res,'param')
            nb_inv = length(panels(i).inv_res);
            for j = 1:nb_inv
                if ~isfield(panels(i).inv_res(j).param,'model_initial')
                    msg{end+1} = ['panels(' num2str(i) ').inv_res(' num2str(j) ').param.model_initial'];
                end
                if ~isfield(panels(i).inv_res(j).param,'tomo_amp')
                    msg{end+1} = ['panels(' num2str(i) ').inv_res(' num2str(j) ').param.tomo_amp'];
                end
                if ~isfield(panels(i).inv_res(j).param,'dv_max')
                    msg{end+1} = ['panels(' num2str(i) ').inv_res(' num2str(j) ').param.dv_max'];
                end
                if ~isfield(panels(i).inv_res(j).param,'nbreitrd')
                    msg{end+1} = ['panels(' num2str(i) ').inv_res(' num2str(j) ').param.nbreitrd'];
                end
                if ~isfield(panels(i).inv_res(j).param,'nbreitrc')
                    msg{end+1} = ['panels(' num2str(i) ').inv_res(' num2str(j) ').param.nbreitrc'];
                end
                if ~isfield(panels(i).inv_res(j).param,'nbreitra')
                    msg{end+1} = ['panels(' num2str(i) ').inv_res(' num2str(j) ').param.nbreitra'];
                end
                if ~isfield(panels(i).inv_res(j).param,'it1_rays_horiz')
                    msg{end+1} = ['panels(' num2str(i) ').inv_res(' num2str(j) ').param.it1_rays_horiz'];
                end
                if ~isfield(panels(i).inv_res(j).param,'plot_xi')
                    msg{end+1} = ['panels(' num2str(i) ').inv_res(' num2str(j) ').param.plot_xi'];
                end
                if ~isfield(panels(i).inv_res(j).param,'radar')
                    msg{end+1} = ['panels(' num2str(i) ').inv_res(' num2str(j) ').param.radar'];
                end
                if ~isfield(panels(i).inv_res(j).param,'pixnonvisite')
                    msg{end+1} = ['panels(' num2str(i) ').inv_res(' num2str(j) ').param.pixnonvisite'];
                end
                if ~isfield(panels(i).inv_res(j).param,'saveInvData')
                    msg{end+1} = ['panels(' num2str(i) ').inv_res(' num2str(j) ').param.saveInvData'];
                end
                if ~isfield(panels(i).inv_res(j).param,'use_cont')
                    msg{end+1} = ['panels(' num2str(i) ').inv_res(' num2str(j) ').param.use_cont'];
                end
                if ~isfield(panels(i).inv_res(j).param,'ival')
                    msg{end+1} = ['panels(' num2str(i) ').inv_res(' num2str(j) ').param.ival'];
                end
                if ~isfield(panels(i).inv_res(j).param,'mode')
                    msg{end+1} = ['panels(' num2str(i) ').inv_res(' num2str(j) ').param.mode'];
                end
                if ~isfield(panels(i).inv_res(j).param,'nbresimu')
                    msg{end+1} = ['panels(' num2str(i) ').inv_res(' num2str(j) ').param.nbresimu'];
                end
            end
        end
        %panels(i).inv_res(i).covar
%         if isfield(panels(i).inv_res,'covar')
%             nb_inv = length(panels(i).inv_res);
%             for j = 1:nb_inv
%                 if ~isfield(panels(i).inv_res(j).covar,'use_c0')
%                     msg{end+1} = ['panels(' num2str(i) ').inv_res(' num2str(j) ').covar.use_c0'];
%                 end
%                 if ~isfield(panels(i).inv_res(j).covar,'model')
%                     msg{end+1} = ['panels(' num2str(i) ').inv_res(' num2str(j) ').covar.model'];
%                 end
%                 if ~isfield(panels(i).inv_res(j).covar,'c')
%                     msg{end+1} = ['panels(' num2str(i) ').inv_res(' num2str(j) ').covar.c'];
%                 end
%                 if ~isfield(panels(i).inv_res(j).covar,'nugget_t')
%                     msg{end+1} = ['panels(' num2str(i) ').inv_res(' num2str(j) ').covar.nugget_t'];
%                 end
%                 if ~isfield(panels(i).inv_res(j).covar,'nugget_l')
%                     msg{end+1} = ['panels(' num2str(i) ').inv_res(' num2str(j) ').covar.nugget_l'];
%                 end
%                 if ~isfield(panels(i).inv_res(j).covar,'aniso')
%                     msg{end+1} = ['panels(' num2str(i) ').inv_res(' num2str(j) ').covar.aniso'];
%                 end
%             end
%         end
        
    end
    
       %panels(i).tt_covar
%     if isfield(panels(i),'tt_covar')
%         
%         if ~isfield(panels(i).tt_covar,'model')
%             msg{end+1} = ['panels(' num2str(i) ').tt_covar.model'];
%         end
%         
%         if ~isfield(panels(i).tt_covar,'c')
%             msg{end+1} = ['panels(' num2str(i) ').tt_covar.c'];
%         end
%         
%         if ~isfield(panels(i).tt_covar,'nugget_t')
%             msg{end+1} = ['panels(' num2str(i) ').tt_covar.nugget_t'];
%         end
%         
%         if ~isfield(panels(i).tt_covar,'nugget_l')
%             msg{end+1} = ['panels(' num2str(i) ').tt_covar.nugget_l'];
%         end
%         
%         if ~isfield(panels(i).tt_covar,'use_c0')
%             msg{end+1} = ['panels(' num2str(i) ').tt_covar.use_c0'];
%         end
%         
%         if ~isfield(panels(i).tt_covar,'aniso')
%             msg{end+1} = ['panels(' num2str(i) ').tt_covar.aniso'];
%         end
%         
%         if ~isfield(panels(i).tt_covar,'model_xi')
%             msg{end+1} = ['panels(' num2str(i) ').tt_covar.model_xi'];
%         end
%         
%         if ~isfield(panels(i).tt_covar,'c_xi')
%             msg{end+1} = ['panels(' num2str(i) ').tt_covar.c_xi'];
%         end
%         
%         if ~isfield(panels(i).tt_covar,'model_th')
%             msg{end+1} = ['panels(' num2str(i) ').tt_covar.model_th'];
%         end
%         
%         if ~isfield(panels(i).tt_covar,'c_th')
%             msg{end+1} = ['panels(' num2str(i) ').tt_covar.c_th'];
%         end
%         
%         if ~isfield(panels(i).tt_covar,'nugget_xi')
%             msg{end+1} = ['panels(' num2str(i) ').tt_covar.nugget_xi'];
%         end
%         
%         if ~isfield(panels(i).tt_covar,'nugget_th')
%             msg{end+1} = ['panels(' num2str(i) ').tt_covar.nugget_th'];
%         end
%     end
    
    %panels(i).amp_covar
%     if isfield(panels(i),'amp_covar')
%         for i = 1:nb_panels
%             if ~isfield(panels(i).amp_covar,'model')
%                 msg{end+1} = ['panels(' num2str(i) ').amp_covar.model'];
%             end
%             
%             if ~isfield(panels(i).amp_covar,'c')
%                 msg{end+1} = ['panels(' num2str(i) ').amp_covar.c'];
%             end
%             
%             if ~isfield(panels(i).amp_covar,'nugget_t')
%                 msg{end+1} = ['panels(' num2str(i) ').amp_covar.nugget_t'];
%             end
%             
%             if ~isfield(panels(i).amp_covar,'nugget_l')
%                 msg{end+1} = ['panels(' num2str(i) ').amp_covar.nugget_l'];
%             end
%             
%             if ~isfield(panels(i).amp_covar,'use_c0')
%                 msg{end+1} = ['panels(' num2str(i) ').amp_covar.use_c0'];
%             end
%             
%             if ~isfield(panels(i).amp_covar,'aniso')
%                 msg{end+1} = ['panels(' num2str(i) ').amp_covar.aniso'];
%             end
%         end
%     end
    
end




%mogs
nb_mogs = length(mogs);
for i = 1:nb_mogs
    if ~isfield(mogs(i),'name')
        msg{end+1} = ['mogs(' num2str(i) ').name'];
    end
end
for i = 1:nb_mogs
    if ~isfield(mogs(i),'data')
        msg{end+1} = ['mogs(' num2str(i) ').data'];
    end
end
for i = 1:nb_mogs
    if ~isfield(mogs(i),'av')
        msg{end+1} = ['mogs(' num2str(i) ').av'];
    end
end
for i = 1:nb_mogs
    if ~isfield(mogs(i),'ap')
        msg{end+1} = ['mogs(' num2str(i) ').ap'];
    end
end
for i = 1:nb_mogs
    if ~isfield(mogs(i),'Tx')
        msg{end+1} = ['mogs(' num2str(i) ').Tx'];
    end
end
for i = 1:nb_mogs
    if ~isfield(mogs(i),'Rx')
        msg{end+1} = ['mogs(' num2str(i) ').Rx'];
    end
end
for i = 1:nb_mogs
    if ~isfield(mogs(i),'tt')
        msg{end+1} = ['mogs(' num2str(i) ').tt'];
    end
end
for i = 1:nb_mogs
    if ~isfield(mogs(i),'et')
        msg{end+1} = ['mogs(' num2str(i) ').et'];
    end
end
for i = 1:nb_mogs
    if ~isfield(mogs(i),'tt_done')
        msg{end+1} = ['mogs(' num2str(i) ').tt_done'];
    end
end
for i = 1:nb_mogs
    if ~isfield(mogs(i),'amp_tmin')
        msg{end+1} = ['mogs(' num2str(i) ').amp_tmin'];
    end
end
for i = 1:nb_mogs
    if ~isfield(mogs(i),'amp_tmax')
        msg{end+1} = ['mogs(' num2str(i) ').amp_tmax'];
    end
end
for i = 1:nb_mogs
    if ~isfield(mogs(i),'amp_done')
        msg{end+1} = ['mogs(' num2str(i) ').amp_done'];
    end
end
for i = 1:nb_mogs
    if ~isfield(mogs(i),'App')
        msg{end+1} = ['mogs(' num2str(i) ').App'];
    end
end
for i = 1:nb_mogs
    if ~isfield(mogs(i),'fcentroid')
        msg{end+1} = ['mogs(' num2str(i) ').fcentroid'];
    end
end
for i = 1:nb_mogs
    if ~isfield(mogs(i),'scentroid')
        msg{end+1} = ['mogs(' num2str(i) ').scentroid'];
    end
end
for i = 1:nb_mogs
    if ~isfield(mogs(i),'tauApp')
        msg{end+1} = ['mogs(' num2str(i) ').tauApp'];
    end
end
for i = 1:nb_mogs
    if ~isfield(mogs(i),'tauApp_et')
        msg{end+1} = ['mogs(' num2str(i) ').tauApp_et'];
    end
end
for i = 1:nb_mogs
    if ~isfield(mogs(i),'tauFce')
        msg{end+1} = ['mogs(' num2str(i) ').tauFce'];
    end
end
for i = 1:nb_mogs
    if ~isfield(mogs(i),'tauFce_et')
        msg{end+1} = ['mogs(' num2str(i) ').tauFce_et'];
    end
end
for i = 1:nb_mogs
    if ~isfield(mogs(i),'tauHyb')
        msg{end+1} = ['mogs(' num2str(i) ').tauHyb'];
    end
end
for i = 1:nb_mogs
    if ~isfield(mogs(i),'tauHyb_et')
        msg{end+1} = ['mogs(' num2str(i) ').tauHyb_et'];
    end
end
for i = 1:nb_mogs
    if ~isfield(mogs(i),'tau_params')
        msg{end+1} = ['mogs(' num2str(i) ').tau_params'];
    end
end
for i = 1:nb_mogs
    if ~isfield(mogs(i),'fw')
        msg{end+1} = ['mogs(' num2str(i) ').fw'];
    end
end
for i = 1:nb_mogs
    if ~isfield(mogs(i),'f_et')
        msg{end+1} = ['mogs(' num2str(i) ').f_et'];
    end
end
for i = 1:nb_mogs
    if ~isfield(mogs(i),'amp_name_Ldc')
        msg{end+1} = ['mogs(' num2str(i) ').amp_name_Ldc'];
    end
end
for i = 1:nb_mogs
    if ~isfield(mogs(i),'type')
        msg{end+1} = ['mogs(' num2str(i) ').type'];
    end
end
for i = 1:nb_mogs
    if ~isfield(mogs(i),'Tx_z_orig')
        msg{end+1} = ['mogs(' num2str(i) ').Tx_z_orig'];
    end
end
for i = 1:nb_mogs
    if ~isfield(mogs(i),'Rx_z_orig')
        msg{end+1} = ['mogs(' num2str(i) ').Rx_z_orig'];
    end
end
for i = 1:nb_mogs
    if ~isfield(mogs(i),'fac_dt')
        msg{end+1} = ['mogs(' num2str(i) ').fac_dt'];
    end
end
for i = 1:nb_mogs
    if ~isfield(mogs(i),'in')
        msg{end+1} = ['mogs(' num2str(i) ').in'];
    end
end
for i = 1:nb_mogs
    if ~isfield(mogs(i),'no_traces')
        msg{end+1} = ['mogs(' num2str(i) ').no_traces'];
    end
end
for i = 1:nb_mogs
    if ~isfield(mogs(i),'TxCosDir')
        msg{end+1} = ['mogs(' num2str(i) ').TxCosDir'];
    end
end
for i = 1:nb_mogs
    if ~isfield(mogs(i),'RxCosDir')
        msg{end+1} = ['mogs(' num2str(i) ').RxCosDir'];
    end
end
for i = 1:nb_mogs
    if ~isfield(mogs(i),'user_fac_dt')
        msg{end+1} = ['mogs(' num2str(i) ').user_fac_dt'];
    end
end
for i = 1:nb_mogs
    if ~isfield(mogs(i),'date')
        msg{end+1} = ['mogs(' num2str(i) ').date'];
    end
end

%boreholes
nb_boreholes = length(boreholes);
for i = 1:nb_boreholes
    if ~isfield(boreholes(i),'name')
        msg{end+1} = ['boreholes(' num2str(i) ').name'];
    end
end
for i = 1:nb_boreholes
    if ~isfield(boreholes(i),'X')
        msg{end+1} = ['boreholes(' num2str(i) ').X'];
    end
end
for i = 1:nb_boreholes
    if ~isfield(boreholes(i),'Y')
        msg{end+1} = ['boreholes(' num2str(i) ').Y'];
    end
end
for i = 1:nb_boreholes
    if ~isfield(boreholes(i),'Z')
        msg{end+1} = ['boreholes(' num2str(i) ').Z'];
    end
end
for i = 1:nb_boreholes
    if ~isfield(boreholes(i),'Xmax')
        msg{end+1} = ['boreholes(' num2str(i) ').Xmax'];
    end
end
for i = 1:nb_boreholes
    if ~isfield(boreholes(i),'Ymax')
        msg{end+1} = ['boreholes(' num2str(i) ').Ymax'];
    end
end
for i = 1:nb_boreholes
    if ~isfield(boreholes(i),'Zmax')
        msg{end+1} = ['boreholes(' num2str(i) ').Zmax'];
    end
end
for i = 1:nb_boreholes
    if ~isfield(boreholes(i),'Z_surf')
        msg{end+1} = ['boreholes(' num2str(i) ').Z_surf'];
    end
end
for i = 1:nb_boreholes
    if ~isfield(boreholes(i),'scont')
        msg{end+1} = ['boreholes(' num2str(i) ').scont'];
    end
end
for i = 1:nb_boreholes
    if ~isfield(boreholes(i),'acont')
        msg{end+1} = ['boreholes(' num2str(i) ').acont'];
    end
end
for i = 1:nb_boreholes
    if ~isfield(boreholes(i),'fdata')
        msg{end+1} = ['boreholes(' num2str(i) ').fdata'];
    end
end
for i = 1:nb_boreholes
    if ~isfield(boreholes(i),'Z_water')
        msg{end+1} = ['boreholes(' num2str(i) ').Z_water'];
    end
end
for i = 1:nb_boreholes
    if ~isfield(boreholes(i),'diam')
        msg{end+1} = ['boreholes(' num2str(i) ').diam'];
    end
end
if ~isempty(msg)
    msgbox([{'The following fields are not present:           '} {''} {''}  msg],'Warning','warn','modal');
end
