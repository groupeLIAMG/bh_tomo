function [boreholes, air, mogs, models] = load_db_h5(filename)
%
%

value = H5F.is_hdf5(filename);
if value <= 0
    fprintf([filename,' is not an HDF5 file\n']);
end

% Boreholes

info = h5info(filename,'/boreholes');
boreholes = Borehole.empty(0, numel(info.Groups));
for n=1:numel(info.Groups)
    name = char(h5readatt(filename, info.Groups(n).Name, 'name'));
    boreholes(n) = Borehole(name);
    boreholes(n).X = h5readatt(filename, info.Groups(n).Name, 'X');
    boreholes(n).Y = h5readatt(filename, info.Groups(n).Name, 'Y');
    boreholes(n).Z = h5readatt(filename, info.Groups(n).Name, 'Z');
    boreholes(n).Xmax = h5readatt(filename, info.Groups(n).Name, 'Xmax');
    boreholes(n).Ymax = h5readatt(filename, info.Groups(n).Name, 'Ymax');
    boreholes(n).Zmax = h5readatt(filename, info.Groups(n).Name, 'Zmax');
    boreholes(n).Z_surf = h5readatt(filename, info.Groups(n).Name, 'Z_surf');
    boreholes(n).Z_water = h5readatt(filename, info.Groups(n).Name, 'Z_water');
    boreholes(n).fdata = h5read(filename, [info.Groups(n).Name, '/fdata'])';
    boreholes(n).scont = h5read(filename, [info.Groups(n).Name, '/scont'])';
    boreholes(n).acont = h5read(filename, [info.Groups(n).Name, '/acont'])';
end

% Air shots

info = h5info(filename,'/air_shots');
air = AirShots.empty(0, numel(info.Groups));
for n=1:numel(info.Groups)
    name = char(h5readatt(filename, info.Groups(n).Name, 'name'));
    air(n) = AirShots(name);
    air(n).fac_dt = h5readatt(filename, info.Groups(n).Name, 'fac_dt');
    air(n).method = char(h5readatt(filename, info.Groups(n).Name, 'method'));
    air(n).d_TxRx = h5read(filename, [info.Groups(n).Name, '/d_TxRx'])';
    air(n).tt = h5read(filename, [info.Groups(n).Name, '/tt'])';
    air(n).et = h5read(filename, [info.Groups(n).Name, '/et'])';
    air(n).tt_done = h5read(filename, [info.Groups(n).Name, '/tt_done'])';
    
    tmp = h5read(filename, [info.Groups(n).Name, '/in_vect'])';
    air(n).in = false(1, numel(tmp));
    for nn=1:numel(tmp)
        if strcmp(tmp{nn}, 'TRUE')
            air(n).in(nn) = 1;
        end
    end

    % mog data
    air(n).data = MogData();
    air(n).data.RxOffset = h5readatt(filename, [info.Groups(n).Name, '/data'], 'RxOffset');
    air(n).data.TxOffset = h5readatt(filename, [info.Groups(n).Name, '/data'], 'TxOffset');
    air(n).data.antennas = char(h5readatt(filename, [info.Groups(n).Name, '/data'], 'antennas'));
    air(n).data.comment = char(h5readatt(filename, [info.Groups(n).Name, '/data'], 'comment'));
    air(n).data.csurvmod = char(h5readatt(filename, [info.Groups(n).Name, '/data'], 'csurvmod'));
    air(n).data.cunits = char(h5readatt(filename, [info.Groups(n).Name, '/data'], 'cunits'));
    air(n).data.date = char(h5readatt(filename, [info.Groups(n).Name, '/data'], 'date'));
    air(n).data.nptsptrc = h5readatt(filename, [info.Groups(n).Name, '/data'], 'nptsptrc');
    air(n).data.ntrace = h5readatt(filename, [info.Groups(n).Name, '/data'], 'ntrace');
    air(n).data.rnomfreq = h5readatt(filename, [info.Groups(n).Name, '/data'], 'rnomfreq');
    air(n).data.rstepsz = h5readatt(filename, [info.Groups(n).Name, '/data'], 'rstepsz');
    air(n).data.synthetique = h5readatt(filename, [info.Groups(n).Name, '/data'], 'synthetique');
    air(n).data.timec = h5readatt(filename, [info.Groups(n).Name, '/data'], 'timec');
    air(n).data.tunits = char(h5readatt(filename, [info.Groups(n).Name, '/data'], 'tunits'));
    air(n).data.rdata = h5read(filename, [info.Groups(n).Name, '/data/rdata'])';
    air(n).data.timestp = h5read(filename, [info.Groups(n).Name, '/data/timestp'])';
    air(n).data.Rx_x = h5read(filename, [info.Groups(n).Name, '/data/Rx_x'])';
    air(n).data.Rx_y = h5read(filename, [info.Groups(n).Name, '/data/Rx_y'])';
    air(n).data.Rx_z = h5read(filename, [info.Groups(n).Name, '/data/Rx_z'])';
    air(n).data.Tx_x = h5read(filename, [info.Groups(n).Name, '/data/Tx_x'])';
    air(n).data.Tx_y = h5read(filename, [info.Groups(n).Name, '/data/Tx_y'])';
    air(n).data.Tx_z = h5read(filename, [info.Groups(n).Name, '/data/Tx_z'])';
end

% Mogs

info = h5info(filename,'/mogs');
mogs = Mog.empty(0,numel(info.Groups));
for n=1:numel(info.Groups)
    name = char(h5readatt(filename, info.Groups(n).Name, 'name'));
    
    mogs(n) = Mog(name);
    mogs(n).amp_name_Ldc = char(h5readatt(filename, info.Groups(n).Name, 'amp_name_Ldc'));
    mogs(n).date = char(h5readatt(filename, info.Groups(n).Name, 'date'));
    mogs(n).f_et = h5readatt(filename, info.Groups(n).Name, 'f_et');
    mogs(n).fac_dt = h5readatt(filename, info.Groups(n).Name, 'fac_dt');
    mogs(n).type  = h5readatt(filename, info.Groups(n).Name, 'type')+1;
    mogs(n).useAirShots = h5readatt(filename, info.Groups(n).Name, 'useAirShots');
    mogs(n).user_fac_dt = h5readatt(filename, info.Groups(n).Name, 'user_fac_dt');
    mogs(n).tt = h5read(filename, [info.Groups(n).Name, '/tt'])';
    mogs(n).et = h5read(filename, [info.Groups(n).Name, '/et'])';
    mogs(n).App = h5read(filename, [info.Groups(n).Name, '/App'])';
    mogs(n).RxCosDir = h5read(filename, [info.Groups(n).Name, '/RxCosDir'])';
    mogs(n).TxCosDir = h5read(filename, [info.Groups(n).Name, '/TxCosDir'])';
    mogs(n).Rx_z_orig = h5read(filename, [info.Groups(n).Name, '/Rx_z_orig'])';
    mogs(n).Tx_z_orig = h5read(filename, [info.Groups(n).Name, '/Tx_z_orig'])';
    mogs(n).amp_tmax = h5read(filename, [info.Groups(n).Name, '/amp_tmax'])';
    mogs(n).amp_tmin = h5read(filename, [info.Groups(n).Name, '/amp_tmin'])';
    mogs(n).fcentroid = h5read(filename, [info.Groups(n).Name, '/fcentroid'])';
    mogs(n).scentroid = h5read(filename, [info.Groups(n).Name, '/scentroid'])';
    mogs(n).fw = h5read(filename, [info.Groups(n).Name, '/fw'])';
    mogs(n).tauApp = h5read(filename, [info.Groups(n).Name, '/tauApp'])';
    mogs(n).tauApp_et = h5read(filename, [info.Groups(n).Name, '/tauApp_et'])';
    mogs(n).tauFce = h5read(filename, [info.Groups(n).Name, '/tauFce'])';
    mogs(n).tauFce_et = h5read(filename, [info.Groups(n).Name, '/tauFce_et'])';
    mogs(n).tauHyb = h5read(filename, [info.Groups(n).Name, '/tauHyb'])';
    mogs(n).tauHyb_et = h5read(filename, [info.Groups(n).Name, '/tauHyb_et'])';
    mogs(n).tau_params = h5read(filename, [info.Groups(n).Name, '/tau_params'])';
    mogs(n).ttTx = h5read(filename, [info.Groups(n).Name, '/ttTx'])';
    
    tmp = h5read(filename, [info.Groups(n).Name, '/tt_done'])';
    mogs(n).tt_done =  false(1, numel(tmp));
    for nn=1:numel(tmp)
        if strcmp(tmp{nn}, 'TRUE')
            mogs(n).tt_done(nn) = 1;
        end
    end
    tmp = h5read(filename, [info.Groups(n).Name, '/in_vect'])';
    mogs(n).in =  false(1, numel(tmp));
    for nn=1:numel(tmp)
        if strcmp(tmp{nn}, 'TRUE')
            mogs(n).in(nn) = 1;
        end
    end
    tmp = h5read(filename, [info.Groups(n).Name, '/ttTx_done'])';
    mogs(n).ttTx_done =  false(1, numel(tmp));
    for nn=1:numel(tmp)
        if strcmp(tmp{nn}, 'TRUE')
            mogs(n).ttTx_done(nn) = 1;
        end
    end
    tmp = h5read(filename, [info.Groups(n).Name, '/amp_done'])';
    mogs(n).amp_done =  false(1, numel(tmp));
    for nn=1:numel(tmp)
        if strcmp(tmp{nn}, 'TRUE')
            mogs(n).amp_done(nn) = 1;
        end
    end

    % mog data
    mogs(n).data = MogData();
    mogs(n).data.RxOffset = h5readatt(filename, [info.Groups(n).Name, '/data'], 'RxOffset');
    mogs(n).data.TxOffset = h5readatt(filename, [info.Groups(n).Name, '/data'], 'TxOffset');
    mogs(n).data.antennas = char(h5readatt(filename, [info.Groups(n).Name, '/data'], 'antennas'));
    mogs(n).data.comment = char(h5readatt(filename, [info.Groups(n).Name, '/data'], 'comment'));
    mogs(n).data.csurvmod = char(h5readatt(filename, [info.Groups(n).Name, '/data'], 'csurvmod'));
    mogs(n).data.cunits = char(h5readatt(filename, [info.Groups(n).Name, '/data'], 'cunits'));
    mogs(n).data.date = char(h5readatt(filename, [info.Groups(n).Name, '/data'], 'date'));
    mogs(n).data.nptsptrc = h5readatt(filename, [info.Groups(n).Name, '/data'], 'nptsptrc');
    mogs(n).data.ntrace = h5readatt(filename, [info.Groups(n).Name, '/data'], 'ntrace');
    mogs(n).data.rnomfreq = h5readatt(filename, [info.Groups(n).Name, '/data'], 'rnomfreq');
    mogs(n).data.rstepsz = h5readatt(filename, [info.Groups(n).Name, '/data'], 'rstepsz');
    mogs(n).data.synthetique = h5readatt(filename, [info.Groups(n).Name, '/data'], 'synthetique');
    mogs(n).data.timec = h5readatt(filename, [info.Groups(n).Name, '/data'], 'timec');
    mogs(n).data.tunits = char(h5readatt(filename, [info.Groups(n).Name, '/data'], 'tunits'));
    mogs(n).data.rdata = h5read(filename, [info.Groups(n).Name, '/data/rdata'])';
    mogs(n).data.timestp = h5read(filename, [info.Groups(n).Name, '/data/timestp'])';
    mogs(n).data.Rx_x = h5read(filename, [info.Groups(n).Name, '/data/Rx_x'])';
    mogs(n).data.Rx_y = h5read(filename, [info.Groups(n).Name, '/data/Rx_y'])';
    mogs(n).data.Rx_z = h5read(filename, [info.Groups(n).Name, '/data/Rx_z'])';
    mogs(n).data.Tx_x = h5read(filename, [info.Groups(n).Name, '/data/Tx_x'])';
    mogs(n).data.Tx_y = h5read(filename, [info.Groups(n).Name, '/data/Tx_y'])';
    mogs(n).data.Tx_z = h5read(filename, [info.Groups(n).Name, '/data/Tx_z'])';
    
    tmp = char(h5readatt(filename, [info.Groups(n).Name, '/Rx'], 'name'));
    for nn=1:numel(boreholes)
        if strcmp(tmp, boreholes(nn).name)
            mogs(n).Rx = nn;
        end
    end
    tmp = char(h5readatt(filename, [info.Groups(n).Name, '/Tx'], 'name'));
    for nn=1:numel(boreholes)
        if strcmp(tmp, boreholes(nn).name)
            mogs(n).Tx = nn;
        end
    end
    tmp = char(h5readatt(filename, [info.Groups(n).Name, '/av'], 'name'));
    for nn=1:numel(air)
        if strcmp(tmp, air(nn).name)
            mogs(n).av = nn;
        end
    end
    tmp = char(h5readatt(filename, [info.Groups(n).Name, '/ap'], 'name'));
    for nn=1:numel(air)
        if strcmp(tmp, air(nn).name)
            mogs(n).ap = nn;
        end
    end
end

mogs(1).no_traces = 1:mogs(1).data.ntrace;
for n=2:numel(mogs)
    mogs(n).no_traces = mogs(n-1).no_traces(end) + (1:mogs(n).data.ntrace);
end

% TODO: complete models
info = h5info(filename,'/models');
models = Model.empty(0,numel(info.Groups));
for n=1:numel(info.Groups)
    name = char(h5readatt(filename, info.Groups(n).Name, 'name'));
    
    models(n) = Model(name);
    list_mogs = h5info(filename, [info.Groups(n).Name, '/_list_mogs']);
    models(n).mogs = zeros(1, numel(list_mogs.Groups));
    models(n).boreholes = [];
    for nm = 1:numel(list_mogs.Groups)
        name = char(h5readatt(filename, list_mogs.Groups(nm).Name, 'name'));
        for mo = 1:numel(mogs)
            if strcmp(mogs(mo).name, name)
                models(n).mogs(nm) = mo;
                models(n).boreholes = [models(n).boreholes mogs(mo).Tx mogs(mo).Rx];
            end
        end
    end
    models(n).boreholes = unique(models(n).boreholes);
    
end