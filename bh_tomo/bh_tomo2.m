function bh_tomo2

file = '';
rep = '';


fs = 11;
vScale = 1;
if ispc
    fs = 9;
    vScale = 0.81;
end

width = 400;
height = 625*vScale;
% get screen size
su = get(groot,'Units');
set(groot,'Units','points')
scnsize = get(groot,'ScreenSize');
pos = [scnsize(3)/2-width/2 scnsize(4)/2-height/2 width height];
set(groot,'Units',su)       % Restore default root screen units


f = figure('Visible','off',...
    'Units','points',...
    'Position',pos,...
    'Tag','fig_bh_tomo2',...
    'Name','bh_tomo',...
    'NumberTitle','off',...
    'Resize','off',...
    'ToolBar','none',...,
    'MenuBar','None',...
    'CloseRequestFcn',@closeWindow);


uicontrol('Style','pushbutton',...
    'String','Magnetic Nano Fluid Saturation',...
    'Units','points',...
    'FontSize',fs+1,...
    'HorizontalAlignment','center',...
    'Position',[40 40*vScale 320 35*vScale],...
    'Enable','off',...
    'Callback',@nanofluid,...
    'Parent',f);

uicontrol('Style','pushbutton',...
    'String','Time-lapse Visualization',...
    'Units','points',...
    'FontSize',fs+1,...
    'HorizontalAlignment','center',...
    'Position',[40 2*40*vScale 320 35*vScale],...
    'Enable','off',...
    'Callback',@timelapse,...
    'Parent',f);

uicontrol('Style','pushbutton',...
    'String','Time-lapse Inversion',...
    'Units','points',...
    'FontSize',fs+1,...
    'HorizontalAlignment','center',...
    'Position',[40 3*40*vScale 320 35*vScale],...
    'Enable','off',...
    'Callback',@tltomo,...
    'Parent',f);

uicontrol('Style','pushbutton',...
    'String','Interpretation (GPR)',...
    'Units','points',...
    'FontSize',fs+1,...
    'HorizontalAlignment','center',...
    'Position',[40 4*40*vScale 320 35*vScale],...
    'Enable','on',...
    'Callback',@interp,...
    'Parent',f);

uicontrol('Style','pushbutton',...
    'String','Inversion',...
    'Units','points',...
    'FontSize',fs+1,...
    'HorizontalAlignment','center',...
    'Position',[40 5*40*vScale 320 35*vScale],...
    'Enable','on',...
    'Callback',@inv,...
    'Parent',f);

uicontrol('Style','pushbutton',...
    'String','Covariance Model',...
    'Units','points',...
    'FontSize',fs+1,...
    'HorizontalAlignment','center',...
    'Position',[40 6*40*vScale 320 35*vScale],...
    'Enable','on',...
    'Callback',@fitCovar,...
    'Parent',f);

uicontrol('Style','pushbutton',...
    'String','Manual Amplitude Picking',...
    'Units','points',...
    'FontSize',fs+1,...
    'HorizontalAlignment','center',...
    'Position',[40 7*40*vScale 320 35*vScale],...
    'Enable','off',...
    'Callback',@amp,...
    'Parent',f);

uicontrol('Style','pushbutton',...
    'String','Manual Traveltime Picking',...
    'Units','points',...
    'FontSize',fs+1,...
    'HorizontalAlignment','center',...
    'Position',[40 8*40*vScale 320 35*vScale],...
    'Enable','on',...
    'Callback',@tt,...
    'Parent',f);

uicontrol('Style','pushbutton',...
    'String','Automatic Traveltime Picking (AIC-CWT)',...
    'Units','points',...
    'FontSize',fs+1,...
    'HorizontalAlignment','center',...
    'Position',[40 9*40*vScale 320 35*vScale],...
    'Enable','off',...
    'Callback',@phase_pick,...
    'Parent',f);

uicontrol('Style','pushbutton',...
    'String','Semi-Automatic Traveltime Picking (x-corr)',...
    'Units','points',...
    'FontSize',fs+1,...
    'HorizontalAlignment','center',...
    'Position',[40 10*40*vScale 320 35*vScale],...
    'Enable','on',...
    'Callback',@pick,...
    'Parent',f);

uicontrol('Style','pushbutton',...
    'String','Database',...
    'Units','points',...
    'FontSize',fs+1,...
    'HorizontalAlignment','center',...
    'Position',[40 11*40*vScale 320 35*vScale],...
    'Callback',@db,...
    'Parent',f);

hdbname = uicontrol('Style','text',...
    'String','',...
    'Units','points',...
    'FontSize',fs+1,...
    'HorizontalAlignment','center',...
    'Position',[40 12.5*40*vScale 320 35*vScale],...
    'BackgroundColor',[1 1 1],...
    'Parent',f);

uicontrol('Style','text',...
    'String',{'BH TOMO','Borehole Radar/Seismic Data Processing Center'},...
    'Units','points',...
    'FontSize',fs+2,...
    'FontWeight','bold',...
    'HorizontalAlignment','center',...
    'Position',[40 14*40*vScale 320 35*vScale],...
    'ForegroundColor',[1 0 0],...
    'Parent',f);

hmenu = uimenu(f,'Label','File');
uimenu(hmenu,'Label','Choose Database ...',...
    'Accelerator','O',...
    'Callback',@chooseDB);
uimenu(hmenu,'Label','Close',...
    'Separator','on',...
    'Accelerator','W',...
    'Callback',@closeWindow);

hmenu2 = uimenu(f,'Label','Edit');
uimenu(hmenu2,'Label','Convert Database ...',...
    'Accelerator','C',...
    'Callback',@convertDB);


f.Visible = 'on';

    function closeWindow(varargin)
        delete(f)
    end

    function chooseDB(varargin)
        [file, rep] = uigetfile('*.mat','Open Database');
        if isequal(file,0)
            return
        end
        hdbname.String = char(file);
    end
    function db(varargin)
        bh_tomo2_db(rep,file,fs);
    end
    function pick(varargin)
        bh_tomo_pick( 'UserData', [rep,file] )
    end
    function phase_pick(varargin)
    end
    function tt(varargin)
        bh_tomo_tt( 'UserData', [rep,file] )
    end
    function amp(varargin)
    end
    function fitCovar(varargin)
        bh_tomo2_covar(rep,file,fs);
    end
    function inv(varargin)
        bh_tomo2_inv(rep,file,fs);
    end
    function interp(varargin)
        bh_tomo2_interp(rep,file,fs);
    end
    function tltomo(varargin)
    end
    function timelapse(varargin)
    end
    function nanofluid(varargin)
    end

    function convertDB(varargin)
        if isempty(file)
            warndlg('Choose a Database First')
            return
        end
        
        fields_v1 = {'mogs','names_mog','air','boreholes','panels',...
            'par_decime','auto_pick','model3d'};
        tmp = load([rep,file]);
        fields_db = fieldnames(tmp);
        if numel(fields_v1)~=numel(fields_db)
            errordlg('Database not in right format')
            return
        end
        test = 1;
        for n=1:numel(fields_v1)
            test = test&& strcmp(fields_v1{n},fields_db{n});
        end
        if test==0
            errordlg('Database not in right format')
            return
        end
        
        % air shots
        air = AirShots.empty(0,numel(tmp.air));
        try
            for n=1:numel(tmp.air)
                air(n) = AirShots(tmp.air(n).name);
                air(n).data = MogData();
                air(n).tt = tmp.air(n).tt;
                air(n).et = tmp.air(n).et;
                air(n).tt_done = tmp.air(n).tt_done;
                air(n).d_TxRx = tmp.air(n).d_TxRx;
                air(n).fac_dt = tmp.air(n).fac_dt;
                air(n).in = tmp.air(n).in;
                air(n).method = tmp.air(n).method;
                air(n).data.ntrace = tmp.air(n).data.ntrace;
                air(n).data.nptsptrc = tmp.air(n).data.nptsptrc;
                air(n).data.rstepsz = tmp.air(n).data.rstepsz;
                air(n).data.cunits = tmp.air(n).data.cunits;
                air(n).data.rnomfreq = tmp.air(n).data.rnomfreq;
                air(n).data.csurvmod = tmp.air(n).data.csurvmod;
                air(n).data.timec = tmp.air(n).data.timec;
                air(n).data.rdata = tmp.air(n).data.rdata;
                air(n).data.timestp = tmp.air(n).data.timestp;
                air(n).data.Tx_x = tmp.air(n).data.Tx_x;
                air(n).data.Tx_y = [];
                air(n).data.Tx_z = tmp.air(n).data.Tx_z;
                air(n).data.Rx_x = tmp.air(n).data.Rx_x;
                air(n).data.Rx_y = [];
                air(n).data.Rx_z = tmp.air(n).data.Rx_z;
                air(n).data.antennas = tmp.air(n).data.antennas;
                air(n).data.synthetique = tmp.air(n).data.synthetique;
                air(n).data.tunits = tmp.air(n).data.tunits;
                air(n).data.TxOffset = tmp.air(n).data.TxOffset;
                air(n).data.RxOffset = tmp.air(n).data.RxOffset;
                air(n).data.comment = tmp.air(n).data.comment;
                air(n).data.date = tmp.air(n).data.cdate;
            end
        catch ME
            errordlg({'Database not in right format:',ME.message})
            return
        end
        
        %boreholes
        boreholes = Borehole.empty(0,numel(tmp.boreholes));
        try
            for n=1:numel(tmp.boreholes)
                boreholes(n) = Borehole(tmp.boreholes(n).name);
                boreholes(n).X = tmp.boreholes(n).X;
                boreholes(n).Y = tmp.boreholes(n).Y;
                boreholes(n).Z = tmp.boreholes(n).Z;
                boreholes(n).Xmax = tmp.boreholes(n).Xmax;
                boreholes(n).Ymax = tmp.boreholes(n).Ymax;
                boreholes(n).Zmax = tmp.boreholes(n).Zmax;
                boreholes(n).Z_surf = tmp.boreholes(n).Z_surf;
                boreholes(n).scont = tmp.boreholes(n).scont;
                boreholes(n).acont = tmp.boreholes(n).acont;
                boreholes(n).diam = tmp.boreholes(n).diam;
                boreholes(n).Z_water = tmp.boreholes(n).Z_water;
                boreholes(n).fdata = tmp.boreholes(n).fdata;
            end
        catch ME
            errordlg({'Database not in right format:',ME.message})
            return
        end
        
        
        % MOGs
        mogs = MOG.empty(0,numel(tmp.mogs));
        try
            for n=1:numel(tmp.mogs)
                mogs(n) = Mog(tmp.mogs(n).name);
                mogs(n).date = tmp.mogs(n).date;
                mogs(n).data = MogData();
                mogs(n).av = tmp.mogs(n).av;
                mogs(n).ap = tmp.mogs(n).ap;
                mogs(n).Tx = tmp.mogs(n).Tx;
                mogs(n).Rx = tmp.mogs(n).Rx;
                mogs(n).tt = tmp.mogs(n).tt;
                mogs(n).et = tmp.mogs(n).et;
                mogs(n).tt_done = tmp.mogs(n).tt_done;
                mogs(n).ttTx = tmp.mogs(n).ttTx;
                mogs(n).ttTx_done = tmp.mogs(n).ttTx_done;
                mogs(n).amp_tmin = tmp.mogs(n).amp_tmin;
                mogs(n).amp_tmax = tmp.mogs(n).amp_tmax;
                mogs(n).amp_done = tmp.mogs(n).amp_done;
                mogs(n).App = tmp.mogs(n).App;
                mogs(n).fcentroid = tmp.mogs(n).fcentroid;
                mogs(n).scentroid = tmp.mogs(n).scentroid;
                mogs(n).tauApp = tmp.mogs(n).tauApp;
                mogs(n).tauApp_et = tmp.mogs(n).tauApp_et;
                mogs(n).tauFce = tmp.mogs(n).tauFce;
                mogs(n).tauFce_et = tmp.mogs(n).tauFce_et;
                mogs(n).tauHyb = tmp.mogs(n).tauHyb;
                mogs(n).tauHyb_et = tmp.mogs(n).tauHyb_et;
                mogs(n).tau_params = tmp.mogs(n).tau_params;
                mogs(n).fw = tmp.mogs(n).fw;
                mogs(n).f_et = tmp.mogs(n).f_et;
                mogs(n).amp_name_Ldc = tmp.mogs(n).amp_name_Ldc;
                mogs(n).type = tmp.mogs(n).type;
                mogs(n).Tx_z_orig = tmp.mogs(n).Tx_z_orig;
                mogs(n).Rx_z_orig = tmp.mogs(n).Rx_z_orig;
                mogs(n).fac_dt = tmp.mogs(n).fac_dt;
                mogs(n).user_fac_dt = tmp.mogs(n).user_fac_dt;
                mogs(n).in = tmp.mogs(n).in;
                mogs(n).no_traces = tmp.mogs(n).no_traces;
                mogs(n).sorted = tmp.mogs(n).sorted;
                mogs(n).TxCosDir = tmp.mogs(n).TxCosDir;
                mogs(n).RxCosDir = tmp.mogs(n).RxCosDir;
                % pruneParams (par_decime) not converted
                
                mogs(n).data.ntrace = tmp.mogs(n).data.ntrace;
                mogs(n).data.nptsptrc = tmp.mogs(n).data.nptsptrc;
                mogs(n).data.rstepsz = tmp.mogs(n).data.rstepsz;
                mogs(n).data.cunits = tmp.mogs(n).data.cunits;
                mogs(n).data.rnomfreq = tmp.mogs(n).data.rnomfreq;
                mogs(n).data.csurvmod = tmp.mogs(n).data.csurvmod;
                mogs(n).data.timec = tmp.mogs(n).data.timec;
                mogs(n).data.rdata = tmp.mogs(n).data.rdata;
                mogs(n).data.timestp = tmp.mogs(n).data.timestp;
                mogs(n).data.Tx_x = tmp.mogs(n).data.Tx_x;
                if isfield(tmp.mogs(n).data,'Tx_y')
                    mogs(n).data.Tx_y = tmp.mogs(n).data.Tx_y;
                else
                    mogs(n).data.Tx_y = zeros(size(tmp.mogs(n).data.Tx_x));
                end
                mogs(n).data.Tx_z = tmp.mogs(n).data.Tx_z;
                mogs(n).data.Rx_x = tmp.mogs(n).data.Rx_x;
                if isfield(tmp.mogs(n).data,'Rx_y')
                    mogs(n).data.Rx_y = tmp.mogs(n).data.Rx_y;
                else
                    mogs(n).data.Rx_y = zeros(size(tmp.mogs(n).data.Rx_x));
                end
                mogs(n).data.Rx_z = tmp.mogs(n).data.Rx_z;
                mogs(n).data.antennas = tmp.mogs(n).data.antennas;
                mogs(n).data.synthetique = tmp.mogs(n).data.synthetique;
                mogs(n).data.tunits = tmp.mogs(n).data.tunits;
                mogs(n).data.TxOffset = tmp.mogs(n).data.TxOffset;
                mogs(n).data.RxOffset = tmp.mogs(n).data.RxOffset;
                mogs(n).data.comment = tmp.mogs(n).data.comment;
            end
        catch ME
            errordlg({'Database not in right format:',ME.message})
            return
        end
        
        %
        % Models
        models = Model.empty; %#ok<NASGU>
        % TODO

        auto_pick = tmp.auto_pick; %#ok<NASGU>
        
        [file2, rep2] = uiputfile('*.mat','Save New Database');
        if isequal(file2,0)
            return
        end
        file=file2;
        rep=rep2;
        save([rep,file],'names_mog','mogs','air','boreholes','models','auto_pick')
        hdbname.String = char(file);
    end

end

