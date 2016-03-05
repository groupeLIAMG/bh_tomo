classdef Mog < handle
    %MOG Class to manage MOGs
    
    properties
        name
        date
        data        % raw data
        av
        ap
        Tx
        Rx
        tt          % traveltime data
        et          % traveltime standard deviation
        tt_done     %
        ttTx
        ttTx_done
        amp_tmin
        amp_tmax
        amp_done
        App
        fcentroid
        scentroid
        tauApp
        tauApp_et
        tauFce
        tauFce_et
        tauHyb
        tauHyb_et
        tau_params
        fw
        f_et
        amp_name_Ldc
        type
        Tx_z_orig
        Rx_z_orig
        fac_dt
        user_fac_dt
        in
        no_traces
        sorted
        TxCosDir
        RxCosDir
        pruneParams
    end
    properties (SetAccess=private)
        ID
    end
    
    methods(Static,Access=private)
        function id = getID(varargin)
            persistent counter;
            if nargin==1
                if isempty(counter)
                    counter = varargin{1};
                else
                    if counter<varargin{1}
                        counter = varargin{1};
                    end
                end
            end
            if isempty(counter)
                counter = 1;
            else
                counter = counter + 1;
            end
            id = counter;
        end
    end
    methods
        function obj = Mog(n)
            obj.name = n;
            obj.date = '';
            obj.data = MogData.empty;
            obj.av = [];        % no du tir aerien pre acquisition
            obj.ap = [];        % no du tir aerien post acquisition
            obj.Tx = 1;         % no du Tx
            obj.Rx = 1;         % no du Rx
            obj.tau_params = [];
            obj.fw = [];                            % donnees filtrees par transf. ondelettes
            obj.f_et = 1;
            obj.amp_name_Ldc = {};
            obj.type = 1;                           % X-hole (1) ou VRP (2)
            obj.fac_dt = 1;
            obj.user_fac_dt = 0;
            obj.pruneParams.sautTx = 0;
            obj.pruneParams.sautRx = 0;
            obj.pruneParams.arrondi = 0;
            obj.pruneParams.use_SB = 0;
            obj.pruneParams.seuil_SB = 0;
            
            obj.ID = Mog.getID();
        end
        function set.name(obj, n)
            if ischar(n)
                obj.name = n;
            else
                error('MOG name must be a string')
            end
        end
        function set.data(obj, d)
            if isa(d, 'MogData')
                obj.data = d;
            else
                error('Data must be a MogData object')
            end
            obj.initialize();
        end
    end
    
    methods (Access=private)
        function initialize(obj)
            if isempty(obj.data)
                return
            end
            obj.date = obj.data.date;
            
            obj.tt = -1*ones(1,obj.data.ntrace);        % temps d'arrivee
            obj.et = -1*ones(1,obj.data.ntrace);        % ecart-type du temps d'arrivee
            obj.tt_done = false(1,obj.data.ntrace);     % temps d'arrivee determine (booleen)
            if isempty(obj.data.tdata)
                obj.ttTx = [];
                obj.ttTx_done = [];
            else
                obj.ttTx = zeros(1,obj.data.ntrace);        % temps d'arrivee
                obj.ttTx_done = false(1,obj.data.ntrace);   % temps d'arrivee determine (booleen)
            end
            obj.amp_tmin = -1*ones(1,obj.data.ntrace);  % t min fenetre de determination Amplitude
            obj.amp_tmax = -1*ones(1,obj.data.ntrace);  % t max fenetre de determination Amplitude
            obj.amp_done = false(1,obj.data.ntrace);    % fenetre de determination Amplitude determinee (booleen)
            obj.App = zeros(1,obj.data.ntrace);         % Amp pic a pic
            obj.fcentroid = zeros(1,obj.data.ntrace);   % freq centroide
            obj.scentroid = zeros(1,obj.data.ntrace);   % variance centroide
            obj.tauApp = -1*ones(1,obj.data.ntrace);       % amplitudes corrigees - amplitude ratio
            obj.tauApp_et = -1*ones(1,obj.data.ntrace);    % ecart-type des amplitudes corrigees
            obj.tauFce = -1*ones(1,obj.data.ntrace);       % amplitudes corrigees - freq centroide
            obj.tauFce_et = -1*ones(1,obj.data.ntrace);    % ecart-type des amplitudes corrigees
            obj.tauHyb = -1*ones(1,obj.data.ntrace);       % amplitudes corrigees - meth. hybride
            obj.tauHyb_et = -1*ones(1,obj.data.ntrace);    % ecart-type des amplitudes corrigees
            obj.Tx_z_orig = obj.data.Tx_z;
            obj.Rx_z_orig = obj.data.Rx_z;
            obj.in = true(1,obj.data.ntrace);
            obj.pruneParams.zmin = min([obj.data.Tx_z obj.data.Rx_z]);
            obj.pruneParams.zmax = max([obj.data.Tx_z obj.data.Rx_z]);
        end
    end
    
    methods (Static)
        function obj = loadobj(a)
            obj = a;
            Mog.getID(obj.ID);  % we must update counter
        end
    end
    
end

