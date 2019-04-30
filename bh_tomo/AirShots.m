classdef AirShots < handle
    %UNTITLED Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        name
        data
        tt
        et
        tt_done
        d_TxRx
        fac_dt
        in
        method
        traces
    end
    
    methods
        function obj = AirShots(n)
            if isstring(n) || ischar(n)
                obj.name = n;
                obj.data = MogData.empty;
            elseif isstruct(n)
                if isfield(n, 'name')
                    obj.name = n.name;
                elseif isfield(n, 'nom')
                    obj.name = n.nom;
                else
                    error('Invalid input')
                end
                obj.data = MogData(n.data);
                obj.traces = n.data.rdata;
                obj.tt = n.tt;
                obj.et = n.et;
                obj.tt_done = n.tt_done;
                obj.d_TxRx = n.d_TxRx;
                obj.fac_dt = n.fac_dt;
                obj.in = n.in;
                if isfield(n, 'method')
                    obj.method = n.method;
                else
                    obj.method = 'fixed_antenna';
                end
                obj.traces = [];
            elseif isa(n, 'AirShots')
                obj.name = n.name;
                obj.data = MogData(n.data);
                obj.traces = n.traces;
                obj.tt = n.tt;
                obj.et = n.et;
                obj.tt_done = n.tt_done;
                obj.d_TxRx = n.d_TxRx;
                obj.fac_dt = n.fac_dt;
                obj.in = n.in;
                obj.method = n.method;
            else
                error('Invalid input')
            end
        end
        function set.name(obj, n)
            if ischar(n)
                obj.name = n;
            else
                error('name must be a string')
            end
        end
        function set.data(obj, d)
            if isa(d, 'MogData')
                obj.data = d;
            else
                error('data must be a MogData object')
            end
            obj.initialize();
        end
    end
    
    methods (Access=private)
        function initialize(obj)
            if isempty(obj.data)
                return
            end
            
            obj.tt = -1*ones(1,obj.data.ntrace);        % temps d'arrivee
            obj.et = -1*ones(1,obj.data.ntrace);        % ecart-type du temps d'arrivee
            obj.tt_done = false(1,obj.data.ntrace);     % temps d'arrivee determine (booleen)
        end
    end
    
end

