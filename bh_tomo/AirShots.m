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
        
    end
    
    methods
        function obj = AirShots(n)
            obj.name = n;
            obj.data = MogData.empty;
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

