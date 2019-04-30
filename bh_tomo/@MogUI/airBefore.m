function airBefore(obj,varargin)

old_rep = pwd;
if ~isempty(obj.data_rep)
    cd( obj.data_rep );
end
[file, rep] = uigetfile('*.rad;*.RAD','RAMAC file *.rad',...
    'Open t0 air shot before survey');
cd( old_rep );
if isequal(file,0) || isequal(rep,0)
    return
end
obj.data_rep = rep;
no = obj.handles.listMogs.Value;
if no>0 && no<=length(obj.mogs)
    
    name = file(1:end-4);
    found = false;
    for n=1:length(obj.air)
        if strcmp( char(obj.air(n).name), char( name ) )
            obj.mogs(no).av = n;
            found = true;
            break
        end
    end
    
    if ~found
        n=length(obj.air)+1;
        
        data = MogData();
        data.readRAMAC([rep,name]);
        
        l = inputdlg('Enter distance between Tx and Rx');
        try
            d = eval(l{1});
        catch
            try
                d = eval(['[',l{1},']']);
                if length(d) ~= data.ntrace
                    ME = MException('MATLAB:INVFMT', 'Number of positions inconsistent with number of traces');
                    throw(ME);
                end
            catch Error
                rethrow(Error);
            end
        end
        
        obj.air(n) = AirShots(char(name));
        obj.air(n).data = data;
        obj.air(n).traces = data.rdata;
        obj.air(n).tt = -1*ones(1,data.ntrace);        % temps d'arrivee
        obj.air(n).et = -1*ones(1,data.ntrace);        % ecart-type du temps d'arrivee
        obj.air(n).tt_done = false(1,data.ntrace);     % temps d'arrivee determine (booleen)
        obj.air(n).d_TxRx = d;
        obj.air(n).fac_dt = 1;
        obj.air(n).in = true(1,data.ntrace);
        if length(d)==1
            obj.air(n).method = 'fixed_antenna';           % 'fixed_antenna' ou 'walkaway'
        else
            obj.air(n).method = 'walkaway';
        end
        obj.handles.textAirBefore.String = obj.air(n).name;
        obj.mogs(no).av = n;
    end
    
end
end

