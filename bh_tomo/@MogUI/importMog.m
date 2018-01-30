function importMog(obj,varargin)

[no_mog,db_file] = chooseMOG();
if no_mog==0
    return
end
tmp = load(db_file,'mogs','boreholes','air');
if ~isa(tmp.mogs,'Mog')
    error('imported mogs are not instances of class Mog')
end
if ~isa(tmp.boreholes,'Borehole')
    error('imported boreholes are not instances of class Borehole')
end
if ~isa(tmp.air,'AirShots')
    error('imported air shot data are not instances of class AirShots')
end
mog = tmp.mogs(no_mog);
for n=1:length(names_mog)
    if strcmp(obj.mogs(n).name, mog.name)
        ButtonName=questdlg('MOG already exists.  Add?');
        switch lower(ButtonName)
            case 'yes'
                mog.name=[mog.name,'_2'];
                warndlg(['Imported MOG name changed to ',mog.name])
            otherwise
                return
        end
    end
end

% check if air shots are present
for n=1:length(obj.air)
    if strcmp(obj.air(n).name, tmp.air(mog.av).name)
        tmp.air(mog.av).name = [tmp.air(mog.av).name,'_2'];
        warndlg(['Air shot of imported MOG name changed to ',tmp.air(mog.av).name])
    end
    if strcmp(obj.air(n).name, tmp.air(mog.ap).name)
        tmp.air(mog.ap).name = [tmp.air(mog.ap).name,'_2'];
        warndlg(['Air shot of imported MOG name changed to ',tmp.air(mog.ap).name])
    end
end
obj.air = [obj.air tmp.air(mog.av)];
mog.av = length(obj.air);
obj.air = [obj.air tmp.air(mog.ap)];
mog.ap = length(obj.air);

% check if boreholes are defined already
addTx = true;
addRx = true;
mogTx = mog.Tx;
mogRx = mog.Rx;
for n=1:length(obj.boreholes)
    if strcmp(obj.boreholes(n).name, tmp.boreholes(mogTx).name)
        addTx = false;
        mog.Tx = n;
    end
    if strcmp(obj.boreholes(n).name, tmp.boreholes(mogRx).name)
        addRx = false;
        mog.Rx = n;
    end
end

if addTx
    obj.boreholes = [obj.boreholes tmp.boreholes(mog.Tx)];
    mog.Tx = length(obj.boreholes);
end
if addRx
    obj.boreholes = [obj.boreholes tmp.boreholes(mog.Rx)];
    mog.Rx = length(obj.boreholes);
end

if isempty(obj.mogs)
    mog.no_traces = 1:mog.data.ntrace;
else
    mog.no_traces = obj.mogs(end).no_traces(end) + (1:mog.data.ntrace);
end

obj.mogs = [obj.mogs mog];

obj.updateList()
obj.updateUIfields()
obj.notify('mogAdded')
end
