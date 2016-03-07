function addMog(obj,varargin)

if strcmp(obj.handles.popupTx.String,' ')
    warndlg('Define boreholes first.')
    return
end
old_rep = pwd;
if ~isempty(obj.data_rep)
    cd( obj.data_rep );
end
[file, rep, filterindex] = uigetfile({
    '*.rad;*.RAD', 'RAMAC file *.rad';...
    '*.mat', 'Magnetic NanoMappers data in Matlab format (*.mat)';...
    '*.hd;*.HD', 'EKKO file *.hd';...
    '*.mat', 'RAMAC file in Matlab format *.mat';...
    '*.dat', 'MSIS file *.dat';...
    '*.su;*.SU', 'Seismic Un*x file *.su';...
    '*.sgy;*.segy;*.SGY;*.SEGY', 'SEG Y file'},...
    'Open data file');
cd( old_rep );
if isequal(file,0) || isequal(rep,0)
    return
end
obj.data_rep = rep;

no = length(obj.mogs)+1;
switch filterindex
    case 1
        name = file(1:end-4);
        obj.mogs(no) = Mog(name);
        obj.mogs(no).data = MogData();
        obj.mogs(no).data.readRAMAC([rep,name]);
    case 2
        name = file(1:end-4);
        obj.mogs(no) = Mog(name);
        obj.mogs(no).data = MogData();
        obj.mogs(no).data.readNanoMapper([rep,name]);
    case 3
        name = file(1:end-3);
        obj.mogs(no) = Mog(name);
        obj.mogs(no).data = MogData();
        obj.mogs(no).data.readEKKO([rep,name]);
    case 4
        name = file(1:end-4);
        obj.mogs(no) = Mog(name);
        obj.mogs(no).data = MogData();
        obj.mogs(no).data.readMAT([rep,name]);
    case 5
        name = file(1:end-4);
        obj.mogs(no) = Mog(name);
        obj.mogs(no).data = MogData();
        obj.mogs(no).data.readMSIS([rep,name]);
    case 6
        name = file(1:end-4);
        obj.mogs(no) = Mog(name);
        obj.mogs(no).data = MogData();
        obj.mogs(no).data.readSU([rep,name]);
    case 7
        name = file(1:end-4);
        obj.mogs(no) = Mog(name);
        obj.mogs(no).data = MogData();
        obj.mogs(no).data.readSEGY([rep,name]);
end
obj.mogs(no).initialize();

debut=0;
if no>1
    debut = obj.mogs(no-1).no_traces( length(obj.mogs(no-1).no_traces) );
end
obj.mogs(no).no_traces = debut + (1:obj.mogs(no).data.ntrace);
obj.mogs(no).sorted = false;

obj.updateCoords(no);

obj.updateList()
obj.updateUIfields()
obj.notify('mogAdded')
end
