classdef MogData < handle
    %MOGDATA  Multi-offset data

    properties
        ntrace                % number of traces in MOG
        nptsptrc              % number of sample per trace
        rstepsz               % theoretical spatial step size between traces
        cunits                % spatial units
        rnomfreq              % nominal frequency of source
        csurvmod              % type of survey
        timec                 % sampling period
        rdata                 % raw data
        timestp               % time vector
        Tx_x                  % x coordinate of source
        Tx_y                  % y coordinate of source
        Tx_z                  % z coordinate of source
        Rx_x                  % x coordinate of receiver
        Rx_y                  % y coordinate of receiver
        Rx_z                  % z coordinate of receiver
        antennas              % type of antenna
        synthetique           % true of synthetic data
        tunits                % time units
        TxOffset              % Tx antenna feedpoint offset
        RxOffset              % Rx antenna feedpoint offset
        comment
        date                  % date of survey
        tdata
    end

    methods
        function obj = MogData
            obj.ntrace = 0;
            obj.nptsptrc = 0;
            obj.rstepsz = 0;
            obj.cunits = '';
            obj.rnomfreq = 0;
            obj.csurvmod = '';
            obj.timec = 0;
            obj.rdata = 0;
            obj.timestp = 0;
            obj.Tx_x = 0;
            obj.Tx_y = 0;
            obj.Tx_z = 0;
            obj.Rx_x = 0;
            obj.Rx_y = 0;
            obj.Rx_z = 0;
            obj.antennas = '';
            obj.synthetique = 0;
            obj.tunits = '';
            obj.TxOffset = 0;
            obj.RxOffset = 0;
            obj.comment = '';
            obj.date = '';
            obj.tdata = [];
        end
        
        % Make a copy of a handle object.
        function new = copy(obj)
            % Instantiate new object of the same class.
            new = feval(class(obj));
            
            % Copy all non-hidden properties.
            p = properties(obj);
            for i = 1:length(p)
                new.(p{i}) = obj.(p{i});
            end
        end

        function s = getSubset(obj,ind)
            s = MogData;
            s.ntrace = sum(ind);
            s.nptsptrc = obj.nptsptrc;
            s.rstepsz = obj.rstepsz;
            s.cunits = obj.cunits;
            s.rnomfreq = obj.rnomfreq;
            s.csurvmod = obj.csurvmod;
            s.timec = obj.timec;
            s.rdata = obj.rdata(:,ind);
            s.timestp = obj.timestp;
            s.Tx_x = obj.Tx_x(ind);
            s.Tx_y = obj.Tx_y(ind);
            s.Tx_z = obj.Tx_z(ind);
            s.Rx_x = obj.Rx_x(ind);
            s.Rx_y = obj.Rx_y(ind);
            s.Rx_z = obj.Rx_z(ind);
            s.antennas = obj.antennas;
            s.synthetique = obj.synthetique;
            s.tunits = obj.tunits;
            s.TxOffset = obj.TxOffset;
            s.RxOffset = obj.RxOffset;
            s.comment = obj.comment;
            s.date = obj.date;
            s.tdata = obj.tdata;
        end
        function readRAMAC(obj, basename)
            obj.readRAD( basename )
            obj.readRD3( basename )
            obj.readTLF( basename )

            % Offset between antenna feed point and top of antenna
            obj.TxOffset = 0;
            obj.RxOffset = 0;
            if obj.synthetique == 0
                if obj.rnomfreq == 100
                    obj.TxOffset = 0.665;
                    obj.RxOffset = 0.665;
                elseif obj.rnomfreq == 250
                    obj.TxOffset = 0.325;
                    obj.RxOffset = 0.365;
                end
            end

            obj.Tx_z = obj.Tx_z(1:obj.ntrace);
            obj.Rx_z = obj.Rx_z(1:obj.ntrace);

            obj.Rx_x = zeros(1,obj.ntrace);
            obj.Rx_y = zeros(1,obj.ntrace);
            obj.Tx_x = zeros(1,obj.ntrace);
            obj.Tx_y = zeros(1,obj.ntrace);
            obj.cunits = 'm';
            obj.tunits = 'ns';
        end
        function readNanoMapper(obj, filename)
            load(filename);
            if exist('data','var') ~= 1
                error('The file does not contain a variable named ''data''.')
            end
            if ~isfield(data, 'comment') %#ok<NODEF>
                data.comment = '';
            end
            if strcmp(data.tunits, 's') == 1
                data.timec = data.timec*1e9;
                data.timestp = data.timestp*1e9;
                data.tunits = 'ns';
            end
            obj.ntrace = data.ntrace;
            obj.nptsptrc = data.nptsptrc;
            obj.rstepsz = data.rstepsz;
            obj.cunits = data.cunits;
            obj.rnomfreq = data.rnomfreq;
            obj.csurvmod = data.csurvmod;
            obj.timec = data.timec;
            obj.rdata = data.rdata;
            obj.timestp = data.timestp;
            obj.Tx_x = data.Tx_x;
            obj.Tx_y = zeros(1,obj.ntrace);
            obj.Tx_z = data.Tx_z;
            obj.Rx_x = data.Rx_x;
            obj.Rx_y = zeros(1,obj.ntrace);
            obj.Rx_z = data.Rx_z;
            obj.antennas = data.antennas;
            obj.synthetique = data.synthetique;
            obj.tunits = data.tunits;
            if isfield(data, 'date')
                obj.date = data.date;
            end
            if isfield(data, 'tdata')
                obj.tdata = data.tdata;
            end
        end
        function readEKKO(obj, basename)
            obj.readHD( basename );
            obj.readDT1( basename );
        end
        function readMAT(obj, basename)
            data = load(basename);
            obj.ntrace = data.ntrace;
            obj.nptsptrc = data.nptsptrc;
            obj.rstepsz = data.rstepsz;
            obj.cunits = data.cunits;
            obj.rnomfreq = data.rnomfreq;
            obj.csurvmod = data.csurvmod;
            obj.timec = data.timec;
            obj.rdata = data.rdata;
            obj.timestp = data.timestp;
            obj.Tx_x = data.Tx_x;
            obj.Tx_y = zeros(1,obj.ntrace);
            obj.Tx_z = data.Tx_z;
            obj.Rx_x = data.Rx_x;
            obj.Rx_y = zeros(1,obj.ntrace);
            obj.Rx_z = data.Rx_z;
            obj.antennas = data.antennas;
            obj.synthetique = data.synthetique;
            obj.tunits = data.tunits;
            if isfield(data, 'date')
                obj.date = data.date;
            end
            if isfield(data, 'tdata')
                obj.tdata = data.tdata;
            end
        end
        function readMSIS(obj, basename)

            fid=fopen(basename,'rt');
            if fid < 2
                fichier = [basename,'.dat'];
                fid=fopen(fichier,'rt');
                if fid < 2
                    fichier = [basename,'.DAT'];
                    fid=fopen(fichier,'rt');
                    if fid < 2
                        error(['Cannot open MSIS file: ',basename])
                    end
                end
            end

            obj.antennas = 'Microsis';
            obj.tunits = 'ms';
            obj.cunits = 'm';

            while 1
                tline = fgetl(fid);
                if ~ischar(tline), break, end
                if ~isempty(strfind(tline,'Number of channel     :	 '))
                    [debut,fin] = regexp(tline,'\d+');
                    obj.ntrace = str2double( tline(debut:fin) );
                elseif ~isempty(strfind(tline,'Sampling rate         :	 '))
                    [debut,fin] = regexp(tline,'\d+');
                    obj.timec = str2double( tline(debut:fin) );  % MHz
                elseif ~isempty(strfind(tline,'Number of saved pts  :	 '))
                    [debut,fin] = regexp(tline,'\d+');
                    obj.nptsptrc = str2double( tline(debut:fin) );
                elseif ~isempty(strfind(tline,'***************************'))
                    position=ftell(fid);
                end
            end
            fseek (fid, position, 'bof');
            obj.rdata = fscanf(fid, '%f',[obj.ntrace inf])';
            fclose(fid);

            obj.timec = 1e-3/obj.timec;
            obj.timestp = obj.timec*(1:obj.nptsptrc);

            readTLF( basename );
        end
        function readSU(obj, basename)
            if exist('ReadSu','file')==0
                error('You need to install SegyMAT to read SU files')
            end
            [obj.rdata, tr_header, header] = ReadSu(basename);
            obj.ntrace = length(tr_header);
            obj.nptsptrc = header.ns;
            obj.cunits = 'm';
            obj.timec = 1e-3*header.dt;
            obj.antennas = 'Seismic Un*x';
            obj.tunits = 'ms';
            obj.comment='true positions';

            obj.timestp = obj.timec*(1:obj.nptsptrc);

            obj.Tx_x = zeros(1,obj.ntrace);
            obj.Tx_y = zeros(1,obj.ntrace);
            obj.Tx_z = zeros(1,obj.ntrace);
            obj.Rx_x = zeros(1,obj.ntrace);
            obj.Rx_y = zeros(1,obj.ntrace);
            obj.Rx_z = zeros(1,obj.ntrace);
            for n=1:obj.ntrace
                obj.Tx_x(n) = tr_header(n).SourceX;
                obj.Tx_y(n) = tr_header(n).SourceY;
                obj.Tx_z(n) = tr_header(n).SourceSurfaceElevation;
                obj.Rx_x(n) = tr_header(n).GroupX;
                obj.Rx_y(n) = tr_header(n).GroupY;
                obj.Rx_z(n) = tr_header(n).ReceiverGroupElevation;
            end

        end
        function readSEGY(obj, basename)
            
            [fields, trueCoord] = MogData.getSEGYfields();
            if isempty(fields)
                error('User aborted')
            end
                        
            s = read_segy(basename, [], fields);

            obj.ntrace = size( s.data, 2 );
            obj.nptsptrc = double(s.bh.hns);
            obj.rstepsz = 0;
            obj.cunits = 'm';
            if s.bh.mfeet == 2
                obj.cunits = 'ft';
            end
            obj.rnomfreq = 0;
            obj.csurvmod = '';
            obj.timec = 0.001 * double(s.bh.hdt);
            obj.rdata = s.data;
            obj.timestp = 0;
            obj.Tx_x = zeros(1,obj.ntrace);
            obj.Tx_y = zeros(1,obj.ntrace);
            obj.Tx_z = zeros(1,obj.ntrace);
            obj.Rx_x = zeros(1,obj.ntrace);
            obj.Rx_y = zeros(1,obj.ntrace);
            obj.Rx_z = zeros(1,obj.ntrace);
            obj.antennas = 'SEG Y';
            obj.tunits = 'ms';
            obj.timestp = obj.timec*(1:obj.nptsptrc);

            % Source Z
            field = fields{3};
            if strcmp(field,'gelev')==1 || strcmp(field,'selev')==1 || ...
                    strcmp(field,'sdepth')==1 || strcmp(field,'gdel')==1 || ...
                    strcmp(field,'sdel')==1 || strcmp(field,'swdep')==1 || ...
                    strcmp(field,'gwdep')==1
                % we must scale 
                ind = s.th.scalel>0;
                if any(ind)
                    eval(['obj.Tx_z(ind) = double(s.th.scalel(ind)) .* double(s.th.',field,'(ind));']);
                end
                ind = ~ind;
                if any(ind)
                    eval(['obj.Tx_z(ind) = double(s.th.',field,'(ind)) ./ abs( double(s.th.scalel(ind)) );']);
                end
            else
                eval(['obj.Tx_z = reshape(double(s.th.',field,'),1,obj.ntrace);']);
            end
            
            % Receiver Z
            field = fields{6};
            if strcmp(field,'gelev')==1 || strcmp(field,'selev')==1 || ...
                    strcmp(field,'sdepth')==1 || strcmp(field,'gdel')==1 || ...
                    strcmp(field,'sdel')==1 || strcmp(field,'swdep')==1 || ...
                    strcmp(field,'gwdep')==1
                % we must scale 
                ind = s.th.scalel>0;
                if any(ind)
                    eval(['obj.Rx_z(ind) = double(s.th.scalel(ind)) .* double(s.th.',field,'(ind));']);
                end
                ind = ~ind;
                if any(ind)
                    eval(['obj.Rx_z(ind) = double(s.th.',field,'(ind)) ./ abs( double(s.th.scalel(ind)) );']);
                end
            else
                eval(['obj.Rx_z = reshape(double(s.th.',field,'),1,obj.ntrace);']);
            end
            
            if trueCoord==true
                obj.comment='true positions';
                
                % Source X
                field = fields{1};
                if strcmp(field,'sx')==1 || strcmp(field,'sy')==1 || ...
                        strcmp(field,'gx')==1 || strcmp(field,'gy')==1 || ...
                        strcmp(field,'xcdp')==1 || strcmp(field,'ycdp')==1
                    % we must scale
                    ind = s.th.scalco>0;
                    if any(ind)
                        eval(['obj.Tx_x(ind) = double(s.th.scalco(ind)) .* double(s.th.',field,'(ind));']);
                    end
                    ind = ~ind;
                    if any(ind)
                        eval(['obj.Tx_x(ind) = double(s.th.',field,'(ind)) ./ abs( double(s.th.scalco(ind)) );']);
                    end
                else
                    eval(['obj.Tx_x = reshape(double(s.th.',field,'),1,obj.ntrace);']);
                end
                
                % Source Y
                field = fields{2};
                if strcmp(field,'sx')==1 || strcmp(field,'sy')==1 || ...
                        strcmp(field,'gx')==1 || strcmp(field,'gy')==1 || ...
                        strcmp(field,'xcdp')==1 || strcmp(field,'ycdp')==1
                    % we must scale
                    ind = s.th.scalco>0;
                    if any(ind)
                        eval(['obj.Tx_y(ind) = double(s.th.scalco(ind)) .* double(s.th.',field,'(ind));']);
                    end
                    ind = ~ind;
                    if any(ind)
                        eval(['obj.Tx_y(ind) = double(s.th.',field,'(ind)) ./ abs( double(s.th.scalco(ind)) );']);
                    end
                else
                    eval(['obj.Tx_y = reshape(double(s.th.',field,'),1,obj.ntrace);']);
                end
                
                % Receiver X
                field = fields{4};
                if strcmp(field,'sx')==1 || strcmp(field,'sy')==1 || ...
                        strcmp(field,'gx')==1 || strcmp(field,'gx')==1 || ...
                        strcmp(field,'xcdp')==1 || strcmp(field,'ycdp')==1
                    % we must scale
                    ind = s.th.scalco>0;
                    if any(ind)
                        eval(['obj.Rx_x(ind) = double(s.th.scalco(ind)) .* double(s.th.',field,'(ind));']);
                    end
                    ind = ~ind;
                    if any(ind)
                        eval(['obj.Rx_x(ind) = double(s.th.',field,'(ind)) ./ abs( double(s.th.scalco(ind)) );']);
                    end
                else
                    eval(['obj.Rx_x = reshape(double(s.th.',field,'),1,obj.ntrace);']);
                end
                
                % Receiver Y
                field = fields{5};
                if strcmp(field,'sx')==1 || strcmp(field,'sy')==1 || ...
                        strcmp(field,'gx')==1 || strcmp(field,'gx')==1 || ...
                        strcmp(field,'xcdp')==1 || strcmp(field,'ycdp')==1
                    % we must scale
                    ind = s.th.scalco>0;
                    if any(ind)
                        eval(['obj.Rx_y(ind) = double(s.th.scalco(ind)) .* double(s.th.',field,'(ind));']);
                    end
                    ind = ~ind;
                    if any(ind)
                        eval(['obj.Rx_y(ind) = double(s.th.',field,'(ind)) ./ abs( double(s.th.scalco(ind)) );']);
                    end
                else
                    eval(['obj.Rx_y = reshape(double(s.th.',field,'),1,obj.ntrace);']);
                end
            end
        end

    end

    methods (Access=private)
        function readRAD(obj, basename)
            fid=fopen(basename,'rt');
            if fid < 2
                fichier = [basename,'.rad'];
                fid=fopen(fichier,'rt');
                if fid < 2
                    fichier = [basename,'.RAD'];
                    fid=fopen(fichier,'rt');
                    if fid < 2
                        error(['Cannot open RAD file: ',basename])
                    end
                end
            end

            while 1
                tline = fgetl(fid);
                if ~ischar(tline), break, end
                if ~isempty(strfind(tline,'SAMPLES:'))
                    obj.nptsptrc = sscanf(tline,'%*8c%d',1);
                elseif ~isempty(strfind(tline,'FREQUENCY:'))
                    obj.timec = sscanf(tline,'%*10c%f',1);
                elseif ~isempty(strfind(tline,'OPERATOR:'))
                    obj.synthetique = ~isempty(regexp(tline,'MoRad', 'once')) || ...
                        ~isempty(regexp(tline,'synthetic', 'once'));
                elseif ~isempty(strfind(tline,'ANTENNAS:'))
                    [debut,fin] = regexp(tline,'\d+');
                    obj.rnomfreq = str2double( tline(debut:fin) );
                    obj.antennas = tline(10:length(tline));
                elseif ~isempty(strfind(tline,'LAST TRACE'))
                    obj.ntrace = sscanf(tline,'%*s %*6c%d',1);
                end
            end
            obj.timec = 1000/obj.timec;
            obj.timestp = obj.timec*(0:obj.nptsptrc-1);
            if obj.synthetique==0
                obj.antennas = strcat(obj.antennas, ' - Ramac');
            end

            fclose(fid);
        end
        function readRD3(obj, basename)
            fid=fopen(basename,'r');
            if fid < 2
                fichier = [basename,'.rd3'];
                fid=fopen(fichier,'r');
                if fid < 2
                    fichier = [basename,'.RD3'];
                    fid=fopen(fichier,'r');
                    if fid < 2
                        error(['Cannot open RD3 file: ',basename])
                    end
                end
            end
            obj.rdata = fread(fid, [obj.nptsptrc obj.ntrace], 'int16','ieee-le');

            fclose(fid);
        end
        function readTLF(obj, basename)
            fid=fopen(basename,'rt');
            if fid < 2
                fichier = [basename,'.tlf'];
                fid=fopen(fichier,'rt');
                if fid < 2
                    fichier = [basename,'.TLF'];
                    fid=fopen(fichier,'rt');
                    if fid < 2
                        
                        button = positionInputdlg(obj);
                        switch button
                            case 'ok'
                                return % params are entered
                            case 'cancel'
                                error(['Cannot open TLF file: ',basename])
                        end
                    end
                end
            end
            fgetl(fid);

            obj.Tx_z = [];
            obj.Rx_z = [];
            while ~feof(fid)
                tnd = fscanf(fid,'%d', 1);
                tnf = fscanf(fid,'%d', 1);
                nt = tnf-tnd+1;

                Rxd = fscanf(fid,'%f', 1);
                Rxf = fscanf(fid,'%f', 1);
                if nt == 1
                    dRx = 1;
                    if Rxd>Rxf, Rxd = Rxf; end
                else
                    dRx = (Rxf-Rxd)/(nt-1);
                end
                Tx  = fscanf(fid,'%f', 1);

                if nt > 0, obj.Tx_z = [obj.Tx_z Tx*ones(1, nt)]; end
                obj.Rx_z = [obj.Rx_z Rxd:dRx:Rxf];
                fgetl(fid);
            end

            fclose(fid);
        end
        function readHD(obj, basename)
            fid=fopen(basename,'rt');
            if fid < 2
                fichier = [basename,'.hd'];
                fid=fopen(fichier,'rt');
                if fid < 2
                    fichier = [basename,'.HD'];
                    fid=fopen(fichier,'rt');
                    if fid < 2
                        error(['Cannot open HD file: ',basename])
                    end
                end
            end

            obj.antennas = 'Ekko';
            obj.synthetique = 0;
            obj.tunits = 'ns';

            fgetl(fid);
            obj.comment = fgetl(fid);
            tline = fgetl(fid);
            tTx_z = 0;
            tRx_z = 0;
            if ~isempty(tline)
                if strcmpi(tline(1:3),'VRP'),
                    %vrp = 1;
                    ic = strfind(tline,'=');
                    count = length(tline);
                    sscanf(tline(ic(1)+1:count),'%f');
                    sscanf(tline(ic(2)+1:count),'%f');
                    tRx_z = sscanf(tline(ic(3)+1:count),'%f');
                elseif strcmpi(tline(1:3),'ZOP'),
                elseif strcmpi(tline(1:3),'MOG'),
                    %mog = 1;
                    ic = strfind(tline,'=');
                    count = length(tline);
                    tTx_z = sscanf(tline(ic(1)+1:count),'%f');
                    tRx_z = sscanf(tline(ic(2)+1:count),'%f');
                end
            end

            obj.date = fgetl(fid);

            while 1
                tline = fgetl(fid);
                if ~ischar(tline), break, end
                if ~isempty(strfind(tline,'NUMBER OF TRACES'))
                    obj.ntrace = sscanf(tline,'%*20c %d',1);
                elseif ~isempty(strfind(tline,'NUMBER OF PTS/TRC'))
                    obj.nptsptrc = sscanf(tline,'%*20c %d',1);
                elseif ~isempty(strfind(tline,'TOTAL TIME WINDOW'))
                    nttwin = sscanf(tline,'%*20c %f',1);
                elseif ~isempty(strfind(tline,'STEP SIZE USED'))
                    obj.rstepsz = sscanf(tline,'%*20c %f',1);
                elseif ~isempty(strfind(tline,'POSITION UNITS'))
                    obj.cunits = sscanf(tline,'%*20c %s',1);
                elseif ~isempty(strfind(tline,'NOMINAL FREQUENCY'))
                    obj.rnomfreq = sscanf(tline,'%*20c %f',1);
                elseif ~isempty(strfind(tline,'ANTENNA SEPARATION'))
                    rantsep = sscanf(tline,'%*20c %f',1);
                elseif ~isempty(strfind(tline,'SURVEY MODE'))
                    obj.csurvmod = tline(21:end);
                end
            end

            fclose(fid);

            obj.timec = nttwin/obj.nptsptrc;
            obj.timestp = obj.timec*(1:obj.nptsptrc);

            obj.Tx_x = zeros(1,obj.ntrace);
            obj.Tx_z = tTx_z*ones(1,obj.ntrace);
            obj.Rx_x = rantsep*ones(1,obj.ntrace);
            obj.Rx_z = tRx_z + obj.rstepsz*(0:(obj.ntrace-1));

        end
        function readDT1(obj, basename)
            fid=fopen(basename,'rt');
            if fid < 2
                fichier = [basename,'.dt1'];
                fid=fopen(fichier,'rt');
                if fid < 2
                    fichier = [basename,'.DT1'];
                    fid=fopen(fichier,'rt');
                    if fid < 2
                        error(['Cannot open DT1 file: ',basename])
                    end
                end
            end
            obj.rdata = [];
            for n=1:nt
                header = fread(fid, 25, 'float32','ieee-le');
                fread(fid, 28, 'char','ieee-le');
                obj.rdata = [obj.rdata fread(fid, header(3), 'int16','ieee-le')];
            end
            fclose(fid);
        end
    end

    methods (Static)
        function [fields, trueCoord] = getSEGYfields()
           
            fnamesSEG = {'tracl'  % 1
            'tracr'
            'fldr'
            'tracf'
            'ep'
            'cdp'
            'cdpt'
            'trid'
            'nvs'
            'nhs'
            'duse'  % 11
            'offset'
            'gelev'
            'selev'
            'sdepth'
            'gdel'
            'sdel'
            'swdep'
            'gwdep'
            'scalel'
            'scalco'  % 21
            'sx'
            'sy'
            'gx'
            'gy'
            'counit'
            'wevel'
            'swevel'
            'sut'
            'gut'
            'sstat' % 31
            'gstat'
            'tstat'
            'laga'
            'lagb'
            'delrt'
            'muts'
            'mute'
            'ns'
            'dt'
            'gain'  % 41
            'igc'
            'igi'
            'corr'
            'sfs'
            'sfe'
            'slen'
            'styp'
            'stas'
            'stae'
            'tatyp'  % 51
            'afilf'
            'afils'
            'nofilf'
            'nofils'
            'lcf'
            'hcf'
            'lcs'
            'hcs'
            'year'
            'day'  % 61
            'hour'
            'minute'
            'sec'
            'timbas'
            'trwf'
            'grnors'
            'grnofr'
            'grnlof'
            'gaps'
            'otrav'  % 71 %  names below arbitrarily given
            'xcdp'   % 72 - X coord of ensemble (CDP) position of this trace
            'ycdp'   % 73 - Y coord of ensemble (CDP) position of this trace
            'ilineno'
            'clineno'
            'shotno'   % 76 - shotpoint number
            'scalsn'
            'tvmunit'   % 78 - trace value measurement units
            'tdcst'
            'tdunit'
            'trid'   % 81 - device/trace identifier
            'scalt'
            'styp'   % 83 - source type/orientation
            'sdir'
            'smeas'
            'smunit'};
        
        wordLengthSEG = [4 4 4 4 4 4 4 2 2 2 2 4 4 4 4 4 4 4 4 2 2 4 4 4 4 2 2 2 2 2 ...
            2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 ...
            2 4 4 4 4 4 2 2 6 2 2 2 2 6 6 2]';
        
        
            % get header fields for coordinates
            trueCoord = true;
            fields = '';
            
            byte1 = cumsum(wordLengthSEG)-wordLengthSEG+1;
            byte2 = cumsum(wordLengthSEG);
            desc = cell(size(fnamesSEG));
            for n=1:numel(desc)
                desc{n} = [num2str(byte1(n),'%3d'),'-',num2str(byte2(n),'%3d'),': ',fnamesSEG{n}];
            end

            fs = 11;
            if ispc
                fs = 9;
            end

            nLines=11;
            vSizeTot = nLines*22 + (nLines+1)*5;

            hSize=350;
            FigPos=[0 0 hSize vSizeTot];
            vSize = 22/vSizeTot;
            vSpace = 5/vSizeTot;
            
            d = dialog('Name','SEGY Trace Header Fields',...
                'Visible','off',...
                'Units','points');
            set(d,'Position',mygetnicedialoglocation(FigPos,get(d,'Units')));
            
            lineNo = 10;
            coordPopup = uicontrol('Style','popupmenu',...
                'FontSize',fs,...
                'String',{'Trace Header Data Contain Abolute Coordinates',...
                'Z Fields Contain Depth Along Borehole'},...
                'Units','normalized',...
                'Position',[0.05 (lineNo-1)*vSize+lineNo*vSpace 0.9 vSize],...
                'Callback',@setTrueCoord,...
                'Parent',d);
            lineNo = 8.5;
            uicontrol('Style','text',...
                'FontSize',fs,...
                'String','Coordinates ',...
                'Units','normalized',...
                'HorizontalAlignment','right',...
                'Position',[0.05 (lineNo-1)*vSize+lineNo*vSpace 0.3 vSize],...
                'Parent',d);
            uicontrol('Style','text',...
                'FontSize',fs,...
                'String','Bytes in Trace Header ',...
                'Units','normalized',...
                'HorizontalAlignment','left',...
                'Position',[0.4 (lineNo-1)*vSize+lineNo*vSpace 0.55 vSize],...
                'Parent',d);
            lineNo = 7.5;
            uicontrol('Style','text',...
                'FontSize',fs,...
                'String','Source X ',...
                'Units','normalized',...
                'HorizontalAlignment','right',...
                'Position',[0.05 (lineNo-1)*vSize+lineNo*vSpace 0.3 vSize],...
                'Parent',d);
            sxPopup = uicontrol('Style','popupmenu',...
                'FontSize',fs,...
                'String',desc,...
                'Value',22,...
                'Units','normalized',...
                'Position',[0.4 (lineNo-1)*vSize+lineNo*vSpace 0.55 vSize],...
                'Callback',@setTrueCoord,...
                'Parent',d);
            lineNo = 6.5;
            uicontrol('Style','text',...
                'FontSize',fs,...
                'String','Source Y ',...
                'Units','normalized',...
                'HorizontalAlignment','right',...
                'Position',[0.05 (lineNo-1)*vSize+lineNo*vSpace 0.3 vSize],...
                'Parent',d);
            syPopup = uicontrol('Style','popupmenu',...
                'FontSize',fs,...
                'String',desc,...
                'Value',23,...
                'Units','normalized',...
                'Position',[0.4 (lineNo-1)*vSize+lineNo*vSpace 0.55 vSize],...
                'Callback',@setTrueCoord,...
                'Parent',d);
            lineNo = 5.5;
            uicontrol('Style','text',...
                'FontSize',fs,...
                'String','Source Z ',...
                'Units','normalized',...
                'HorizontalAlignment','right',...
                'Position',[0.05 (lineNo-1)*vSize+lineNo*vSpace 0.3 vSize],...
                'Parent',d);
            szPopup = uicontrol('Style','popupmenu',...
                'FontSize',fs,...
                'String',desc,...
                'Value',14,...
                'Units','normalized',...
                'Position',[0.4 (lineNo-1)*vSize+lineNo*vSpace 0.55 vSize],...
                'Callback',@setTrueCoord,...
                'Parent',d);

            lineNo = 4.5;
            uicontrol('Style','text',...
                'FontSize',fs,...
                'String','Receiver X ',...
                'Units','normalized',...
                'HorizontalAlignment','right',...
                'Position',[0.05 (lineNo-1)*vSize+lineNo*vSpace 0.3 vSize],...
                'Parent',d);
            gxPopup = uicontrol('Style','popupmenu',...
                'FontSize',fs,...
                'String',desc,...
                'Value',24,...
                'Units','normalized',...
                'Position',[0.4 (lineNo-1)*vSize+lineNo*vSpace 0.55 vSize],...
                'Callback',@setTrueCoord,...
                'Parent',d);
            lineNo = 3.5;
            uicontrol('Style','text',...
                'FontSize',fs,...
                'String','Receiver Y ',...
                'Units','normalized',...
                'HorizontalAlignment','right',...
                'Position',[0.05 (lineNo-1)*vSize+lineNo*vSpace 0.3 vSize],...
                'Parent',d);
            gyPopup = uicontrol('Style','popupmenu',...
                'FontSize',fs,...
                'String',desc,...
                'Value',25,...
                'Units','normalized',...
                'Position',[0.4 (lineNo-1)*vSize+lineNo*vSpace 0.55 vSize],...
                'Callback',@setTrueCoord,...
                'Parent',d);
            lineNo = 2.5;
            uicontrol('Style','text',...
                'FontSize',fs,...
                'String','Receiver Z ',...
                'Units','normalized',...
                'HorizontalAlignment','right',...
                'Position',[0.05 (lineNo-1)*vSize+lineNo*vSpace 0.3 vSize],...
                'Parent',d);
            gzPopup = uicontrol('Style','popupmenu',...
                'FontSize',fs,...
                'String',desc,...
                'Value',13,...
                'Units','normalized',...
                'Position',[0.4 (lineNo-1)*vSize+lineNo*vSpace 0.55 vSize],...
                'Callback',@setTrueCoord,...
                'Parent',d);

            lineNo = 1.25;
            uicontrol('Style','pushbutton',...
                'FontSize',fs,...
                'String','Cancel',...
                'Units','normalized',...
                'Position',[0.1 (lineNo-1)*vSize+lineNo*vSpace 0.35 vSize],...
                'Callback','delete(gcf)',...
                'Parent',d);
            uicontrol('Style','pushbutton',...
                'FontSize',fs,...
                'String','Done',...
                'Units','normalized',...
                'Position',[0.55 (lineNo-1)*vSize+lineNo*vSpace 0.35 vSize],...
                'Callback',@done,...
                'Parent',d);
            
            d.Visible = 'on';
            uiwait(d)
            
            function setTrueCoord(varargin)
                if coordPopup.Value == 2
                    sxPopup.Enable = 'off';
                    syPopup.Enable = 'off';
                    gxPopup.Enable = 'off';
                    gyPopup.Enable = 'off';
                else
                    sxPopup.Enable = 'on';
                    syPopup.Enable = 'on';
                    gxPopup.Enable = 'on';
                    gyPopup.Enable = 'on';
                end   
            end
            function done(varargin)
                if coordPopup.Value==2
                    trueCoord = false;
                else
                    trueCoord = true;
                end
                fields = {fnamesSEG{sxPopup.Value},...
                    fnamesSEG{syPopup.Value},...
                    fnamesSEG{szPopup.Value},...
                    fnamesSEG{gxPopup.Value},...
                    fnamesSEG{gyPopup.Value},...
                    fnamesSEG{gzPopup.Value},...
                    'scalel','scalco'};
                
                delete(d)
            end
        end
    end
end
