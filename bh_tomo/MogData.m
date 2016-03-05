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
            s = read_segy(basename, [], 'scalco,sx,sy,gx,gy,selev,gelev');

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
            obj.comment='true positions';
            obj.tunits = 'ms';

            obj.timestp = obj.timec*(1:obj.nptsptrc);

            ind = s.th.scalco>0;
            if any(ind)
                obj.Tx_x(ind) = double(s.th.scalco(ind)) .* double(s.th.sx(ind));
                obj.Tx_y(ind) = double(s.th.scalco(ind)) .* double(s.th.sy(ind));
                obj.Tx_z(ind) = double(s.th.scalco(ind)) .* double(s.th.selev(ind));
                obj.Rx_x(ind) = double(s.th.scalco(ind)) .* double(s.th.gx(ind));
                obj.Rx_y(ind) = double(s.th.scalco(ind)) .* double(s.th.gy(ind));
                obj.Rx_z(ind) = double(s.th.scalco(ind)) .* double(s.th.gelev(ind));
            end

            ind = ~ind;
            if any(ind)
                obj.Tx_x(ind) = double(s.th.sx(ind)) ./ abs( double(s.th.scalco(ind)) );
                obj.Tx_y(ind) = double(s.th.sy(ind)) ./ abs( double(s.th.scalco(ind)) );
                obj.Tx_z(ind) = double(s.th.selev(ind)) ./ abs( double(s.th.scalco(ind)) );
                obj.Rx_x(ind) = double(s.th.gx(ind)) ./ abs( double(s.th.scalco(ind)) );
                obj.Rx_y(ind) = double(s.th.gy(ind)) ./ abs( double(s.th.scalco(ind)) );
                obj.Rx_z(ind) = double(s.th.gelev(ind)) ./ abs( double(s.th.scalco(ind)) );
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
                        error(['Cannot open TLF file: ',basename])
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
end
