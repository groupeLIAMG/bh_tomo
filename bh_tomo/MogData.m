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
        Tx_z                  % z coordinate of source
        Rx_x                  % x coordinate of receiver
        Rx_z                  % z coordinate of receiver
        antennas              % type of antenna
        synthetique           % true of synthetic data
        tunits                % time units
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
            obj.Tx_z = 0;
            obj.Rx_x = 0;
            obj.Rx_z = 0;
            obj.antennas = '';
            obj.synthetique = 0;
            obj.tunits = '';
        end
        function readRAMAC(obj, basename )
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
    end
end

