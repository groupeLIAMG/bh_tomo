classdef Model < handle
    %UNTITLED Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        name
        mogs         %
        boreholes    % list of boreholes
        grid
        tt_covar     % covariance model - travel time
        amp_covar    % covariance model - amplitudes
        inv_res      % inversion results
    end
    
    methods
        function obj = Model(n)
            obj.name = n;
            obj.grid = Grid.empty;
            obj.tt_covar = CovarianceModel();
            obj.amp_covar = CovarianceModel();
        end
        function set.name(obj, n)
            if ischar(n)
                obj.name = n;
            else
                error('Model name must be a string')
            end
        end
    end
    methods (Static=true)
        function [data,ind] = getModelData(model,db_file,type,selected_mogs,varargin)
            data = [];
            type2 = '';
            if nargin>=5   %4
                vlim = varargin{1};
            else
                vlim = [];
            end
            if nargin>=6  %5
                type2 = varargin{2};
            end

            load(db_file,'mogs','air')
            
            tt = [];
            et = [];
            in = [];
            ind = [];
            switch lower(type)
                %
                case 'tt'
                    for n=1:length(model.mogs)
                        fac_dt = 1;
                        if isempty( find(n==selected_mogs, 1) )
                            ind = [ind false(size(mogs(model.mogs(n)).tt))]; %#ok<*AGROW>
                        else
                            ind = [ind mogs(model.mogs(n)).tt~=-1];
                        end
                        tt = [tt mogs(model.mogs(n)).getCorrectedTravelTimes(air)];
                        et = [et fac_dt*mogs(model.mogs(n)).et*mogs(model.mogs(n)).f_et];
                        in = [in mogs(model.mogs(n)).in];
                    end
                case 'amp'
                    for n=1:length(model.mogs)
                        if isempty( find(n==selected_mogs, 1) )
                            ind = [ind false(size(mogs(model.mogs(n)).tt))];
                        else
                            ind = [ind mogs(model.mogs(n)).tauApp~=-1];
                        end
                        tt = [tt mogs(model.mogs(n)).tauApp];
                        et = [et mogs(model.mogs(n)).tauApp_et*mogs(model.mogs(n)).f_et];
                        in = [in mogs(model.mogs(n)).in];
                    end
                case 'fce'
                    for n=1:length(model.mogs)
                        if isempty( find(n==selected_mogs, 1) )
                            ind = [ind false(size(mogs(model.mogs(n)).tt))];
                        else
                            ind = [ind mogs(model.mogs(n)).tauFce~=-1];
                        end
                        tt = [tt mogs(model.mogs(n)).tauFce];
                        et = [et mogs(model.mogs(n)).tauFce_et*mogs(model.mogs(n)).f_et];
                        in = [in mogs(model.mogs(n)).in];
                    end
                case 'hyb'
                    for n=1:length(model.mogs)
                        if isempty( find(n==selected_mogs, 1) )
                            ind = [ind false(size(mogs(model.mogs(n)).tt))];
                        else
                            ind = [ind mogs(model.mogs(n)).tauHyb~=-1];
                        end
                        tt = [tt mogs(model.mogs(n)).tauHyb];
                        et = [et mogs(model.mogs(n)).tauHyb_et*mogs(model.mogs(n)).f_et];
                        in = [in mogs(model.mogs(n)).in];
                    end
                case 'ant'
                    [~, ind] = getPanneauData(model,db_file,'tt',selected_mogs);
                    load(db_file,'boreholes')
                    for n=1:length(model.mogs)
                        mog = mogs(model.mogs(n));
                        tt = [tt boreholes(mog.Tx).diam*ones(size(mog.tt))];  % contient diametre Tx
                        et = [et boreholes(mog.Rx).diam*ones(size(mog.tt))];  % contient diametre Rx
                    end
                    in = ind;
                case 'depth'
                    if strcmp(type2,'')
                        return
                    end
                    [~, ind] = Model.getModelData(model,db_file,type2,selected_mogs);
                    for n=1:length(model.mogs)
                        tt = [tt mogs(model.mogs(n)).Tx_z_orig];
                        et = [et mogs(model.mogs(n)).Rx_z_orig];
                        in = [in mogs(model.mogs(n)).in];
                    end
            end
            no = [];
            for n=1:length(model.mogs)
                no = [no mogs(model.mogs(n)).no_traces];
            end
            
            ind = ind & in;
            if ~isempty(vlim)
                
                l = sqrt(sum((model.grid.Tx-model.grid.Rx).^2,2))';
                vapp = l./tt;
                in2 = vapp<vlim;
                disp([num2str(sum(~in2&ind)),' rays with apparent velocity above ',num2str(vlim)])
                ind = ind & in2;
            end
            
            data = [tt(ind)' et(ind)' no(ind)'];

        end
    end
end

