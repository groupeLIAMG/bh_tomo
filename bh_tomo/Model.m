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
        tlinv_res    % time-lapse inversion results
    end
    
    methods
        function obj = Model(n)
            if isstring(n) || ischar(n)
                obj.name = n;
                obj.grid = Grid.empty;
                obj.tt_covar = CovarianceModel();
                obj.amp_covar = CovarianceModel();
            elseif isstruct(n)
                if ~strcmp(n.type, 'normal')
                    error('Super panel not handled')
                end
                if isfield(n, 'name')
                    obj.name = n.name;
                elseif isfield(n, 'nom')
                    obj.name = n.nom;
                else
                    error('Invalid input')
                end
                obj.mogs = n.mogs;
                if isfield(n, 'boreholes')
                    obj.boreholes = n.boreholes;
                elseif isfield(n, 'forages')
                    obj.boreholes = n.forages;
                end
                obj.grid = Grid2D(n.grid.grx, n.grid.grz);
                obj.grid.cont = Constraints(n.grid.cont);
                obj.grid.Tx = n.grid.Tx;
                obj.grid.Rx = n.grid.Rx;
                obj.grid.TxCosDir = n.grid.TxCosDir;
                obj.grid.RxCosDir = n.grid.RxCosDir;
                obj.grid.bord = n.grid.bord;
                if isfield(n.grid, 'Tx_Z_water')
                    obj.grid.Tx_Z_water = n.grid.Tx_Z_water;
                elseif isfield(n.grid, 'Tx_Z_eau')
                    obj.grid.Tx_Z_water = n.grid.Tx_Z_eau;
                end
                if isfield(n.grid, 'Rx_Z_water')
                    obj.grid.Rx_Z_water = n.grid.Rx_Z_water;
                elseif isfield(n.grid, 'Rx_Z_eau')
                    obj.grid.Rx_Z_water = n.grid.Rx_Z_eau;
                end
                obj.grid.in = n.grid.in;
                obj.grid.flip = n.grid.flip;
                obj.grid.type = '2D';
                if isfield(n.grid, 'borehole_x0')
                    obj.grid.borehole_x0 = n.grid.borehole_x0;
                elseif isfield(n.grid, 'forage_x0')
                    obj.grid.borehole_x0 = n.grid.forage_x0;
                end
                obj.grid.x0 = n.grid.x0;
                if isfield(n, 'tt_covar')
                    obj.tt_covar = CovarianceModel(n.tt_covar);
                end
                if isfield(n, 'amp_covar')
                    obj.amp_covar = CovarianceModel();
                end
                
            else
                error('Invalid input')
            end
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
        function [data,ind,delta] = getModelData(model,db_file,type,selected_mogs,varargin)
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
            delta = false;
            
            check = false(1, numel(selected_mogs));
            for n=1:numel(selected_mogs)
                check(n) = mogs(model.mogs(selected_mogs(n))).delta;
            end
            if all(check)
                delta = true;
            elseif any(check)
                errordlg('Parameter delta of selected mogs are not consistent')
                return
            end
            
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
                    [~, ind, ~] = Model.getModelData(model,db_file,type2,selected_mogs);
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

