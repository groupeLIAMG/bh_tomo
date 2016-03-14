function updateCoords(obj,varargin)

updateAll = false;
if nargin >=2
	if islogical( varargin{1} )
		updateAll = varargin{1};
	else
		nos = varargin{1};
	end
end
if exist('nos','var')==0
	if updateAll
		nos = 1:length(obj.mogs);
	else
		nos = get(obj.handles.listMogs,'Value');
	end
end

boreholes = obj.bhUI.boreholes;

for no=nos
    
    mog = obj.mogs(no);
    
	if strcmp(mog.data.comment, 'true positions')
        Tx = [mog.data.Tx_x(:) mog.data.Tx_y(:) mog.data.Tx_z(:)];
        mog.TxCosDir = zeros(size(Tx));
        tmp = unique(Tx,'rows');
        tmp = sort(tmp,1,'descend');
        v = -diff(tmp);
        d = sqrt(sum(v.^2,2));
        l = v./kron(d,[1 1 1]);
        l = [l; l(end,:)]; %#ok<AGROW>
        
        for n=1:size(tmp,1)
            ind = Tx(:,1)==tmp(n,1) & Tx(:,2)==tmp(n,2) & Tx(:,3)==tmp(n,3);
            mog.TxCosDir(ind,1) = l(n,1);
            mog.TxCosDir(ind,2) = l(n,2);
            mog.TxCosDir(ind,3) = l(n,3);
        end
        
        Rx = [mog.data.Rx_x(:) mog.data.Rx_y(:) mog.data.Rx_z(:)];
        mog.RxCosDir = zeros(size(Rx));
        tmp = unique(Rx,'rows');
        tmp = sort(tmp,1,'descend');
        v = -diff(tmp);
        d = sqrt(sum(v.^2,2));
        l = v./kron(d,[1 1 1]);
        l = [l; l(end,:)]; %#ok<AGROW>
        
        for n=1:size(tmp,1)
            ind = Tx(:,1)==tmp(n,1) & Tx(:,2)==tmp(n,2) & Tx(:,3)==tmp(n,3);
            mog.RxCosDir(ind,1) = l(n,1);
            mog.RxCosDir(ind,2) = l(n,2);
            mog.RxCosDir(ind,3) = l(n,3);
        end
        
		continue
	end
	if isempty( mog.Tx ) || isempty( mog.Rx )
		if length(nos)==1
			return
		else
			continue
		end
	end

	if mog.Tx == mog.Rx
		if length(nos)==1
			uiwait(warndlg(['Tx et Rx are in the same well: ', boreholes(mog.Rx).name]))
			return
		else
			continue
		end
	end

	if mog.type == 1 % cross hole
		mog.data.csurvmod = 'SURVEY MODE        = Trans. - MOG';
		% Tx
		if abs(boreholes(mog.Tx).X-boreholes(mog.Tx).Xmax)<1.0e-5 && ...
				abs(boreholes(mog.Tx).Y-boreholes(mog.Tx).Ymax)<1.0e-5
			% forage vertical
			mog.data.Tx_x(:) = boreholes( mog.Tx ).X;
			mog.data.Tx_y(:) = boreholes( mog.Tx ).Y;
			mog.data.Tx_z = boreholes( mog.Tx ).Z - mog.data.TxOffset - ...
				mog.Tx_z_orig;
            mog.TxCosDir = repmat([0 0 1],mog.data.ntrace,1);
			boreholes( mog.Tx ).fdata = [boreholes( mog.Tx ).X boreholes( mog.Tx ).Y boreholes( mog.Tx ).Z;
				boreholes( mog.Tx ).Xmax boreholes( mog.Tx ).Ymax boreholes( mog.Tx ).Zmax];
		else
			[mog.data.Tx_x, mog.data.Tx_y, mog.data.Tx_z, mog.TxCosDir] = ...
				obj.projectBorehole(boreholes(mog.Tx).fdata, ...
				mog.Tx_z_orig+mog.data.TxOffset, ['Tx - ',boreholes(mog.Tx).name]);
		end
		% Rx
		if abs(boreholes(mog.Rx).X-boreholes(mog.Rx).Xmax)<1.0e-5 && ...
				abs(boreholes(mog.Rx).Y-boreholes(mog.Rx).Ymax)<1.0e-5
			mog.data.Rx_x(:) = boreholes( mog.Rx ).X;
			mog.data.Rx_y(:) = boreholes( mog.Rx ).Y;
			mog.data.Rx_z = boreholes( mog.Rx ).Z - mog.data.RxOffset - ...
				mog.Rx_z_orig;
            mog.RxCosDir = repmat([0 0 1],mog.data.ntrace,1);
			boreholes( mog.Rx ).fdata = [boreholes( mog.Rx ).X boreholes( mog.Rx ).Y boreholes( mog.Rx ).Z;
				boreholes( mog.Rx ).Xmax boreholes( mog.Rx ).Ymax boreholes( mog.Rx ).Zmax];
		else
			[mog.data.Rx_x, mog.data.Rx_y, mog.data.Rx_z, mog.RxCosDir] = ...
				obj.projectBorehole(boreholes(mog.Rx).fdata, ...
				mog.Rx_z_orig+mog.data.RxOffset, ['Rx - ',boreholes(mog.Rx).name]);
		end
	else % VRP
		mog.data.csurvmod = 'SURVEY MODE        = Trans. - VRP';
        % Rx
		if abs(boreholes(mog.Rx).X-boreholes(mog.Rx).Xmax)<1.0e-5 || ...
                abs(boreholes(mog.Rx).Y-boreholes(mog.Rx).Ymax)<1.0e-5 %#ok<ALIGN>
    		mog.data.Rx_x(:) = boreholes( mog.Rx ).X;
        	mog.data.Rx_y(:) = boreholes( mog.Rx ).Y;
            mog.data.Rx_z = boreholes( mog.Rx ).Z - mog.data.RxOffset - ...
                mog.Rx_z_orig;
            mog.RxCosDir = repmat([0 0 1],mog.data.ntrace,1);
			boreholes( mog.Rx ).fdata = [boreholes( mog.Rx ).X boreholes( mog.Rx ).Y boreholes( mog.Rx ).Z;
                boreholes( mog.Rx ).Xmax boreholes( mog.Rx ).Ymax boreholes( mog.Rx ).Zmax];
		else
			[mog.data.Rx_x, mog.data.Rx_y, mog.data.Rx_z, mog.RxCosDir] = ...
				obj.projectBorehole(boreholes(mog.Rx).fdata, ...
				mog.Rx_z_orig+mog.data.RxOffset, ['Rx - ',boreholes(mog.Rx).name]);
        end
        % Tx en surface
		theta = atan2( boreholes( mog.Tx ).Y-boreholes( mog.Rx ).Y,...
			boreholes( mog.Tx ).X-boreholes( mog.Rx ).X );
		mog.data.Tx_x = boreholes( mog.Rx ).X + mog.Tx_z_orig*cos(theta);
		mog.data.Tx_y = boreholes( mog.Rx ).Y + mog.Tx_z_orig*sin(theta);
		% z -> on assume que z varie lineairement entre les deux trous
		l = sqrt( (boreholes( mog.Tx ).Y-boreholes( mog.Rx ).Y)^2 + ...
			(boreholes( mog.Tx ).X-boreholes( mog.Rx ).X)^2 );
		dz = boreholes( mog.Tx ).Z_surf-boreholes( mog.Rx ).Z_surf;
		mog.data.Tx_z = boreholes( mog.Rx ).Z_surf + dz*mog.Tx_z_orig/l;
        
        boreholes( mog.Tx ).fdata = [mog.data.Tx_x(1) mog.data.Tx_y(1) mog.data.Tx_z(1); ...
            mog.data.Tx_x(end) mog.data.Tx_y(end) mog.data.Tx_z(end)];
        
        d = sqrt(sum((boreholes( mog.Tx ).fdata(2,1:3)-boreholes( mog.Tx ).fdata(1,1:3)).^2));
        % cosinus directeurs
        l = (boreholes( mog.Tx ).fdata(2,1:3)-boreholes( mog.Tx ).fdata(1,1:3))./d;
        mog.TxCosDir = repmat(l,mog.data.ntrace,1);
	end
end


end