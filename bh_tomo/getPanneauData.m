function [data,ind] = getPanneauData(panneau,db_file,type,selected_mogs, varargin)  %YH
% function [data,ind] = getPanneauData(panneau,db_file,type)
%

% Copyright (C) 2005 Bernard Giroux
%
% This program is free software; you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation; either version 3 of the License, or
% (at your option) any later version.
%
% This program is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details.
%
% You should have received a copy of the GNU General Public License
% along with this program.  If not, see <http://www.gnu.org/licenses/>.
% 
%

type2 = '';
if nargin>=5   %4
    vlim = varargin{1};
else
    vlim = [];
end
if nargin>=6  %5
    type2 = varargin{2};
end

%data = [];
load(db_file,'mogs','air')

tt = [];
et = [];
%fait = [];
t0 = [];
in = [];
ind = [];
data = [];
switch lower(type)
	%replace all mogs(n) in this case by mogs(selected_mogs(n)) 2012-10-16 YH
	case 'tt'
	 for n=1:length(panneau.mogs)
		 fac_dt = 1;
		 if isempty( find(n==selected_mogs, 1) )
			 ind = [ind false(size(mogs(panneau.mogs(n)).tt))];
			 t0 = [t0 zeros(size(mogs(panneau.mogs(n)).tt))];
		 else
			 if strcmp( mogs(panneau.mogs(n)).data.tunits,'ns' )==1
				 if ~isempty(mogs(panneau.mogs(n)).ttTx)
					 t0b = mogs(panneau.mogs(n)).ttTx;
				 else
					 av = air( mogs(panneau.mogs(n)).av );
					 ap = air( mogs(panneau.mogs(n)).ap );
					 if mogs(panneau.mogs(n)).data.synthetique==1
						 t0b = zeros(1,length(mogs(panneau.mogs(n)).tt));
						 fac_dt_av = 1;
						 fac_dt_ap = 1;
					 else
						 [t0b, fac_dt_av, fac_dt_ap] = corr_t0(length(mogs(panneau.mogs(n)).tt), av, ap);
					 end
					 if mogs(panneau.mogs(n)).user_fac_dt == 1
						 fac_dt = mogs(panneau.mogs(n)).fac_dt;
					 else
						 if fac_dt_av~=1 && fac_dt_ap ~= 1
							 fac_dt = 0.5*(fac_dt_av+fac_dt_ap);
						 elseif fac_dt_av~=1
							 fac_dt = fac_dt_av;
						 elseif fac_dt_ap~=1
							 fac_dt = fac_dt_ap;
						 end
					 end
				 end
				 t0 = [t0 fac_dt*t0b];
			 end
			 ind = [ind mogs(panneau.mogs(n)).tt~=-1];
		 end
		 tt = [tt fac_dt*mogs(panneau.mogs(n)).tt];
		 et = [et fac_dt*mogs(panneau.mogs(n)).et*mogs(panneau.mogs(n)).f_et];
		 %fait = [fait mogs(panneau.mogs(n)).tt_fait];
		 in = [in mogs(panneau.mogs(n)).in];
	 end
 case 'amp'
	for n=1:length(panneau.mogs) 
		if isempty( find(n==selected_mogs, 1) )
			ind = [ind false(size(mogs(panneau.mogs(n)).tt))];
		else
			ind = [ind mogs(panneau.mogs(n)).tauApp~=-1];
		end
		tt = [tt mogs(panneau.mogs(n)).tauApp];
		et = [et mogs(panneau.mogs(n)).tauApp_et*mogs(panneau.mogs(n)).f_et];
		%fait = [fait mogs(panneau.mogs(n)).amp_fait];
		in = [in mogs(panneau.mogs(n)).in];
	end
 case 'fce'
	for n=1:length(panneau.mogs) 
		if isempty( find(n==selected_mogs, 1) )
			ind = [ind false(size(mogs(panneau.mogs(n)).tt))];
		else
			ind = [ind mogs(panneau.mogs(n)).tauFce~=-1];
		end
		tt = [tt mogs(panneau.mogs(n)).tauFce];
		et = [et mogs(panneau.mogs(n)).tauFce_et*mogs(panneau.mogs(n)).f_et];
		%fait = [fait mogs(panneau.mogs(n)).amp_fait];
		in = [in mogs(panneau.mogs(n)).in];
	end
 case 'hyb'
	for n=1:length(panneau.mogs) 
		if isempty( find(n==selected_mogs, 1) )
			ind = [ind false(size(mogs(panneau.mogs(n)).tt))];
		else
			ind = [ind mogs(panneau.mogs(n)).tauHyb~=-1];
		end
		tt = [tt mogs(panneau.mogs(n)).tauHyb];
		et = [et mogs(panneau.mogs(n)).tauHyb_et*mogs(panneau.mogs(n)).f_et];
		%fait = [fait mogs(panneau.mogs(n)).amp_fait];
		in = [in mogs(panneau.mogs(n)).in];
	end
 case 'ant'
	[~, ind] = getPanneauData(panneau,db_file,'tt',selected_mogs);
	load(db_file,'forages')
	for n=1:length(panneau.mogs)
		mog = mogs(panneau.mogs(n));
		tt = [tt forages(mog.Tx).diam*ones(size(mog.tt))];  % contient diametre Tx
		et = [et forages(mog.Rx).diam*ones(size(mog.tt))];  % contient diametre Rx
	end
	in = ind;
 case 'depth'
	if strcmp(type2,'')
		return
	end
	[~, ind] = getPanneauData(panneau,db_file,type2,selected_mogs);
	for n=1:length(panneau.mogs)
		tt = [tt mogs(panneau.mogs(n)).Tx_z_orig];
		et = [et mogs(panneau.mogs(n)).Rx_z_orig];
		in = [in mogs(panneau.mogs(n)).in];
	end        
end
no = [];
for n=1:length(panneau.mogs)
	no = [no mogs(panneau.mogs(n)).no_traces];
end

ind = ind & in;
if ~isempty(t0)
	tt = tt-t0;
end
if ~isempty(vlim)
	%%%%%%%%%%%%%%YH
	if isfield(panneau,'grid3d')
		l = sqrt(sum((panneau.grid3d.Tx-panneau.grid3d.Rx).^2,2))';
	else
    %%%%%%%%%%%%%%%%%%%%
		l = sqrt(sum((panneau.grid.Tx-panneau.grid.Rx).^2,2))';
	end
	vapp = l./tt;
	in2 = vapp<vlim;
	disp([num2str(sum(~in2&ind)),' rays with apparent velocity above ',num2str(vlim)])
	ind = ind & in2;
end

data = [tt(ind)' et(ind)' no(ind)'];
