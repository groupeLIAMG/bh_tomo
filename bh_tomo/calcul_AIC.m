function [p,p_loc,AIC]=calcul_AIC(x,ind_max,Fa,delta,SNR)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % radar onset picking                                                   %
%  ---------------------------------------------------                    %
% ======================================================================= %
% Copyright (C) 2008 Abderrezak BOUCHEDDA                                 %
%=====oooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooo====%
%          contact: ---------->\\\////                                    %
%                               |_ _|                                     %
%                               (@ @)                                     %
%               **********oooO***(_)***Oooo**********                     %
%               * -----> Abderrezak BOUCHEDDA<----- *                     %
%               *          Bernard Giroux           *                     %
%               *          giroux@geo.polymtl.ca    *                     %
%               *      bouchedda@geo.polymtl.ca     *                     %
%               *  Ecole polytechnique de Montreal  *                     %
%               *  depart. de géophysique appliquée *                     %
%               *  http://geo.polymtl.ca            *                     %
%               *************************************                     %
%                             |_______|                                   %
%                              |__|__|                                    %
%                               () ()                                     %
%                              ooO Ooo                                    %
% This program is free software; you can redistribute it and/or modify    %
% it under the terms of the GNU General Public License as published by    %
% the Free Software Foundation; either version 3 of the License, or       %
% (at your option) any later version.                                     %
%                                                                         %
% This program is distributed in the hope that it will be useful,         %
% but WITHOUT ANY WARRANTY; without even the implied warranty of          %
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the           %
% GNU General Public License for more details.                            %
%                                                                         %
% You should have received a copy of the GNU General Public License       %
% along with this program.  If not, see <http://www.gnu.org/licenses/>.   %
%                                                                         %
%                                                                         %
%                                                                         %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% p : AIC onset picking



[N,M]=size(x);


%warning off
p=zeros(M,1);
p_loc=p;
h_AIC = waitbar(0,'AIC calculation ...Please wait ...');

inc_wb = round(M/100);
warning off
for i=1:M
	nb_p15=round(1/(Fa(i)*delta)); % number of samples of 1.5 period
    x11=x(1:ind_max(i),i);
	x1=x11;
	N=length(x1);
	if N<10
		p(i) = 1;
		continue
	end
	
	AIC=zeros(size(x1));

	for k=2:N

		AIC(k) = k*log(var(x1(1:k),0,1))+(N-k-1)*log(var(x1(k+1:N),0,1));

	end
	AIC(1)=AIC(2);
	AIC(isinf(AIC)==1)=NaN; % enlever les valeurs inf
	[I,p1]=min(AIC);

	if SNR(i) > 20
		[indmin] = extr(AIC');
		p2=p1-nb_p15; % on cherche les minimas locaux sur 1 periode en arriere.

		if p2 < 1
			p2=1;
		end

		%indmin(AIC(indmin) > AIC(p1)+4)=[];
		ind1=find(indmin >= p2 & indmin <= p1);

		if isempty(ind1)==1
			p(i)=p1(1);
		else
			p_loc(i)=p1(1)-indmin(ind1(1));
            indmax = indmin(ind1(1))+round(0.25*nb_p15);
            if indmax > numel(AIC)
                indmax = numel(AIC);
            end
			diff = AIC(indmin(ind1(1)))-AIC(indmax);
			if  diff < 0 ; % eviter les minimas locaux dus aux bruit transitoires
				p(i)=indmin(ind1(1));
			else
				p(i)=p1(1);
			end
		end

	else
		p(i)=p1(1);
	end
	if rem(i, inc_wb)==0, waitbar(i/M, h_AIC), end
end
close(h_AIC)

warning on