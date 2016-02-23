function [stdp,SNRD]=calcul_stdp(x,C,pr,Fa,delta,RSNR)

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


% picking standard deviation calculation

[N,M]=size(x);

%nb_p4=round(1/(4*delta*(freq*10^6))); % nb d'echantillons sur 1/4 de période
stdp=zeros(M,1); pos=stdp; SNRD=stdp; pr2lob=stdp;
prc = pr;
h = waitbar(0,'picking standard deviation calculation ...Please wait ...');

inc_wb = round(M/100);

for i=1:M

	nb_p=round(1/(delta*Fa(i)));%nb d'echantillons sur 1 période
	std_noise=std(x(N-100:end,i));
	xp=abs(x(pr(i),i));

	% replacer pr(i) si le picking depasse le premier lobe

	if RSNR(i) < 20 && xp > 2*std_noise
		pr(i)= pr(i)-round(nb_p/3);
		if pr(i)< 1;
			pr(i)=1;
		end
	end

	upper = pr(i)+2*nb_p;
	upper(upper>N) = N;
	x1=x(pr(i):upper,i);

	% picking correction
	[indmin,indmax] = extr(x1');
	if isempty(indmin)==0 && indmin(1) < round(nb_p/4) && RSNR(i) > 20  %| indmax(1) < round(nb_p/4)
		corr = indmin(1);
		prc(i)=pr(i)+ corr;
	end
	%--------------------
	
	upper = pr(i)+4*nb_p;
	upper(upper>N) = N;
	
	% localization of second lobe
	x2=angle(C(pr(i)+2:upper,i)); % wavelet phase to detect the second lobe
	[indmin2,indmax2] = extr(x2');
	ind2lob= indmax2(x2(indmax2) > 2.4);
	ind2lob=ind2lob(1);

	% position of first lobe
	%pos1=sort([indmin(:);indmax(:)]);
	%pos1(pos1 < round(nb_p/4))=[];
	pos1=min(indmax);
	if isempty(pos1)
        pos1=1;
    end
    pos(i)=pos1;
	xpos1=x1(pos1);
	%--------------
	pos_std=find(x1(1:pos1)<xpos1/3);

	std_sig=abs(xpos1-x(pr(i),i));
	if std_noise ~= 0
		SNRD(i)= std_sig/std_noise;
	else
		SNRD(i) = inf;
	end

	pr2lob(i)= prc(i)+ 2 + ind2lob;

	if SNRD(i) < 10 && RSNR(i) < 20 || SNRD(i) < 10 && isempty(pos_std)==1 && RSNR(i) < 20
		%prc(i)= prc(i)- 2 + ind2lob - 1.5 * nb_p;
		stdp(i)=round(nb_p/2);
	elseif isempty(pos_std) | pos_std < 3
		stdp(i)=pos1;
	else
		stdp(i)=pos_std(end);
	end
	
	if rem(i, inc_wb), waitbar(i/M, h), end
end

close(h)



