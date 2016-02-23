function [pc]=AIC_correction(x,y,ind_max,p,Fa,delta,SNR)

% x: complex wavelet coeff
% y : denoised data
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

%warning off



[N,M]=size(x);
K = 2^floor(log2(N));
pc=zeros(M,1);% apc=pc;pz=pc;
h_AIC = waitbar(0,'AIC correction calculation ...Please wait ...');

inc_wb = round(M/100);
for i=1:M

    nb_p  = round(1/(delta*Fa(i)));
    nb_p1 = round(1/(3*delta*Fa(i)));
    nb_p2 = round(1.3*nb_p);
    stdn= std(y(K-nb_p:K,i)); % std of noise

    % phase calculation
    maxi=ind_max(i)+30;
    if maxi > N
        maxi=N;
    end

    phi=angle(x(1:maxi,i));
    inc1=p(i)-round(nb_p/4);
    inc2=p(i)+round(nb_p/4);

    if inc1 < 1
        inc1=1;
        disp('warning : onset picking close to start time')
    end

    if inc2 > maxi
        inc2=maxi;
    end

    a=phi(inc1:inc2);
    a=a(:)';
    [indmin,indmax] = extr(a);

    ind1=(indmax+inc1-p(i));

    ind1 = ind1( (ind1 > -nb_p1 & ind1 < nb_p1)==1 );
    ind=find(abs(ind1)==min(abs(ind1)));  % prendre l'indice le plus proche de l'AIC picking

    if isempty(ind1)==1
        indice=0;
    elseif length(ind) > 1
        indice = ind1(ind(1)); % si AIC se trouve a mi-chemin entre un chang. de phase prendre l'indice avant l'AIC picking
    elseif length(ind) == 1
        indice = ind1(ind);
    end

    pc1= p(i)+indice-1; % nouvelle position
	if pc1<1
		pc1=1;
	end
    
    if abs(y(pc1,i)) > abs( y(p(i),i)) % eviter d'aller sur un maximum
        pc1=p(i);
    end

    inc1= pc1-nb_p2;
    inc2=pc1-2;

    if inc1 <=0
        inc1=1;
    end

    a=y(inc1:inc2,i);
    a=a(:)';
    [indmin] = extr(a);

    if SNR(i) < 10 && abs(y(pc1,i)-mean(y(1:pc1-4,i))) > 2*stdn && isempty(indmin)==0;
        pc(i) = pc1-2+(indmin(end)-length(a));
    elseif abs(y(pc1,i)) < 2*stdn && (isempty(indmin)==1 || p(i)-indmin(end) > nb_p1) ;
        pc(i)=p(i);
    else
        pc(i)=pc1;
    end

    if rem(i, inc_wb)==0, waitbar(i/M, h_AIC), end
end
close(h_AIC)

%warning on