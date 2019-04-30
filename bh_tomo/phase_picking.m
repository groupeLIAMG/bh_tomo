function pick=phase_picking(mog, par)


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
%               *        giroux@geo.polymtl.ca      *                     %
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
%                                                                         %
%  input                                                                  %
%  =====                                                                  %
% mog        : multi-offset-gather data structure                         %
%                                                                         %
% TH         : threshold for data selection                               %
%                                                                         %
% output                                                                  %
% ======                                                                  %
%                                                                         %
% pick.pAIC           : AIC picking                                       %
% pick.pr             : corrected AIC picking or final picking            %
% pick.stdp           : standard deviation of picking                     %
% pick.std_pAIC_pr    : difference beteween corrected picking and AIC     %
%                       picking                                           %
% pick.SNRR=SNR       : signal to noise ratio of raw radar data           %
% pick.SNRD           : signal to noise ratio after wavelet denoising and %
%                       median filtering                                  %
% pick.ind_select     : index of selected data 1: for selected,           %
%                                              0 : for not selected       %
% pick.datad          : data after wavelet denoising and median filtering %
% pick.data           : data after detrending without denoising           %
%                                                                         %
% Used functions / fonctions utilisées :                                  %
% ======================================                                  %
% detrend_poly, data_select, calcul_AIC, AIC_correction, calcul_stdp      %
%                                                                         %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%=========================================================================%

if nargin < 2
    par.trace_cut = 300;
    par.L = 100;
    par.medfilt_lth = 10;
    par.TH = 5; % default value for data selection
    par.denoise = true;
	par.AIC_corr = true;
    par.inv_polarite_cwt = false;
end

fac_f = 1.0;
fac_t = 1.0;
if strcmp(mog.data.tunits,'ns')
    % radar: freq nominale donnée en MHz
    fac_f = 10^6;
    fac_t = 10^-9;
elseif strcmp(mog.data.tunits,'ms')
    % sismique: on assume que f dominante est en kHz
    disp('warning : assuming Tx nominal frequency in _kHz_')
    fac_f = 10^3;
    fac_t = 10^-3;
else
    disp('warning : assuming Tx nominal frequency in _Hz_')
    disp('warning : assuming time step in _s_')
end

freq = mog.data.rnomfreq * fac_f;
dt = mog.data.timec * fac_t;



% polynomial data detrend with polynomial of 4th order
%======================================================
%[M,N]=size(mog.data.rdata);
%[data]=detrend_poly(1:M,mog.data.rdata,1);

data = mog.traces(:,mog.in);

no_traces = 1:mog.data.ntrace;
no_traces = no_traces(mog.in);

% data selection
%===============
[SNRR,ind_select,ind_max] = data_select(data,freq,dt,par.L,par.TH,par.medfilt_lth);

SNR=SNRR(ind_select);
ind_max = ind_max(ind_select);
data = data(:,ind_select);
no_traces = no_traces(ind_select);
%-------------------------

[M,N]=size(data); % data size after selection


%wavelet filtering for background noise removal
%===============================================
if par.denoise
	[data2,ind_max] = filtrage_wavelet2(data);
else
	data2 = data;
end
% radar wavelet detection
m=max(ind_max);
m = m+par.trace_cut;
if m > M
    m=M;
end

C=zeros(m,N);
data=data(1:m,:);
data2=data2(1:m,:);


h_wav = waitbar(0,'Wavelet transform calculation ...Please wait ...');
L=length(data(:,1));
frq=0:1/(dt*(L-1)):1/(dt);

% fenêtre de sélection de la fréquence dominante (en Hz)

F1 = 0.2*freq;
F2 = 2.0*freq;

indf=find(frq > F1 & frq < F2);
frq=frq(indf);

signe = 1;
if par.inv_polarite_cwt
    signe = -1;
end

Fa = zeros(N,1);
inc_wb = round(N/100);
for i=1:N

    Fa1 = get_fdom(data2(:,i), dt, freq,F1,F2);
    if Fa1 == -1
        ff=abs(fft(data2(:,i)));
        indFa=find(ff(indf)==max(ff(indf)));
        Fa1 = frq(indFa(1));
    end
	
    a=1/(dt*Fa1);
    c1=cwt(signe*data2(:,i),a,'cmor1.2-1');%'fbsp2-1-.5'
    C(:,i)=c1';
    Fa(i)=Fa1;
    if rem(i, inc_wb)==0, waitbar(i/N, h_wav), end
end

close(h_wav)

%[p, p_loc] = calcul_AIC(data2,ind_max,Fa,dt,SNR); % AIC calculation for denoised data
p = calcul_AIC(data2,ind_max,Fa,dt,SNR); % AIC calculation for denoised data

[pr] = AIC_correction(C,data2,ind_max,p,Fa,dt,SNR);

%[stdp, SNRD] = calcul_stdp(data2,C,pr,Fa,dt,SNR);
stdp = calcul_stdp(data2,C,pr,Fa,dt,SNR);

pick.tt = zeros(1,mog.data.ntrace)-1;
pick.et = zeros(1,mog.data.ntrace)-1;
pick.tt_fait = false(1,mog.data.ntrace);
pick.pAIC = zeros(1,mog.data.ntrace)-1;
pick.SNRR = zeros(1,mog.data.ntrace)-1;
pick.Fa = zeros(1,mog.data.ntrace)-1;

pick.tt(no_traces) = pr * dt/fac_t;
pick.et(no_traces) = stdp * dt/fac_t;
pick.tt_fait(no_traces) = true;
pick.SNRR(no_traces) = SNR;      % signal to noise ratio of raw radar data

pick.pAIC(no_traces)        = p * dt/fac_t;          % AIC picking
pick.no_traces = no_traces;
% pick.pAIC_loc    = p_loc;      % number of samples between AIC onset picking and nearest local minimum.
% pick.pr          = pr;         % corrected AIC picking or final picking
% pick.stdp        = stdp;       % standard deviation of picking
% pick.std_pAIC_pr = pr-p;       % difference beteween corrected picking and AIC picking
% pick.SNRD        = SNRD;       % signal to noise ratio after wavelet denoising and median filtering
% pick.ind_select  = ind_select; % index of selected data 1: for selected, 0 : for not selected
pick.datad       = data2;      % data after wavelet denoising and median filtering
pick.data        = data;       % data after detrending without denoising
% pick.std_env     = std_env;    % trace enveloppe standard deviation
pick.Fa(no_traces) = Fa;

