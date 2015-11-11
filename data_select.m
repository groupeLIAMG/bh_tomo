function [SNR,ind_data_select,ind_max]=data_select(data,freq,dt,L,seuil,medfilt_lth)
% [SNR,ind_data_select,ind_max]=data_select(data,freq,dt,L,seuil,medfilt_lth)
%
% input
%%%%%%%
% data :  data matrix that contain all traces
% freq : nominal frequency in Hertz
% dt : sampling period in s
% L : noise window length
% seuil : threshold
% medfilt_lth : length of median filter
% output  
%%%%%%%%%%
% SNR : radar signal to noise
% ind_data_select : selected data index
% ind_max : index of the maximum

% std_env : stadard deviation of trace enveloppe

if nargin <= 3
    disp('100''s last samples in all traces are considered as noise')
    L=100;
    seuil=5;
    medfilt_lth = 10;
end

if nargin <= 4
    seuil=5;
    medfilt_lth = 10;
end

if nargin <=5
    medfilt_lth = 10;
end

%h_wb = waitbar(0,'Trace selection ...Please wait ...');

%std_env= std(abs(hilbert(data)));

[N,M]=size(data);
std_sig=zeros(M,1);
ind_data_select= false(M,1);
ind_max= zeros(M,1);
nb_p=round(1/(dt*freq));
if medfilt_lth>0
	data = medfilt1(data, medfilt_lth);
end

%inc_wb = round(M/100);
for i=1:M
    
    [Amax, ind1]=max(data(:,i));
    ind_max(i) = ind1;
    ind = ind1-nb_p:ind1+2*nb_p;
    
    if ind(1) < 1
        ind = 1:ind1+60;
    elseif ind(end) > N
        ind = ind1-60:ind1;
    end
    
    std_sig(i,1)=std(data(ind,i));
%	if rem(i, inc_wb)==0, waitbar(i/M, h_wb), end
end
%close(h_wb)

std_noise = std(data(end-L:end,:))';

SNR = std_sig./std_noise;

ind_data_select(SNR > seuil) = true;
