function [st,t,f] = st_bg(timeseries,minfreq,maxfreq,samplingrate,freqsamplingrate)
% Returns the Stockwell Transform of the timeseries.
% Code by Robert Glenn Stockwell.
% DO NOT DISTRIBUTE
% BETA TEST ONLY
% Reference is "Localization of the Complex Spectrum: The S Transform"
% from IEEE Transactions on Signal Processing, vol. 44., number 4, April 1996, pages 998-1001.
%
%-------Inputs Needed------------------------------------------------
%  
%   *****All frequencies in (cycles/(time unit))!******
%	"timeseries" - vector of data to be transformed
%-------Optional Inputs ------------------------------------------------
%
%"minfreq" is the minimum frequency to be sampled (Default=0)
%"maxfreq" is the minimum frequency to be sampled (Default=Nyquist)
%"samplingrate" is the minimum frequency to be sampled (Default=1)
%"freqsamplingrate" is the frequency-sampling interval of the spectrum (Default=1)
%Passing a negative number will give the default ex.  [s,t,f] = st(data,-1,-1,2,2)
%-------Outputs Returned------------------------------------------------
%
% st     - a complex matrix containing the Stockwell transform. 
%			 The rows of STOutput are the frequencies and the 
%            columns are the time values ie each column is 
%            the "local spectrum" for that point in time
%  t     - a vector containing the sampled times
%  f     - a vector containing the sampled frequencies
%--------Additional details-----------------------
%   %  There are several parameters immediately below that
%  the user may change. They are:
%[power]      if true gives the local power spectrum.
%[amplitude]  if true gives the local amplitude spectrum.
%                Both [power] and [amplitude] false then returns the local complex-valued spectrum.
%[removeedge] if true, removes a least squares fit parabola
%                and puts a 5% hanning taper on the edges of the time series.
%                This is usually a good idea.
%[analytic_signal]  if the timeseries is real-valued
%                      this takes the analytic signal and STs it.
%                      This is almost always a good idea.
%[factor]     the width factor of the localizing gaussian
%                ie, a sinusoid of period 10 seconds has a 
%                gaussian window of width factor*10 seconds.
%                I usually use factor=1, but sometimes factor = 3
%                to get better frequency resolution.
%   Copyright (c) by Bob Stockwell
%   $Revision: 1.2 $  $Date: 1997/07/08  $


% This is the S transform wrapper that holds default values for the function.
TRUE = 1;
FALSE = 0;
%%% DEFAULT PARAMETERS  [change these for your particular application]
power = FALSE;           
amplitude = TRUE;       
removeedge= TRUE;
analytic_signal =  TRUE;
factor = 1;
%%% END of DEFAULT PARAMETERS


%%%START OF INPUT VARIABLE CHECK
% First:  make sure it is a valid time_series 
%         If not, return the help message

if nargin == 0 
  st_help
  t=0;,st=-1;,f=0;
  return
end

% Change to column vector
if size(timeseries,2) > size(timeseries,1)
  timeseries=timeseries';	
end

% Make sure it is a 1-dimensional array
if size(timeseries,2) > 1
  error('Please enter a *vector* of data, not matrix')
  return
elseif (size(timeseries)==[1 1]) == 1
  error('Please enter a *vector* of data, not a scalar')
  return
end

% use defaults for input variables

if nargin == 1
  minfreq = 0;
  maxfreq = fix(length(timeseries)/2);
  samplingrate=1;
  freqsamplingrate=1;
elseif nargin==2
  maxfreq = fix(length(timeseries)/2);
  samplingrate=1;
  freqsamplingrate=1;
  [ minfreq,maxfreq,samplingrate,freqsamplingrate] =  check_input(minfreq,maxfreq,samplingrate,freqsamplingrate,timeseries);
elseif nargin==3 
  samplingrate=1;
  freqsamplingrate=1;
  [ minfreq,maxfreq,samplingrate,freqsamplingrate] =  check_input(minfreq,maxfreq,samplingrate,freqsamplingrate,timeseries);
elseif nargin==4   
  freqsamplingrate=1;
  [ minfreq,maxfreq,samplingrate,freqsamplingrate] =  check_input(minfreq,maxfreq,samplingrate,freqsamplingrate,timeseries);
elseif nargin == 5
  [ minfreq,maxfreq,samplingrate,freqsamplingrate] =  check_input(minfreq,maxfreq,samplingrate,freqsamplingrate,timeseries);
else      
  minfreq = 0;
  maxfreq = fix(length(timeseries)/2);
  samplingrate=1;
  freqsamplingrate=1;
end
%END OF INPUT VARIABLE CHECK

% If you want to "hardwire" minfreq & maxfreq & samplingrate & freqsamplingrate do it here

% calculate the sampled time and frequency values from the two sampling rates
t = (0:length(timeseries)-1)*samplingrate;
spe_nelements =ceil((maxfreq - minfreq+1)/freqsamplingrate)   ;
f = (minfreq + [0:spe_nelements-1]*freqsamplingrate)/(samplingrate*length(timeseries));

% The actual S Transform function is here:
st = strans(timeseries,minfreq,maxfreq,samplingrate,freqsamplingrate,removeedge,analytic_signal,factor); 
% this function is below, thus nicely encapsulated

if power
  st = abs(st).^2;
elseif amplitude
  st = abs(st);
end

return


%^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
%^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
%^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
%^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
%^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^


function st = strans(timeseries,minfreq,maxfreq,samplingrate,freqsamplingrate,removeedge,analytic_signal,factor); 
% Returns the Stockwell Transform, STOutput, of the time-series
% Code by R.G. Stockwell.
% Reference is "Localization of the Complex Spectrum: The S Transform"
% from IEEE Transactions on Signal Processing, vol. 44., number 4,
% April 1996, pages 998-1001.
%
%-------Inputs Returned------------------------------------------------
%         - are all taken care of in the wrapper function above
%
%-------Outputs Returned------------------------------------------------
%
%	ST    -a complex matrix containing the Stockwell transform.
%			 The rows of STOutput are the frequencies and the
%			 columns are the time values
%
%
%-----------------------------------------------------------------------
  
% Compute the length of the data.
n=length(timeseries);
original = timeseries;
if removeedge
  ind = [0:n-1]';
  r = polyfit(ind,timeseries,2);
  fit = polyval(r,ind) ;
  timeseries = timeseries - fit;
  sh_len = floor(length(timeseries)/10);
  wn = hanning(sh_len);
  if(sh_len==0)
	sh_len=length(timeseries);
	wn = 1&[1:sh_len];
  end
  % make sure wn is a column vector, because timeseries is
  if size(wn,2) > size(wn,1)
	wn=wn';	
  end
  
  timeseries(1:floor(sh_len/2),1) = timeseries(1:floor(sh_len/2),1).*wn(1:floor(sh_len/2),1);
  timeseries(length(timeseries)-floor(sh_len/2):n,1) = timeseries(length(timeseries)-floor(sh_len/2):n,1).*wn(sh_len-floor(sh_len/2):sh_len,1);
end

% If vector is real, do the analytic signal 

if analytic_signal
  % this version of the hilbert transform is different than hilbert.m
  %  This is correct!
  ts_spe = fft(real(timeseries));
  h = [1; 2*ones(fix((n-1)/2),1); ones(1-rem(n,2),1); zeros(fix((n-1)/2),1)];
  ts_spe(:) = ts_spe.*h(:);
  timeseries = ifft(ts_spe);
end  

% Compute FFT's
vector_fft=fft(timeseries);
vector_fft=[vector_fft,vector_fft];

% Preallocate the STOutput matrix
st=zeros(ceil((maxfreq - minfreq+1)/freqsamplingrate),n);
% Compute the mean
% Compute S-transform value for 1 ... ceil(n/2+1)-1 frequency points
if minfreq == 0
  st(1,:) = mean(timeseries)*(1&[1:1:n]);
else
  st(1,:)=ifft(vector_fft(minfreq+1:minfreq+n).*g_window(n,minfreq,factor));
end

%the actual calculation of the ST
% Start loop to increment the frequency point
for banana=freqsamplingrate:freqsamplingrate:(maxfreq-minfreq)
  st(banana/freqsamplingrate+1,:)=ifft(vector_fft(minfreq+banana+1:minfreq+banana+n).*g_window(n,minfreq+banana,factor));
end
% a fruit loop!   aaaaa ha ha ha ha ha ha ha ha ha ha
% End loop to increment the frequency point

%%% end strans function

%------------------------------------------------------------------------
function gauss=g_window(length,freq,factor)

% Function to compute the Gaussion window for 
% function Stransform. g_window is used by function
% Stransform. Programmed by Eric Tittley
%
%-----Inputs Needed--------------------------
%
%	length-the length of the Gaussian window
%
%	freq-the frequency at which to evaluate
%		  the window.
%	factor- the window-width factor
%
%-----Outputs Returned--------------------------
%
%	gauss-The Gaussian window
%

vector(1,:)=[0:length-1];
vector(2,:)=[-length:-1];
vector=vector.^2;    
vector=vector*(-factor*2*pi^2/freq^2);
% Compute the Gaussion window
gauss=sum(exp(vector));

%-----------------------------------------------------------------------

%^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^%
function [ minfreq,maxfreq,samplingrate,freqsamplingrate] =  check_input(minfreq,maxfreq,samplingrate,freqsamplingrate,timeseries)
% this checks numbers, and replaces them with defaults if invalid

% if the parameters are passed as an array, put them into the appropriate variables
s = size(minfreq);
l = max(s);
if l > 1  
  temp=minfreq;
  minfreq = temp(1);;
  if l > 1  maxfreq = temp(2);,end;
  if l > 2  samplingrate = temp(3);,end;
  if l > 3  freqsamplingrate = temp(4);,end;
end      

if minfreq < 0 | minfreq > fix(length(timeseries)/2);
  minfreq = 0;
end
if maxfreq > length(timeseries)/2  | maxfreq < 0 
  maxfreq = fix(length(timeseries)/2);
end
if minfreq > maxfreq 
  temporary = minfreq;
  minfreq = maxfreq;
  maxfreq = temporary;
  clear temporary;
end
if samplingrate <0
  samplingrate = abs(samplingrate);
end
if freqsamplingrate < 0   % check 'what if freqsamplingrate > maxfreq - minfreq' case
  freqsamplingrate = abs(freqsamplingrate);
end

% bloody odd how you don't end a function

%^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^%
function st_help
disp(' ')
disp('st()  HELP COMMAND')
disp('st() returns  - 1 or an error message if it fails')
disp('USAGE::    [localspectra,timevector,freqvector] = st(timeseries)')
disp('NOTE::   The function st() sets default parameters then calls the function strans()')
disp(' ')  
disp('You can call strans() directly and pass the following parameters')
disp(' **** Warning!  These inputs are not checked if strans() is called directly!! ****')
disp('USAGE::  localspectra = strans(timeseries,minfreq,maxfreq,samplingrate,freqsamplingrate,removeedge,analytic_signal,factor) ')

disp(' ')
disp('Default parameters (available in st.m)')
disp('POWER            - returns the power spectrum')
disp('ABS              - returns the absolute value spectrum')
disp('REMOVEEDGE       - removes the edge with a 5% taper, and takes')
disp('FACTOR           -  the width factor of the localizing gaussian')
disp('                    ie, a sinusoid of period 10 seconds has a ')
disp('                    gaussian window of width factor*10 seconds.')
disp('                    I usually use factor=1, but sometimes factor = 3')
disp('                    to get better frequency resolution.')
disp(' ')
disp('Default input variables')
disp('MINFREQ           - the lowest frequency to be sampled (Default=0)')
disp('MAXFREQ           - the highest frequency to be sampled (Default=nyquist')
disp('SAMPLINGRATE      - the time interval between successive data points (Default = 1)')
disp('FREQSAMPLINGRATE  - the number of frequencies between samples')

% end of st_help procedure   
  

%[verbose]    if true %[power]      if true gives the local power spectrum.
%[amplitude]  if true gives the local amplitude spectrum.
%                Both [power] and [amplitude] false then returns the local complex-valued spectrum.
%[removeedge] if true, removes a least squares fit parabola
%                and puts a 5% hanning taper on the edges of the time series.
%                This is usually a good idea.
%[analytic_signal]  if the timeseries is real-valued
%                      this takes the analytic signal and STs it.
%                      This is almost always a good idea.
%[factor]     %             %              %                I usually use factor=1, but sometimes factor = 3
