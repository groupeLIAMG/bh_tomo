% READ_SEGY_DATA - read the trace data of a SEG-Y file
%  d = read_segy_data(segyfile,  traces)
%
% Input:
%    segyfile (mandatory) : name of SEG-Y file
%    traces (optional)  : vector of desired traces within the file, give
%                         empty array to get all traces
%
% Output:
%    d : matrix of size n_samples x n_traces
%
%
% ---
%
% Bernard Giroux
% INRS-ETE
% 2010-08-04
%