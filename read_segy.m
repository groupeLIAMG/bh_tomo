function s = read_segy(segyfile, varargin)
% READ_SEGY - read the content of a SEG-Y file
%  s = read_segy(segyfile,  traces, fields)
%
% Input:
%    segyfile (mandatory) : name of SEG-Y file
%    traces (optional)  : vector of desired traces within the file, give
%                          empty array to get all traces
%    fields (optional)  : indices or names of trace header word to retreive
%                - indices are given as a vector such as [3 6 8:9]
%                - names are separated by commas such as 'sx,sy,gx,gy' (see
%                  read_segy_tr_headers for details)
%
% Output:
%    s : a struct variable holding
%                - the binary header data (s.bh)
%                - the trace header data (s.th)
%                - the trace data (s.data)
%
% Caveat:
%     Text headers are not read
%     Traces of variable length are not handled (results unpredictable)
%
% ---
%
% Bernard Giroux
% INRS-ETE
% 2010-08-06

traces = [];
fields = [];
if nargin > 1
    traces = varargin{1};
end
if nargin > 2
    fields = varargin{2};
end

s.bh = read_segy_b_header(segyfile);
if isempty( fields )
    s.th = [];
else
    s.th = read_segy_tr_headers(segyfile, traces, fields);
end
s.data = read_segy_data(segyfile, traces);
