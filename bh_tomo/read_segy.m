function s = read_segy(segyfile, varargin)
% READ_SEGY - read the content of a SEG-Y file
%  s = read_segy(segyfile, traces, fields, dict, word_length)
%
% Input:
%    segyfile (mandatory) : name of SEG-Y file
%    traces (optional)  : vector of desired traces within the file, give
%                          empty array to get all traces
%    fields (optional)  : indices or names of trace header word to retreive
%                - indices are given as a vector such as [3 6 8:9]
%                - names are separated by commas such as 'sx,sy,gx,gy' (see
%                  read_segy_tr_headers for details)
%    dict (optional)  : custom dictionary for trace header (cell array
%                       of strings)
%    word_length (optional)  : word length in bytes for trace header
%                              (vector of double)
%                              (mandatory if dict given)
%	               - sum(word_length) must be less than 241
%                - numel(word_length) must equal numel(dict)
%                - 4-byte IBM float words are handled by setting word_length=5
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
dict = [];
word_length = [];

if nargin > 1
    traces = varargin{1};
end
if nargin > 2
    fields = varargin{2};
end
if nargin > 3
    dict = varargin{3};
end
if nargin > 4
    word_length = varargin{4};
end

s.bh = read_segy_b_header(segyfile);
if isempty( fields )
    s.th = [];
elseif isempty( dict );
    s.th = read_segy_tr_headers(segyfile, traces, fields);
else
    s.th = read_segy_tr_headers(segyfile, traces, fields, dict, word_length);
end
s.data = read_segy_data(segyfile, traces);
