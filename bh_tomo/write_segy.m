function write_segy(segyfile, s, varargin)
%WRITE_SEGY - save the content
%   Detailed explanation goes here


s.bh.format = int16(5);  % data sample format code set to 5 (4-byte IEEE floating point)
s.bh.extfh  = int16(0);  % no extended textual file headers

% write textual file header (empty for now)
fid = fopen(segyfile, 'w');
if fid == -1
    error(['Cannot open ',segyfile,' for writing'])
end
fprintf(fid,'C 1                                                                            \n');
fprintf(fid,'C 2                                                                            \n');
fprintf(fid,'C 3                                                                            \n');
fprintf(fid,'C 4                                                                            \n');
fprintf(fid,'C 5                                                                            \n');
fprintf(fid,'C 6                                                                            \n');
fprintf(fid,'C 7                                                                            \n');
fprintf(fid,'C 8                                                                            \n');
fprintf(fid,'C 9                                                                            \n');
fprintf(fid,'C10                                                                            \n');
fprintf(fid,'C11                                                                            \n');
fprintf(fid,'C12                                                                            \n');
fprintf(fid,'C13                                                                            \n');
fprintf(fid,'C14                                                                            \n');
fprintf(fid,'C15                                                                            \n');
fprintf(fid,'C16                                                                            \n');
fprintf(fid,'C17                                                                            \n');
fprintf(fid,'C18                                                                            \n');
fprintf(fid,'C19                                                                            \n');
fprintf(fid,'C20                                                                            \n');
fprintf(fid,'C21                                                                            \n');
fprintf(fid,'C22                                                                            \n');
fprintf(fid,'C23                                                                            \n');
fprintf(fid,'C24                                                                            \n');
fprintf(fid,'C25                                                                            \n');
fprintf(fid,'C26                                                                            \n');
fprintf(fid,'C27                                                                            \n');
fprintf(fid,'C28                                                                            \n');
fprintf(fid,'C29                                                                            \n');
fprintf(fid,'C30                                                                            \n');
fprintf(fid,'C31                                                                            \n');
fprintf(fid,'C32                                                                            \n');
fprintf(fid,'C33                                                                            \n');
fprintf(fid,'C34                                                                            \n');
fprintf(fid,'C35                                                                            \n');
fprintf(fid,'C36                                                                            \n');
fprintf(fid,'C37                                                                            \n');
fprintf(fid,'C38                                                                            \n');
fprintf(fid,'C39                                                                            \n');
fprintf(fid,'C40                                                                            \n');
fclose(fid);

write_segy_b_header(segyfile,s.bh)
write_segy_traces(segyfile,s.th,s.data,varargin{:})
end

