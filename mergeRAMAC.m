function mergeRAMAC(outFile, prefix, nos)
% mergeRAMAC - fusionne plusieurs fichiers de données RAMAC
%
% function mergeRAMAC(outFile, prefix, nos)
%

% Copyright (C) 2007 Bernard Giroux
%
% This program is free software; you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation; either version 3 of the License, or
% (at your option) any later version.
%
% This program is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details.
% 
% You should have received a copy of the GNU General Public License
% along with this program.  If not, see <http://www.gnu.org/licenses/>.
% 
%

rdata = [];
for n=1:length(nos)
    disp(['Je traite ',prefix,num2str(nos(n))]);
    data = lisRAMAC2([prefix,num2str(nos(n))]);
    if ( data.ntrace ~= size(data.rdata, 2) )
        disp(['Erreur, nombre de trace non concordant - ',...
            [prefix,num2str(nos(n))]])
        return
    end
    if ( data.nptsptrc ~= size(data.rdata, 1) )
        disp(['Erreur, nombre de pts/trace non concordant - ',...
            [prefix,num2str(nos(n))]])
        return
    end
    rdata = [rdata data.rdata];
end

Tx_z = [];
Rx_z = [];
for n=1:length(nos)
  pos = lisTLF([prefix,num2str(nos(n))]);
  Tx_z = [Tx_z pos.Tx_z];
  Rx_z = [Rx_z pos.Rx_z];
end

[Tx_z, ind] = sort(Tx_z);
Rx_z = Rx_z(ind);
rdata = rdata(:,ind);


fid=fopen([outFile,'.rd3'],'w','ieee-le');
fwrite(fid,rdata,'int16');
fclose(fid);

fid2 = fopen([outFile,'.tlf'],'wt');
fprintf(fid2,['#First trace     Last trace     First pos     Last pos' ...
                          ' Fixed pos\n']);
fprintf(fid2, '     %d      %d      %5.2f     %5.2f     %5.2f\n', ...
                [(1:length(Tx_z))' (1:length(Tx_z))' Rx_z' Rx_z' Tx_z']');
fclose(fid2);


fid = fopen([prefix,num2str(nos(1)),'.rad'],'rt');
fid2 = fopen([outFile,'.rad'],'wt');
while 1
  tline = fgets(fid);
  if ~ischar(tline), break, end
  if ~isempty(findstr(tline,'LAST TRACE'))
        fprintf(fid2, 'LAST TRACE:%d\n', size(rdata,2));
  else
        fprintf(fid2, '%s', tline);
  end
end
fclose(fid2);
fclose(fid);


