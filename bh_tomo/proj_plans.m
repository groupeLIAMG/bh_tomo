function [p_data,no_plan] = proj_plans(data, plans)
%  p_data = proj_plan(data, x0, a)
%
% data : coordonnees a projeter (nx3)
% plans: vecteur des plans, contenant struct
%        x0   : pt sur le plan (1x3)
%        a    : cosinus directeur de la normale au plan (1x3)
%
% p_data : coordonnees projetees
% no_plan : no du plan de projection
%
%

% Copyright (C) 2005 Bernard Giroux
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

p_data = data;
no_plan = zeros(1,size(data,1));
p = zeros(1,length(plans));
for n=1:size(data,1)
    for nn=1:length(plans)
        r = plans(nn).x0-data(n,:);           % vecteur pointant de data vers x0
        p(nn) = abs(dot(plans(nn).a, r));          % distance entre data et le plan
    end
    [~,no] = min(p);
    p_data(n,:) = data(n,:) + p(no)*plans(no).a;  % coord de data projete sur le plan
    no_plan(n) = no;
end
