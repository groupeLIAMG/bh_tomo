function [i0,i1] = findCommonTraces(no0, no1)
%
%
%

i0 = nan(numel(no0),1);
i1 = nan(numel(no1),1);

for n0=1:numel(no0)
	for n1=1:numel(no1)
		if no0(n0)==no1(n1)
			i0(n0) = n0;
			i1(n1) = n1;
		end
	end
end

i0=i0(~isnan(i0));
i1=i1(~isnan(i1));
