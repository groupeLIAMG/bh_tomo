function [i0, i1] = findCommonTxRx(TxRx0, TxRx1)
%
%  
%

i0 = nan(size(TxRx0,1),1);
i1 = nan(size(TxRx1,1),1);

for n0=1:size(TxRx0,1)
	for n1=1:size(TxRx1,1)
		if TxRx0(n0,:)==TxRx1(n1,:)
			i0(n0) = n0;
			i1(n1) = n1;
		end
	end
end

i0=i0(~isnan(i0));
i1=i1(~isnan(i1));
