function idebut = ltasta(trace,nlta,nsta,threshold)
idebut = 1;
smallnum = 1e-4;
nh = 20;
han = hanning(nh);

npts = numel(trace);
ht = hilbert(trace);
E = abs(ht); % enveloppe

ratios = zeros(size(trace));
lta_tot = sum(E(1:nlta));
sta_tot = sum(E(nlta+(1:nsta)));
lta = lta_tot/nlta;
sta = sta_tot/nsta;
if lta>smallnum && sta>smallnum
    ratios(nlta) = sta/lta;
end
for n=(nlta+1):(npts-nsta)
    lta_tot = lta_tot + E(n) - E(n-nlta);
    sta_tot = sta_tot + E(n+nsta) - E(n);
    lta = lta_tot/nlta;
    sta = sta_tot/nsta;
    if lta>smallnum && sta>smallnum
        ratios(n+nsta) = sta/lta;
    end
end

smooth = conv(ratios,han);
smooth = smooth(nh/2:end-nh/2);
smooth = smooth*max(ratios)/max(smooth);

pick = smooth>threshold;
nos=1:npts;
nos = nos(pick);
if ~isempty(nos)
    idebut = nos(1);
end

