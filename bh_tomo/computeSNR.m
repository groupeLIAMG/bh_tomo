function SNR = computeSNR(data)
    nptsptrc = size(data,1);
    ntrace = size(data,2);
    SNR = ones(size(1,ntrace));
    [~,i] = max(abs(data));

    largeur = 60;

    i1 = i-largeur/2;
    i2 = i+largeur/2;
    i1(i1<1) = 1;
    i2(i2>nptsptrc) = nptsptrc;
    for n=1:ntrace
        SNR(n) = std(data(i1(n):i2(n), n))/std(data(1:largeur, n));
        if isnan(SNR(n))
            SNR(n) = 1e-9;
        end
    end
end

