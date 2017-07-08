function Interference = Interference_MBS(FBS, BS, sigma2, NumRealization) % inputs are dBm, output is not in dB
    % compute power of signal at FUE received from MBS        
    f=0.9; % Frequency in GHz
    PLi = -1.8*f^2+10.6*f-5.5;
    H = abs((1/sqrt(2)) * (randn(1,NumRealization)+1i*randn(1,NumRealization)));
    hbs = H.^2;
    sigma = 10^((sigma2-30)/10);

    d = sqrt((BS.X-FBS.FUEX).^2+(BS.Y-FBS.FUEY).^2);
    PL0 = 62.3+32.*log10(d/5);
    PL_FB = PL0 + PLi;
    pbs = 10.^((BS.P-PL_FB-30)/10);

    Interference = sum(pbs.*hbs+sigma)/NumRealization;
end