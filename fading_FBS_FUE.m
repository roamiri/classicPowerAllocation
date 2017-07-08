function fading = fading_FBS_FUE(FBS, NumRealization) % inputs are dBm, output is not in dB
    % compute fading of signal at FUE received from FBS        
    H = abs((1/sqrt(2)) * (randn(1,NumRealization)+1i*randn(1,NumRealization)));
    h = H.^2;

    d = sqrt((FBS.X-FBS.FUEX).^2+(FBS.Y-FBS.FUEY).^2);
    PL0 = 62.3+40*log10(d/5);
    l = 10.^((PL0-30)/10);

    fading = sum(h./l)/NumRealization;
end