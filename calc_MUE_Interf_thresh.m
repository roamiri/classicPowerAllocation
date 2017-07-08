function I_th = calc_MUE_Interf_thresh(MBS, mue, R_MUE, sigma2, NumRealization) % inputs are dBm, output is not in dB
    % compute Interference threshold which is derived from MUE required
    % rate
    H = abs((1/sqrt(2)) * (randn(1,NumRealization)+1i*randn(1,NumRealization)));
    h = H.^2;

    d = sqrt((MBS.X-mue.X).^2+(MBS.Y-mue.Y).^2);
    PL0 = 62.3+40*log10(d/5);
    p = 10.^((MBS.P-PL0-30)/10);
    sigma = 10.^((sigma2-30)/10);
    I_th = sum((p.*h)./(2^R_MUE-1))/NumRealization - sigma;
end