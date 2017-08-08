function [fading, Loss] = fading_FBS_MUE(FBS, mue, NumRealization) % output is not in dB
    % compute power of signal at UE received from BS        
    xx = mue.X;
    yy = mue.Y;

    H = abs((1/sqrt(2)) * (randn(1,NumRealization)+1i*randn(1,NumRealization)));
    alpha = H.^2; %small fading

    % compute power of interference at MUE received from Femtocell
    f=0.9; % Frequency in GHz
    PLi = -1.8*f^2+10.6*f-5.5;
    d = sqrt((FBS.X-xx).^2+(FBS.Y-yy).^2);
    PL0 = 62.3+32.*log10(d/5);
    PL_FB = PL0 + PLi;
%     Loss = PL_FB;
    Loss = 10.^((PL_FB)/10); %large scale fading

    fading = sum(alpha)/NumRealization;
end