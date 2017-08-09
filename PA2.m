%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Power Allocation in a femtocell network based on dual decomposition
% This file is written based on solution in page 2 of June 19-25, 2017
% notes
% In this optimization problem Pmin is not considered in the constraints.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function Q = PA2(fbsCount, NumRealization,MaxIteration, showPlot)
format long;
%% Initialization
%clear;
% clc;
% format short
% format compact

%% Parameters
Pmin = -20; %dBm
pmin = 10^((Pmin-30)/10);

Pmax = 25; %dBm
pmax = 10^((Pmax-30)/10);
ppmax = pmax/(2e7);

PBS = 50 ; %dBm
% Ith = 1; % MUE maximum interference threshold which represents its desired QoS
% R = 1; % FUE minimum required rate which represents QoS
% K = 16;  % Number of femtocells

% lambda = cell(1,K); % Lagrange variable for Pmax
% gamma = cell(1,K);  % Lagrange variable for Pmin
% Mu = cell(1,K);    % Lagrange variable for Rate
% eta = [];   % Lagrange variable for Ith
% MaxIteration = 200000;
R_MUE = 8.0;
R_FUE =4.0; %FUE Rate requirements
%% Initialize MBS and MUE
MBS = BaseStation(0 , 0 , PBS); % (0,0) is the location and P_MBS = 50 dBm
mue = UE(-200, 0);
CC = calc_MUE_Cap_NoInterf(MBS, mue, -120, NumRealization);
%% Initialize FBSs
%Generate fbsCount=16 FBSs
FBS_Max = cell(1,16);
for i=1:3
%     if i<= fbsCount
        FBS_Max{i} = FemtoStation(180+(i-1)*35,150, MBS, mue, 10);
%     end
end

for i=1:3
%     if i+3<= fbsCount
        FBS_Max{i+3} = FemtoStation(165+(i-1)*30,180, MBS, mue, 10);
%     end
end

for i=1:4
%     if i+6<= fbsCount
        FBS_Max{i+6} = FemtoStation(150+(i-1)*35,200, MBS, mue, 10);
%     end
end

for i=1:3
%     if i+10<= fbsCount
        FBS_Max{i+10} = FemtoStation(160+(i-1)*35,240, MBS, mue, 10);
%     end
end

for i=1:3
%     if i+13<= fbsCount
        FBS_Max{i+13} = FemtoStation(150+(i-1)*35,280, MBS, mue, 10);
%     end
end
%%
% 
% fbsCount = K; % notice the K and fbsCount varibales when making a function
FBS = cell(1,fbsCount);

if fbsCount>=1, FBS{1} = FBS_Max{1}; end
if fbsCount>=2, FBS{2} = FBS_Max{3}; end
if fbsCount>=3, FBS{3} = FBS_Max{14}; end
if fbsCount>=4, FBS{4} = FBS_Max{16}; end
if fbsCount>=5, FBS{5} = FBS_Max{9}; end
if fbsCount>=6, FBS{6} = FBS_Max{4}; end
if fbsCount>=7, FBS{7} = FBS_Max{2}; end
if fbsCount>=8, FBS{8} = FBS_Max{15}; end
if fbsCount>=9, FBS{9} = FBS_Max{10}; end
if fbsCount>=10, FBS{10} = FBS_Max{12}; end
if fbsCount>=11, FBS{11} = FBS_Max{5}; end
if fbsCount>=12, FBS{12} = FBS_Max{7}; end
if fbsCount>=13, FBS{13} = FBS_Max{11}; end
if fbsCount>=14, FBS{14} = FBS_Max{6}; end
if fbsCount>=15, FBS{15} = FBS_Max{8}; end
if fbsCount>=16, FBS{16} = FBS_Max{13}; end

%% Main Loop
textprogressbar(sprintf('calculating outputs:'));
I_th = 1e6*calc_MUE_Interf_thresh(MBS, mue, R_MUE, -120, NumRealization);
interf = [];
for i=1:size(FBS,2)
    fbs = FBS{i};
    fbs = fbs.setR_FUE(R_FUE);
    FBS{i} = fbs;
end

for i=1:size(FBS,2)
    fbs = FBS{i};
    [gMF_k, lMF_k] = fading_FBS_MUE(fbs, mue, 100);
    fbs = fbs.setGMF(gMF_k, lMF_k);
    [gF_k, lF_k] = fading_FBS_FUE(fbs,100);
    fbs = fbs.setGF(gF_k, lF_k);
    I_k = Interference_MBS(fbs, MBS, -120, NumRealization);
    fbs = fbs.setInterf(I_k);
    FBS{i} = fbs;
end

interf = [];
for iter=1:MaxIteration
    textprogressbar((iter/MaxIteration)*100);
    eta = MBS.getLagrangeVar();
    Total_MUE_Interf = 0.0;
    for i=1:size(FBS,2)
        % Power Allocation
        fbs = FBS{i};
        [lambda ,gamma ,mu] = fbs.getLagrangeVars();
        gMF_k = fbs.gmf;% fading_FBS_MUE(fbs, mue, 100);
%         fbs = fbs.setGMF(gMF_k);
        gF_k = fbs.gf;%fading_FBS_FUE(fbs,100);
        I_k = fbs.If;%Interference_MBS(fbs, MBS, -120, NumRealization);
        p = (1+mu/log(2))*(lambda+eta*gMF_k)^(-1)-(I_k/gF_k);
        fbs = fbs.setPower(max(p,0.0));
        fbs = fbs.setCapacity(log2(1+fbs.P*gF_k/I_k));
        Total_MUE_Interf = Total_MUE_Interf + fbs.P * gMF_k;
        FBS{i} = fbs;
    end
    interf = [interf Total_MUE_Interf];
    % Update Lagrange Vars
    beta=0.01;
    for i=1:size(FBS,2)
        fbs = FBS{i};
        [lambda ,gamma ,mu] = fbs.getLagrangeVars();
%          beta = gss(FBS,i, pmax, pmin, eta, I_th, 1);
        lambda = max(lambda + (0.3) * (fbs.P - ppmax), 0.0);
%          lambda = max(beta, 0);
%          beta = gss(FBS,i, pmax, pmin, eta, I_th, 3);
        mu = max(mu+(0.005)*(fbs.R_FUE-fbs.C_FUE), 0.0);
        fbs = fbs.updateLagrangeVars(lambda, gamma, mu);
        FBS{i} = fbs;
    end
%     beta = gss(FBS,i, pmax, pmin, eta, I_th, 4);
    eta = max(eta -0.3 * (I_th - Total_MUE_Interf), 0.0);
    MBS = MBS.updateLagrangeVar(eta);
    mue_C = calc_MUE_Capacity(MBS, mue, -120, Total_MUE_Interf,NumRealization);
    mue = mue.setCapacity(mue_C);
end

fprintf('total interf = %4.4f\n', 1e4*Total_MUE_Interf);
fprintf('   Threshold = %4.4f\n', 1e4*I_th);
cc_mue = mue.C_profile;
ss = size(cc_mue,2);
cc_mean = sum(cc_mue(0.8*ss:ss))/(0.2*ss);
Q.CMUE = cc_mean;

sum_C = 0;
min_C = inf;
for i=1:size(FBS,2)
    cc_fue = FBS{i}.C_profile;
    ss = size(cc_fue,2);
    cc_mean = sum(cc_fue(0.8*ss:ss))/(0.2*ss);
    if min_C > cc_mean
        min_C = cc_mean;
    end
    sum_C = sum_C + cc_mean;
end
Q.min = min_C;
Q.sum = sum_C;

if showPlot==1
    figure;
    hold on;
    for i=1:size(FBS,2)
        plot(FBS{i}.powerProfile);
        title('Power profile');
    end
    figure;
    plot(MBS.eta);
    title('MBS Eta');

    figure;
    hold on;
    plot(FBS{1}.lambda);
    plot(FBS{2}.lambda);
    title('Lambda');

    figure;
    hold on;
    plot(FBS{1}.mu);
    plot(FBS{2}.mu);
    title('Mu');

    figure;
    plot(mue.C_profile);
    title('MUE Capacity');

    figure;
    plot(interf);
    title('Interference at MUE.');
end
end