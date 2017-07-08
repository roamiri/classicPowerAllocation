% Simulation of  reference 2 
% "Cross Layer Scheduling for OFDMA Amplify and Forward Relay Networks"
% Author  Ruhollah Amiri    amiriruhollah@gmail.com
%% 
clear all
close all
clc
addpath ./winner
MaxIter = 2;
%% Defining the Parameters
% Constants
K = 15;  % Number of Users in the cell
M = 3;   % Number of Relays
A = 3;   % Number of sectors    (Noe that since A = M so in each sector there is one Relay)
nF = 8; % Number of Subcarriers

%temporary matrices
temp1 = zeros(nF, nF, K/A);
temp2 = zeros(nF, nF, K/A);
x1 = zeros(nF, nF, K/A);
temp4 = zeros(nF, nF, K/A);
temp5 = zeros(nF, nF, K/A);

%Channel Matrices
HSR = zeros(M , nF) ;   % channel matrice between BS and Relays  ( HSR~ CN(sqrt(k/(1+k)), 1/(1+k)))  k=Rician factor
HRD = zeros(K/A , nF, M) ;  % channel matrice between Relay one and the MS that are supported by this relay
                                      % HRiD ~ CN(0, sigmaE)  sigmaE = error variance of imperfect CSI estimation at Relays from Ms to Relay channel
sigmaE = 0.01 ; % error variance of imperfect CSI estimation at Relays from Ms to Relay channel

% Transmit Powers
% PTotal = 10 ^ 1.6 ; % Total Power Constraint
PSR = zeros(K/A, nF, M) ; % transmit power of BS to Relays in first time slot
PRD = zeros(K/A, nF, M) ; % transmit power of Relays to MS in second time slot

% Large Scale Loss 
scm = scmparset;
scm.ScmOptions = 'los' ;
l = linkparset(M, 500);
PathLoss= 0.1 * pathloss(scm,l) ; % loss between BS and Relay
LSR = 10.^ (-1 .* PathLoss') ;
scm.ScmOptions = 'none' ;
l = linkparset(K/A, 1000);
LRD = ones(K/A, 1, M) ;
for relay = 1:M
    PathLoss= 0.1 * pathloss(scm,l) ; 
    LRD(:, 1, relay) =10.^ (-1 .* PathLoss') ;  % loss between Relay and MS
end

% Allocation parameters
Wk = ones(1, K/A, M) ;  % weights given by MAC layer based on Prportional fairness or max-min fairness or ... .
Wk(1, 1, :) = 2;      Wk(1, 2, :) = 1.5;     Wk(1, 3, :) = 1.1;         Wk(1, 4:5, :) = 0.2;


Epsilonk = 0.01*ones(1, K/A, M) ;  % Outage probability for each user

Carrier = cell(1, M) ; % Carrier assignment each element of cell is a (nF, nF, K/A) array.
Rate = cell(1, M) ;  % Rate allocation on each subchannel each element of cell is a (nF, nF, K/A) array.

Finv = zeros(K/A, nF, M) ;   % Inverse function of chi-square with 2 degree of freedom and noncentrally parameter of abs(HRD)^2/sigmaE
                                          % ncx2inv
Gamma = cell(1, M) ;  % equivalent SNR in calculating the rates from equation (13)  % each element is (nF, nF, K/A) array

%% Channel large scale loss calculated from SCM channel model
%  Version 1.2, Jan 11, 2005

%% Resource allocation based on proposed distributed scheduling algorithm given in figure3.

% Relays obtain CSI of BS-to-Relay links and Relay-to-user links. 
% Bs initializes all lagrange multipliers.
k = 10 ^ 0.6;
HSR =  sqrt(1/(1+k)/2)*(randn(M,nF)+1j*randn(M, nF) ) + (sqrt(k/(1+k)) + 1j* sqrt(k/(1+k)) )  * ones(M, nF) ;  % Notes: the channel gain should be derived from frequency selective channel in next steps
HRD = sqrt(1/2) * (randn(K/A, nF, M)+1j*randn(K/A, nF, M) ) ; 
HRD = HRD + sqrt(sigmaE/2) * (randn(K/A, nF, M)+1j*randn(K/A, nF, M) ) ; % Imperfect CSI in Relay
landa = zeros(1, MaxIter);
beta = zeros(nF, MaxIter);
gamma = zeros(nF, MaxIter);
landa(1 ,1) = 2.2; beta(:, 1) = rand(nF,1); gamma(:, 1) = rand(nF,1);

for iter = 1 : MaxIter
 %Each relay solves the subproblem in (17) based on its local CSI

%     matlabpool(2) ;
%     Phi = cell(1, M) ;
%     Omega = cell(1, M) ;
    for m = 1 : M
        for k = 1 : K/A
            display(k)
            for i = 1 : nF
                for j = 1 : nF
                    e = Epsilonk(1, k, m)  ;
                    e = ncx2inv(e, 2, abs(HRD(k, j, m))^2/sigmaE) ;
                    a = abs(HSR(m, i))^2 ;  % a = LSR(1, m) * abs(HSR(m, i))^2 ;
                    b =  e ;  % b = LRD(k, 1, m) * e ;
                    temp1(i, j, k) = ((sqrt(a) + sqrt(b))^2) / (a * b) ;
                    temp2(i, j, k) = sqrt(a/b) ;
                    PSR(k, i, m) = max(0 , Wk(1, k, m) * (1-Epsilonk(1, k, m))/landa(1, iter)  - temp1(i, j, k) ) / (1 + temp2(i, j, k) ) ; 
                    PRD(k, j, m) = max(0 , Wk(1, k, m) * (1-Epsilonk(1, k, m))/landa(1, iter)  - temp1(i, j, k) ) / (1 + 1/temp2(i, j, k) ) ;
                    x1 = (PSR(k , i, m) * a * PRD(k, j, m) * b) / (0.01 + PSR(k , i, m) * a + PRD(k, j, m) * b) ;
                    temp4(i, j, k) = 0.5 * log2(1 + x1) ;
                    compare = log2(1 + x1) -  x1/( 1 + x1) ;
                   
%                     if ( PSR(k, i, m) ~= 0 ) && ( PRD(k, j, m) ~= 0)
%                         display(compare);
%                         display( 2*( gamma(j, iter) + beta(i, iter) ) / ( Wk(1, k, m) * (1-Epsilonk(1, k, m)) ));
%                     end
                   
                    if compare >= 2*( gamma(j, iter) + beta(i, iter) ) / ( Wk(1, k, m) * (1-Epsilonk(1, k, m)) ) 
                        temp5(i, j, k) = 1;
                        display( 'Carrier =1')
                    else
                        temp5(i, j, k) = 0;
                    end
                end
            end
        end
%         Phi{1, m} = temp1 ;
%         Omega{1, m} = temp2 ;
%         Gamma{1, m} = temp3 ;
        Rate{1, m} = temp4;
        Carrier{1, m} = temp5 ;
    end
    
 % The BS updates the Lagrange multipliers using the demand and supply
 % adjustment rule in (25)
 counter = zeros(1, nF);

 for j = 1: nF
      for m=1 : M
         temp1 = Carrier{1, m} ;
         counter(1, j) = sum(sum(temp1(:, j, :)) );
     end
     gamma(j, iter+1) =  max ( 0, gamma(j, iter) - 2 *( 1- counter(1, j)) )  ;
 end
 
for i = 1: nF
      for m=1 : M
         temp1 = Carrier{1, m} ;
         counter(1, i) = sum(sum(temp1(i, :, :)) );
     end
     beta(i, iter+1) =  max ( 0, beta(i, iter) - 2 * ( 1 - counter(1, i) ) ) ;
end    
pconsumed = 0 ;
 for m= 1 : M
     temp1 = Carrier{1, m} ;
     for k =1 : k
         for i = 1 : nF
             for j = 1: nF
                 if temp1(i, j, k) == 1
                     pconsumed = pconsumed + PSR(k, i, m) + PRD(k, j, m);
                 end
             end
         end
     end
 end
landa(1, iter +1) =max(0, landa(1, iter) - 2 *  (PTotal - pconsumed) ) ;

if landa(1, iter +1) == 0
    landa(1, iter +1) = landa(1, iter) - 0.1 ; 
end

end

Goodput =  0;
for m =1:M
    temp1 = Carrier{1, m} ;
    temp2 = Rate{1,m} ;
      for k =1 : k
         for i = 1 : nF
             for j = 1: nF
                 if temp1(i, j, k) == 1
                     Goodput = Goodput+ sum(sum(sum(temp2))) ;
                 end
             end
         end
      end
 
end





