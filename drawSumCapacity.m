%
% Draw sum capacity of network with two femtocells
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Initialize MBS and MUE
PBS = 50 ; %dBm
MBS = BaseStation(0 , 0 , PBS); % (0,0) is the location and P_MBS = 50 dBm
mue = UE(-200, 0);
%% Initialize FBSs
%Generate fbsCount=16 FBSs
FBS_Max = cell(1,16);
for i=1:3
%     if i<= fbsCount
        FBS_Max{i} = FemtoStation(180+(i-1)*35,150, MBS, mue, 10);
%     end
end
FBS = cell(1,2);

FBS{1} = FBS_Max{1};
FBS{2} = FBS_Max{3};

p = linspace(0,1,100);
% figure;
% hold on;
% grid on;
for i=0:100
    for j=0:100
        p1 = i*0.01;
        p2 = j*0.01;
        gF_1 = fading_FBS_FUE(FBS{1},100);
        gF_2 = fading_FBS_FUE(FBS{2},100);
        I_1 = Interference_MBS(FBS{1}, MBS, -120, 100);
        I_2 = Interference_MBS(FBS{2}, MBS, -120, 100);
        C = log2(1+(p1*gF_1)/(I_1)) + log2(1+(p2*gF_2)/(I_2));
        Z(i+1,j+1) = C;
    end
end
[X,Y] = meshgrid(0:0.01:1);
mesh(X,Y,Z);
% for p1=p
%     for p2=p
%         gF_1 = fading_FBS_FUE(FBS{1},100);
%         gF_2 = fading_FBS_FUE(FBS{2},100);
%         I_1 = Interference_MBS(FBS{1}, MBS, -120, 100);
%         I_2 = Interference_MBS(FBS{2}, MBS, -120, 100);
%         C = log2(1+(p1*gF_1)/(I_1)) + log2(1+(p2*gF_2)/(I_2));
%         plot3(p1,p2,C);
%     end
% end