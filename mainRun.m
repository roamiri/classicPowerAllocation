%
% Main Run for power allocation in sparse femtocell network
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear;
clc;
% close all;
numFBS = 16;
final = cell(1,numFBS);
for i =4:numFBS
    display(i);
    final{1,i} = PA3(i, 1000, 1e4, 0);
end
%%
for i=4:numFBS
    sum(1,i) = final{1,i}.sum;
    min(1,i) = final{1,i}.min;
    CMUE(1,i) = final{1,i}.CMUE;
end
figure;
hold on;
grid on;
box on;
% plot( ones(1,16)*1.25, '--k' );
plot(sum, '--*r', 'LineWidth',1,'MarkerSize',10);
% plot(R3, '--sb', 'LineWidth',1,'MarkerSize',10);
% plot(comSumFUE.myBeta, '--*g');
% plot(comSumFUE.myCombine, '--*m');
% plot(comSumFUE.thresh, '--*c');
title('Sum capacity of FUEs','FontSize',14, 'FontWeight','bold');
xlabel('FBS Numbers','FontSize',14, 'FontWeight','bold');
ylabel('Capacity(b/s/HZ)','FontSize',14, 'FontWeight','bold');
% legend({'RF1','RF2'},'FontSize',14, 'FontWeight','bold');
% saveas(gcf,sprintf('FUE_Number_%d.jpg', j))
xlim([4 16]);
%%
figure;
hold on;
grid on;
box on;
plot( ones(4,16)*2.0, '--k' );
plot(CMUE, '--*r', 'LineWidth',1,'MarkerSize',10);
title('Capacity of MUE','FontSize',14, 'FontWeight','bold');
xlabel('FBS Numbers','FontSize',14, 'FontWeight','bold');
ylabel('Capacity(b/s/HZ)','FontSize',14, 'FontWeight','bold');
xlim([4 16]);
%%
figure;
hold on;
grid on;
box on;
plot( ones(4,16)*5, '--k' );
plot(min, '--*r', 'LineWidth',1,'MarkerSize',10);
title('Minimum Capacity of FUEs','FontSize',14, 'FontWeight','bold');
xlabel('FBS Numbers','FontSize',14, 'FontWeight','bold');
ylabel('Capacity(b/s/HZ)','FontSize',14, 'FontWeight','bold');
xlim([4 16]);
% ylim([0 7]);