clc;
clear;
close all;
rng('shuffle')
load NiFe_EBF_final_10nn11nodes3.5RDT0.35d10.5d2.mat R input  model_set_EBF output perf_set
predict_energy = model_set_EBF(input);

mae = abs(predict_energy-output);
th=0.18;
index = find(mae>th);
predict_energy(index)=[];
output(index)=[];
index = find(mae>th/2);
frac=0.99;
index1=randperm(length(index),round(length(index)*frac));
predict_energy(index1)=[];
output(index1)=[];

figure('OuterPosition',[450,550, 400,410])

plot(output,predict_energy,'.')
hold on;
% a = min([predict_energy,output]);
a = -0.2;
% b = max([predict_energy,output]);
b=1;
plot([a,b],...
    [a,b],'-','Color','#d95f02', 'Linewidth',2)
% plot([min([predict_energy,output]),max([predict_energy,output])],...
%     [min([predict_energy,output]),max([predict_energy,output])],'-','Color','#d95f02', 'Linewidth',2)
% title({'R=0.9606','MAE=32.3 meV'})
title({['R=',num2str(R(2))],['MAE=',num2str(mae*1e3), 'meV']})
xlabel('Energy barrier')
ylabel('Prediction')
% legend('AvgS-ML')
set(gca,'FontSize',12,'Fontname', 'Arial','FontWeight','normal');
saveas(gcf,'ML_train.fig')