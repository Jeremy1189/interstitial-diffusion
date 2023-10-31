clc;
clear; 
close all;
% addpath(genpath(pwd));
load ../Interstitial_all_ratio_900K.mat ratio f_cell JF_cell k_tot_avg_cell

for i=1:length(f_cell)
    f(i)=f_cell{i}(1);
    JF(i)=JF_cell{i}(1);
    k_tot_avg= k_tot_avg_cell{i};
    JF(i)= 1/mean(1./k_tot_avg);
end

figure('OuterPosition',[450,550, 400,410])
plot(ratio(:,2),f,'-o','LineWidth',2,'Color','#1b9e77')
xlabel('Fe_xNi_{1-x}')
ylabel('f_{tr}')
legend('ML-E_m')
set(gca,'FontSize',12,'Fontname', 'Arial','FontWeight','bold');
saveas(gcf,'f_ML_Em.fig')
%%
figure('OuterPosition',[150,550, 400,410])
semilogy(ratio(:,2),JF,'-o','LineWidth',2,'Color','#1b9e77')
% hold on;
% semilogy(ratio(:,2),JumFre,'-o','LineWidth',2,'Color','r')
xlabel('Fe_xNi_{1-x}')
ylabel('\nu')
legend('ML-E_m')
set(gca,'FontSize',12,'Fontname', 'Arial','FontWeight','bold');
saveas(gcf,'v_ML_Em.fig')