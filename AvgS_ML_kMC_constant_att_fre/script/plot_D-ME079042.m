clc;
clear;
close all;
% load Interstitial_all_ratio_500K_no_dif.mat tracks_cell T D_cell ratio
% load test_code_all_ratio_900K_no_dif.mat tracks_cell T D_cell ratio
% load test_code_all_ratio_900K_2w.mat tracks_cell T D_cell ratio
% load Interstitial_all_ratio_1100K_20e4.mat tracks_cell T D_cell ratio
load ../Interstitial_all_ratio_900K.mat tracks_cell T D_cell ratio
plot_flag=1;
for i=1:length(ratio)
%    T=T_set(i);
% for i=5
if plot_flag==1
%     figure('OuterPosition',[350+i*100,450, 450,450])
         tracks=tracks_cell{i};
        L_tr=length(tracks);
        t_set  =tracks(1:L_tr,1);
        MSD_set  =tracks(1:L_tr,10);       
        x=t_set;
%         figure;
%         scatter(x,MSD_set,'.')
        P=polyfit(x,MSD_set,1);
        a=0:max(x)/1e3:max(x);
        f = polyval(P,a);
%         hold on;
%         plot(a,f,'r')
        D=P(1)./6;
        D_set(i)=D*4.001e-17;
        disp(D)
%         xlabel('time (s)')
% %         ylabel('MSD ($\AA$)','interpreter','latex')
%          ylabel('MSD ')
%         title(['Fe_{',num2str(ratio(i,2)),'}Ni_{',num2str(ratio(i,1)),'}-T=',num2str(T),'K'])
%         set(gca,'FontSize',12,'Fontname', 'Arial','FontWeight','bold',"LineWidth",2)
       
end  
% D_set(i)=D_cell{i}(1);
end
fold = 1;
figure('OuterPosition',[350,450, 450,450])
plot(ratio(:,2), D_set.*fold,'-o','Color','b',"LineWidth",2);
xlabel('Fe_xNi_{1-x}')
%         ylabel('MSD ($\AA$)','interpreter','latex')
ylabel('D^*(m^2/s)')
title(['NiFe-',num2str(T),'K'])
set(gca,'FontSize',12,'Fontname', 'Arial','FontWeight','bold',"LineWidth",2,'YScale','log')
saveas(gcf,['NiFe-',num2str(T),'K.fig'])
