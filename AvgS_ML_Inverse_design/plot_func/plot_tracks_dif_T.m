clc;
clear;
close all;
% load dif_m_all_ratio_900K_no_dif.mat tracks_cell T D_cell
load test_code_all_ratio_900K_no_dif.mat tracks_cell T D_cell

for i=1:length(D_cell)
     T=T_set(i);
    figure('OuterPosition',[350,450, 450,450])
         tracks=tracks_cell{i};
        L_tr=length(tracks);
        t_set  =tracks(1:L_tr,1);
        MSD_set  =tracks(1:L_tr,10)*1e-20;       
        x=t_set;
%         figure;
        scatter(x,MSD_set,'.')
        P=polyfit(x,MSD_set,1);
        a=0:max(x)/1e3:max(x);
        f = polyval(P,a);
        hold on;
        plot(a,f,'r')
        D=P(1)./6;
        disp(D)
        xlabel('time (s)')
%         ylabel('MSD ($\AA$)','interpreter','latex')
         ylabel('MSD (m^2)')
        title(['T=',num2str(T),'K'])
        set(gca,'FontSize',12,'Fontname', 'Arial','FontWeight','bold',"LineWidth",2)
    
end
  figure('OuterPosition',[350,450, 450,450])
  plot(1000./T_set, D_cell,'-o','Color','b',"LineWidth",2);
   xlabel('1000/T (1/K)')
%         ylabel('MSD ($\AA$)','interpreter','latex')
         ylabel('log(D^*) (m^2)')
        title(['NiFe-1e15'])
        set(gca,'FontSize',12,'Fontname', 'Arial','FontWeight','bold',"LineWidth",2,'YScale','log')

