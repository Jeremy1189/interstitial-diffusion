clc;
clear;
close all;
% load Interstitial_all_ratio_500K_no_dif.mat tracks_cell T D_cell ratio
% load test_code_all_ratio_900K_no_dif.mat tracks_cell T D_cell ratio
% load test_code_all_ratio_900K_2w.mat tracks_cell T D_cell ratio
load ../Interstitial_all_ratio_500_1200K_30w.mat tracks_cell T D_cell ratio  NN1_mig_ID_count_cell...
      NN1_mig_energy_cell interstial_struct_cell mig_index_cell
plot_flag=1;
jump_prob_cell= cell(length(ratio),3);
% for i=1:length(ratio)
figure('OuterPosition',[150,450,450,450])
y_label='Jump Frequency';
shapes='ovphsdv*<^>';
x=1:8;
for i=10:-1:2    
       NN1_mig_ID_set=NN1_mig_ID_count_cell{i};
       NN1_mig_energy=NN1_mig_energy_cell{i};
        tracks=tracks_cell{i};
        per=interstial_struct_cell{i};
        types=per(:,2);
        L_tr=length(tracks);
        I1  =tracks(1:L_tr,2);%back
        I1_types_set= types(I1); 
        I2  =tracks(1:L_tr,6);   % mig
        I2_types_set= types(I2); 
       
       
       mig_index=mig_index_cell{i};
       mig_energy=zeros(length(mig_index),1);
       mig_ID_set=zeros(length(mig_index),1);
       for s=1:length(mig_index)
           mig_ID_set(s)=NN1_mig_ID_set(s,mig_index(s));
           mig_energy(s)=NN1_mig_energy(s,mig_index(s)); 
           if mig_index(s)<=4
              temp_type= I1_types_set(s);
              I1_types_set(s) = I2_types_set(s);
              I2_types_set(s) = temp_type;
           end
       end
  
        I3_types_set= types(mig_ID_set);        

        types_pair_set=[I1_types_set,I2_types_set,I3_types_set];
        EBF=mig_energy;
        index_Ni_Ni_Ni = find(types_pair_set(:,1)==1 & types_pair_set(:,2)==1 & types_pair_set(:,3)==1 );        
        index_Ni_Fe_Ni=find(types_pair_set(:,1)==1 & types_pair_set(:,2)==2 & types_pair_set(:,3)==1 );
        index_Ni_Ni_Fe=find(types_pair_set(:,1)==1 & types_pair_set(:,2)==1 & types_pair_set(:,3)==2);
        index_Fe_Ni_Ni=find(types_pair_set(:,1)==2 & types_pair_set(:,2)==1 & types_pair_set(:,3)==1 );
        index_Ni_Fe_Fe=find(types_pair_set(:,1)==1 & types_pair_set(:,2)==2 & types_pair_set(:,3)==2 );
        index_Fe_Fe_Fe=find(types_pair_set(:,1)==2 & types_pair_set(:,2)==2 & types_pair_set(:,3)==2 );
        index_Fe_Ni_Fe=find(types_pair_set(:,1)==2 & types_pair_set(:,2)==1 & types_pair_set(:,3)==2 );
        index_Fe_Fe_Ni=find(types_pair_set(:,1)==2 & types_pair_set(:,2)==2 & types_pair_set(:,3)==1 );
        jump_prob_set=zeros(8,1);
        jump_prob_set(1)=length(index_Ni_Ni_Ni)/length(EBF);
        jump_prob_set(2)=length(index_Ni_Fe_Ni)/length(EBF);
        jump_prob_set(3)=length(index_Ni_Ni_Fe)/length(EBF);
        jump_prob_set(4)=length(index_Fe_Ni_Ni)/length(EBF);
        jump_prob_set(5)=length(index_Ni_Fe_Fe)/length(EBF);
        jump_prob_set(6)=length(index_Fe_Fe_Fe)/length(EBF);
        jump_prob_set(7)=length(index_Fe_Ni_Fe)/length(EBF);
        jump_prob_set(8)=length(index_Fe_Fe_Ni)/length(EBF);
     
        str_name=['Ni_{',num2str(ratio(i,1)),'}Fe_{',num2str(ratio(i,2)),'}'];
        jump_prob_cell{i,1}=str_name;
        legend_set={'Ni-Ni-Ni','Ni-Fe-Ni','Ni-Ni-Fe','Fe-Ni-Ni','Ni-Fe-Fe','Fe-Fe-Fe','Fe-Ni-Fe',...
            'Fe-Fe-Ni'};
        jump_prob_cell{i,2}=legend_set;
        jump_prob_cell{i,3}=jump_prob_set;
        plot(x,jump_prob_set,['--',shapes(i)],'LineWidth',2);
        hold on;
end
xticks(x);
%xticklabels({'Ni-Ni-Ni','Ni-Fe-Ni','Ni-Ni-Fe','Fe-Ni-Ni','Ni-Fe-Fe','Fe-Fe-Fe','Fe-Ni-Fe',...
%     'Fe-Fe-Ni'});
legend(jump_prob_cell(10:-1:2,1))
xticklabels(legend_set);
xlabel('Interstitial Pairs')
ylabel(y_label);
set(gca,'FontSize',12,'Fontname', 'Arial','FontWeight','bold')
saveas(gcf,'Jump_prob_ML_kMC.fig')
save Jump_prop_ML_kMC.mat jump_prob_set

% figure('OuterPosition',[350,450, 450,450])
% plot(ratio(:,2), D_set*0.72,'-o','Color','b',"LineWidth",2);
% xlabel('Fe_xNi_{1-x}')
% %         ylabel('MSD ($\AA$)','interpreter','latex')
% ylabel('D^*(m^2/s)')
% title(['NiFe-',num2str(T),'K'])
% set(gca,'FontSize',12,'Fontname', 'Arial','FontWeight','bold',"LineWidth",2,'YScale','log')
