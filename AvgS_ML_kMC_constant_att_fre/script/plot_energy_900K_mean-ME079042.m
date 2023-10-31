clc;
clear;
close all;
% load Interstitial_all_ratio_500K_no_dif.mat tracks_cell T D_cell ratio
% load test_code_all_ratio_900K_no_dif.mat tracks_cell T D_cell ratio
% load test_code_all_ratio_900K_2w.mat tracks_cell T D_cell ratio
load ../Interstitial_all_ratio_900K.mat tracks_cell T D_cell ratio  NN1_mig_ID_count_cell...
     NN1_mig_energy_cell interstial_struct_cell
plot_flag=1;
mu_sigma_cell= cell(length(ratio),3);

% for i=1:length(ratio)
for i=1:11

       NN1_mig_ID_set=NN1_mig_ID_count_cell{i};
       NN1_mig_energy=NN1_mig_energy_cell{i}';
        tracks=tracks_cell{i};
        per=interstial_struct_cell{i};
        types=per(:,2);
        L_tr=length(tracks);
        I1  =tracks(1:L_tr,2);
        I2  =tracks(1:L_tr,6);       
        I3_types_set= types(NN1_mig_ID_set');
        
        I1_set=repmat(types(I1),1,8);
        I1_set=I1_set';
        I2_set=repmat(types(I2),1,8);
        I2_set=I2_set';
         % the first 4 energy in NN1_energy is recorded the I1 migration, the remain 4 is the I2 migration
         % then, we exchange the 4 columes of I1 and I2 to make I2 is the
         % migration type matrix, I1 is the back types
        temp=I1_set(1:4,:);
        I1_set(1:4,:)=I2_set(1:4,:);
        I2_set(1:4,:)=temp;
        types_pair_set=[I1_set(:),I2_set(:),I3_types_set(:)];
        EBF=NN1_mig_energy(:);
        index_Ni_Ni_Ni = find(types_pair_set(:,1)==1 & types_pair_set(:,2)==1 & types_pair_set(:,3)==1 );        
        index_Ni_Fe_Ni=find(types_pair_set(:,1)==1 & types_pair_set(:,2)==2 & types_pair_set(:,3)==1 );
        index_Ni_Ni_Fe=find(types_pair_set(:,1)==1 & types_pair_set(:,2)==1 & types_pair_set(:,3)==2);
        index_Fe_Ni_Ni=find(types_pair_set(:,1)==2 & types_pair_set(:,2)==1 & types_pair_set(:,3)==1 );
        index_Ni_Fe_Fe=find(types_pair_set(:,1)==1 & types_pair_set(:,2)==2 & types_pair_set(:,3)==2 );
        index_Fe_Fe_Fe=find(types_pair_set(:,1)==2 & types_pair_set(:,2)==2 & types_pair_set(:,3)==2 );
        index_Fe_Ni_Fe=find(types_pair_set(:,1)==2 & types_pair_set(:,2)==1 & types_pair_set(:,3)==2 );
        index_Fe_Fe_Ni=find(types_pair_set(:,1)==2 & types_pair_set(:,2)==2 & types_pair_set(:,3)==1 );
        mu_sigma_set=zeros(8,1);
       
        count=0;
        count=count+1;
        mu_sigma_set(count) =  mean(EBF(index_Ni_Ni_Ni));

        count=count+1;
        mu_sigma_set(count) =  mean(EBF(index_Ni_Ni_Fe));
        count=count+1;
        mu_sigma_set(count) =  mean(EBF(index_Ni_Fe_Ni));
        count=count+1;
        mu_sigma_set(count) =  mean(EBF(index_Ni_Fe_Fe));
        count=count+1;
        mu_sigma_set(count) =  mean(EBF(index_Fe_Ni_Ni));
        count=count+1;
        mu_sigma_set(count) =  mean(EBF(index_Fe_Ni_Fe));
        count=count+1;
        mu_sigma_set(count) =  mean(EBF(index_Fe_Fe_Ni));
        count=count+1;
        mu_sigma_set(count) =  mean(EBF(index_Fe_Fe_Fe));
        str_name=['Ni_{',num2str(ratio(i,1)),'}Fe_{',num2str(ratio(i,2)),'}'];
        mu_sigma_cell{i,1}=str_name;
        legend_set={'Ni-Ni-Ni','Ni-Ni-Fe','Ni-Fe-Ni','Ni-Fe-Fe','Fe-Ni-Ni',...
            'Fe-Ni-Fe','Fe-Fe-Ni','Fe-Fe-Fe'};
        mu_sigma_cell{i,2}=legend_set;
        mu_sigma_cell{i,3}=mu_sigma_set;
       
end
save mean_ML_kMC_dif_ratio_900K.mat mu_sigma_cell

% figure('OuterPosition',[350,450, 450,450])
% plot(ratio(:,2), D_set*0.72,'-o','Color','b',"LineWidth",2);
% xlabel('Fe_xNi_{1-x}')
% %         ylabel('MSD ($\AA$)','interpreter','latex')
% ylabel('D^*(m^2/s)')
% title(['NiFe-',num2str(T),'K'])
% set(gca,'FontSize',12,'Fontname', 'Arial','FontWeight','bold',"LineWidth",2,'YScale','log')
