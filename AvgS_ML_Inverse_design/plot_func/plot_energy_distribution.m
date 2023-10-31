clc;clear;close all;
load Interstitial_all_ratio_900K_no_dif.mat  NN1_mig_energy_cell NN1_mig_ID_count_cell...
    interstial_struct_cell ratio tracks_cell
ratio_set=ratio(:,1);

N_size=20000;
si=58000-1;
fi=si+N_size;
num_s=2;
tracks=tracks_cell{num_s}(si:fi,:);
interstitial_IDs=tracks(:,[2,6]);
interstitial_struct = interstial_struct_cell{num_s};
per_types=interstitial_struct(:,2);
back_mig_types=per_types(interstitial_IDs(:,1));
mig_mig_types=per_types(interstitial_IDs(:,2));
NN1_mig_ID=NN1_mig_ID_count_cell{num_s}(si:fi,:);
NN1_types_set=per_types(NN1_mig_ID);
back_NN1_types_set=NN1_types_set(:,1:4);
mig_NN1_types_set=NN1_types_set(:,5:8);
% NN1_types_set=vac_mig_NN1_type_set{1}(si:fi,3:end);
% new-interstitial
%migration atom divid
back_mig_Ni_index=find(back_mig_types==1);
back_mig_Fe_index=find(back_mig_types==2);
back_Ni_mig_NN1_types= back_NN1_types_set(back_mig_Ni_index,:);
back_Fe_mig_NN1_types= back_NN1_types_set(back_mig_Fe_index,:);
%% NN1 dvidid
back_NN1_Ni_index=find(back_Ni_mig_NN1_types==1);
back_NN1_Fe_index=find(back_Ni_mig_NN1_types==2);

mig_NN1_Ni_index=find(back_Fe_mig_NN1_types==1);
mig_NN1_Fe_index=find(back_Fe_mig_NN1_types==2);
%% energy
ML_NN1_mig_energy =NN1_mig_energy_cell{num_s}(si:fi,:);
back_ML_NN1_mig_energy= ML_NN1_mig_energy(:,1:4);
mig_ML_NN1_mig_energy= ML_NN1_mig_energy(:,5:8);


back_Ni_mig_NN1_energy =  back_ML_NN1_mig_energy(back_mig_Ni_index,:);
back_Fe_mig_NN1_energy =  back_ML_NN1_mig_energy(back_mig_Fe_index,:);


ML_NN1_energy_Ni1 = back_Ni_mig_NN1_energy(back_NN1_Ni_index);
ML_NN1_energy_Fe1 = back_Ni_mig_NN1_energy(back_NN1_Fe_index);

ML_NN1_energy_Ni2 = back_Fe_mig_NN1_energy(mig_NN1_Ni_index);
ML_NN1_energy_Fe2 = back_Fe_mig_NN1_energy(mig_NN1_Fe_index);


% ML_NN1_energy_Cr =  ML_NN1_mig_energy(Cr_index);
Alpha=0.5;
Alpha1=0.7;
s=0.03;
edges=0:s:1.5;
edges1=0.5:0.08:1.8;
FontSize=12;
LW=2;
N=10e4;
figure('OuterPosition',[350+0*100,450, 450,450])
%     Figure1=figure;
    h0=histogram(back_ML_NN1_mig_energy,"BinEdges",edges,...
        "EdgeAlpha",0,'FaceAlpha',Alpha,"FaceColor",[231,41,138]./255)
  title(['Fe_{',num2str(1-ratio_set(num_s)),'}Ni_{',num2str(ratio_set(num_s)),'}']);
%     title(['1350K-NN1-start-',num2str(si)],'FontSize',FontSize,"FontWeight","bold")
    xlabel('Energy barrier (eV)')
    ylabel('Count')
    set(gca,'FontSize',FontSize, 'FontWeight' ,"bold",'FontName','Arial','LineWidth',LW)


figure('OuterPosition',[350+1*100,450, 450,450])
%     Figure1=figure;
    h1=histogram(ML_NN1_energy_Ni1,"BinEdges",edges,...
        "EdgeAlpha",0,'FaceAlpha',Alpha,"FaceColor",[231,41,138]./255);
    hold on
%     figure
    h2=histogram(ML_NN1_energy_Fe1,"BinEdges",edges,...
        "EdgeAlpha",0,'FaceAlpha',Alpha,"FaceColor",[117,112,179]./255);
%     figure
%     h3=histogram(ML_NN1_energy_Cr,"BinEdges",edges,...
%         "EdgeAlpha",0,'FaceAlpha',Alpha1,"FaceColor",[27,158,119]./255); 
   
    x=sort([edges,edges+ s],'ascend');
    Y=[h1.Values,0;h1.Values,0;];
    y=Y(:);
    %          plot(x,y,'-','Color',[231,41,138]./255,'LineWidth',LW)
    plot(x,y,'-','Color',[231,41,138]./255,'LineWidth',LW)
    Y=[h2.Values,0;h2.Values,0;];
    y=Y(:);
    %          plot(x,y,'-','Color',[117,112,179]./255,'LineWidth',LW)
    plot(x,y,'-','Color',[117,112,179]./255,'LineWidth',LW)
%     Y=[h3.Values,0;h3.Values,0;];
%     y=Y(:);
%     plot(x,y,'-','Color',[27,158,119]./255,'LineWidth',LW)
    xlabel('Energy barrier (eV)')
    ylabel('Count')
  
%         legend('Ni','Fe','Cr','box','off')
        legend('Ni-Ni','Ni-Fe','box','off')
  
    %     title(['The migration type shares in ',Ratio{j}])
     title(['Fe_{',num2str(1-ratio_set(num_s)),'}Ni_{',num2str(ratio_set(num_s)),'}']);
%     title(['1350K-NN1-start-',num2str(si)],'FontSize',FontSize,"FontWeight","bold")
    %
    set(gca,'FontSize',FontSize, 'FontWeight' ,"bold",'FontName','Arial','LineWidth',LW)
    
  
    %%
    figure('OuterPosition',[350+2*100,450, 450,450])
%     Figure1=figure;
    h1=histogram(ML_NN1_energy_Ni2,"BinEdges",edges,...
        "EdgeAlpha",0,'FaceAlpha',Alpha,"FaceColor",[231,41,138]./255);
    hold on
%     figure
    h2=histogram(ML_NN1_energy_Fe2,"BinEdges",edges,...
        "EdgeAlpha",0,'FaceAlpha',Alpha,"FaceColor",[117,112,179]./255);
%     figure
%     h3=histogram(ML_NN1_energy_Cr,"BinEdges",edges,...
%         "EdgeAlpha",0,'FaceAlpha',Alpha1,"FaceColor",[27,158,119]./255); 
   
    x=sort([edges,edges+ s],'ascend');
    Y=[h1.Values,0;h1.Values,0;];
    y=Y(:);
    %          plot(x,y,'-','Color',[231,41,138]./255,'LineWidth',LW)
    plot(x,y,'-','Color',[231,41,138]./255,'LineWidth',LW)
    Y=[h2.Values,0;h2.Values,0;];
    y=Y(:);
    %          plot(x,y,'-','Color',[117,112,179]./255,'LineWidth',LW)
    plot(x,y,'-','Color',[117,112,179]./255,'LineWidth',LW)
%     Y=[h3.Values,0;h3.Values,0;];
%     y=Y(:);
%     plot(x,y,'-','Color',[27,158,119]./255,'LineWidth',LW)
    xlabel('Energy barrier (eV)')
    ylabel('Count')
  
%         legend('Ni','Fe','Cr','box','off')
        legend('Fe-Ni','Fe-Fe','box','off')
  
    %     title(['The migration type shares in ',Ratio{j}])
     title(['Fe_{',num2str(1-ratio_set(num_s)),'}Ni_{',num2str(ratio_set(num_s)),'}']);
%     title(['1350K-NN1-start-',num2str(si)],'FontSize',FontSize,"FontWeight","bold")
    %
    set(gca,'FontSize',FontSize, 'FontWeight' ,"bold",'FontName','Arial','LineWidth',LW)