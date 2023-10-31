load NiCoCr_ML_1350k.mat mig_type_set mig_index_set NN1_mig_energy_set1
N_size=2000;
si=1;
fi=si+N_size;
mig_energy_set=NN1_mig_energy_set1(i,mig_index_set(i));
count=0;
for i=si:fi
    count=count+1;
    mig_energy_set(count)=NN1_mig_energy_set1(i,mig_index_set(i));
end
% for i=1:length(mig_type_set)
%     mig_energy_set(i)=NN1_mig_energy_set1(i,mig_index_set(i));
% end
mig_type=mig_type_set(si:fi,2);
Ni_index=find(mig_type==1);
Co_index=find(mig_type==2);
Cr_index=find(mig_type==3);

mig_energy_Ni =  mig_energy_set(Ni_index);
mig_energy_Co =  mig_energy_set(Co_index);
mig_energy_Cr =  mig_energy_set(Cr_index);

Alpha=0.5;
Alpha1=0.7;
s=0.005;
edges=0:s:1.5;
edges1=0.5:0.08:1.8;
FontSize=12;
LW=2;
% figure;
% 
% % for i=1:10000
% %     
% %     switch mig_type(i)
% %         case 1 
% %             scatter(i,mig_energy_set(i),'bo')
% %             hold on
% %         case 2 
% %             scatter(i,mig_energy_set(i),'rv')
% %             hold on
% %             case 3 
% %             scatter(i,mig_energy_set(i),'gs')
% %             hold on
% %     end
% % end
% legend('Ni','Co','Cr','box','off')
% set(gca,'FontSize',FontSize, 'FontWeight' ,"bold",'FontName','Arial','LineWidth',LW)
% figure; 
% N=5000;
% index_1=find(Ni_index>N,1);
% index_2=find(Co_index>N,1);
% index_3=find(Cr_index>N,1);
% scatter(Ni_index(1:index_1),mig_energy_Ni(1:index_1),'bo');
% hold on;
% plot(Co_index(1:index_2),mig_energy_Co(1:index_2),'rv');
% plot(Cr_index(1:index_3),mig_energy_Cr(1:index_3),'gs');
%  legend('Ni','Co','Cr','box','off')
%  set(gca,'FontSize',FontSize, 'FontWeight' ,"bold",'FontName','Arial','LineWidth',LW)
figure;
N=0.5e4;
%     Figure1=figure;
    h1=histogram(mig_energy_Ni,"BinEdges",edges,...
        "EdgeAlpha",0,'FaceAlpha',Alpha,"FaceColor",[231,41,138]./255);
    hold on
%     figure
    h2=histogram(mig_energy_Co,"BinEdges",edges,...
        "EdgeAlpha",0,'FaceAlpha',Alpha,"FaceColor",[117,112,179]./255);
%     figure
    h3=histogram(mig_energy_Cr,"BinEdges",edges,...
        "EdgeAlpha",0,'FaceAlpha',Alpha1,"FaceColor",[27,158,119]./255); 
   
    x=sort([edges,edges+ s],'ascend');
    Y=[h1.Values,0;h1.Values,0;];
    y=Y(:);
    %          plot(x,y,'-','Color',[231,41,138]./255,'LineWidth',LW)
    plot(x,y,'-','Color',[231,41,138]./255,'LineWidth',LW)
    Y=[h2.Values,0;h2.Values,0;];
    y=Y(:);
    %          plot(x,y,'-','Color',[117,112,179]./255,'LineWidth',LW)
    plot(x,y,'-','Color',[117,112,179]./255,'LineWidth',LW)
    Y=[h3.Values,0;h3.Values,0;];
    y=Y(:);
    plot(x,y,'-','Color',[27,158,119]./255,'LineWidth',LW)
   
  
        legend('Ni','Co','Cr','box','off')
  
    %     title(['The migration type shares in ',Ratio{j}])
    %     title(['Ni(',num2str(ratio_set(j)),')Co(',num2str(1-ratio_set(j)),')']);
    title(['1350K-Mig-E-start-',num2str(si)],'FontSize',FontSize,"FontWeight","bold")
    %
    set(gca,'FontSize',FontSize, 'FontWeight' ,"bold",'FontName','Arial','LineWidth',LW)