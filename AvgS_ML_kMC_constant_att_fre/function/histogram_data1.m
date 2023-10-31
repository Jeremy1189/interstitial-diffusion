function h=histogram_data1(data,bin_nums,x_label,y_label,str,Facecolor,edges,Alpha)
if isempty(data)
    data=1e-2.*ones(1,5);
elseif length(unique(data))==1
    data(1:2)=repmat(unique(data),[1,2])-1e-3;
else
    disp('normal')
end
% figure('OuterPosition',[450,450,450,450])
h=histogram(data,bin_nums,'EdgeColor','none',"BinEdges",edges,...
        "EdgeAlpha",0,'FaceAlpha',Alpha,"FaceColor",Facecolor);
xlabel(x_label);
ylabel(y_label);
title(str)
set(gca,'FontSize',12,'Fontname', 'Arial','FontWeight','bold',"LineWidth",2)
