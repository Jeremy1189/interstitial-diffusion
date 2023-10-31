function  plot_histog_edge(edges,s,h,color,LW)
    x=sort([edges,edges+ s],'ascend');
    Y=[h.Values,0;h.Values,0;];
    y=Y(:);
    %          plot(x,y,'-','Color',[231,41,138]./255,'LineWidth',LW)
    plot(x,y,'-','Color',color,'LineWidth',LW)
end