function plot_EBF_RDT_d1_d2(flag,EBF,RDT,dis_set)

if flag==1
%         close all;
        figure;
        histogram(EBF,40)  
        title('EBF')
%         histogram(RDT,40)
        figure;
%         histogram(EBF,40)  
        histogram(RDT,40)
        title('RDT')
        figure;
        histogram(dis_set(:,1),100)
        hold on;
        histogram(dis_set(:,2),100)
        histogram(dis_set(:,3),100)
        histogram(dis_set(:,4),100)
        title('d1')
        
        
%       
        figure;
        
        histogram(dis_set(:,5),100)
        hold on;
        histogram(dis_set(:,6),100)
         title('d2')
end



end