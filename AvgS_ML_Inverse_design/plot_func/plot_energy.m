
figure;

x=0:1:9;
y=energy_count_set(1:10,:);
plot(x,y(:,1),'-o','Color','#d95f02', 'LineWidth',2);
hold on;
plot(x,y(:,2),'-s','Color','#1b9e77', 'LineWidth',2);
plot(x,y(:,3),'-v','Color','#7570b3', 'LineWidth',2);

xlabel('timesteps(*1e4)')
ylabel('energy')
legend('perfect','vac-initial','vac-final','box','off')
title([num2str(T),'K'])
set(gca,'FontSize',12,'FontName','Arial','LineWidth',2,'FontWeight','b')

