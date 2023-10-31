function D_set = cal_diffusion_coef(tracks)
% tracks: NX4 [T, SD, SD_1, SD_2]
[~,L]=size(tracks);
t_set=tracks(:,1);
MSD_set=tracks(:,2:end);
D_set=zeros(1,L-1);
 figure
for i =2:L
    p= polyfit(t_set, MSD_set(:,i-1),1);
    %% plot
     hold on;
    scatter(t_set, MSD_set(:,i-1),'o')
    x0=0:max(t_set)/1e3:max(t_set);
    y0= polyval(p,x0);
    plot(x0,y0,'--')

    D_set(i-1) = p(1)/6;
end
% change the unit of diffusion cefficient from angstrom^2/s to m^2/s
% D_set=D_set.*1e-20;
end