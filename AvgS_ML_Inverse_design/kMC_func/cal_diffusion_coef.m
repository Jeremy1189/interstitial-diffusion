function D_set = cal_diffusion_coef(tracks)
% tracks: NX4 [T, SD, SD_1, SD_2]
t_set=tracks(:,1);
MSD_set=tracks(:,2);
p=polyfit(t_set,MSD_set,1);
D_set = p(1)/6;
if D_set<=0
    D_set=1e20;
end
% change the unit of diffusion cefficient from angstrom^2/s to m^2/s
D_set=D_set.*1e-20;

end