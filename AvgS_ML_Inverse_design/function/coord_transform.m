function transformed_coords =coord_transform(interstitial_coords,initia_I1_ID_back,initia_I2_ID_mig,...
    final_mig_I3_ID,lattice_constant,supersize)

%x-axis
I1_coord=interstitial_coords(initia_I1_ID_back,:);
I2_coord=interstitial_coords(initia_I2_ID_mig,:);
L=supersize*lattice_constant;
boundary=[0,L;0,L;0,L];
dif_coords1=I2_coord-I1_coord;
dif_coords1=dis_boundary_check(dif_coords1, boundary);%
%  dif_coords1 = boundary_check(dif_coords1,boundary);% x-axis

%origin
central_position_mig_former=I1_coord+dif_coords1/2;% the original position
central_position_mig_former=boundary_check(central_position_mig_former, boundary);
% x-axis
x_axis_temp= I2_coord-central_position_mig_former;
x_axis_temp= dis_boundary_check(x_axis_temp,boundary);
x_axis=x_axis_temp./norm(x_axis_temp);% the distrance between two interstitial
% central_position_mig_former=I1_coord;
% y-axis
I3_coord= interstitial_coords(final_mig_I3_ID,:);
dif_coords2=I3_coord-central_position_mig_former;
vector_temp= dis_boundary_check(dif_coords2, boundary);% 
y_axis_temp=cross(x_axis,vector_temp);%the distrance between one interstitial and centre
y_axis= y_axis_temp./norm(y_axis_temp);
%z-axis
z_axis= cross(x_axis,y_axis);

% relative coordinates
origin_coord = central_position_mig_former;
interstitial_str_coords=interstitial_coords;
L = length(interstitial_str_coords);
relative_coords = interstitial_str_coords - repmat(origin_coord,[L,1]);
relative_coords= dis_boundary_check(relative_coords, boundary);% 
% for i=1:length(relative_coords)
%     relative_coords(i,:) = boundary_check(relative_coords(i,:),boundary);% x-axis
% end

fix_point_trans=1e10;
relative_coords= round(relative_coords.*fix_point_trans)./fix_point_trans;

%% the projection length vector
x_new =relative_coords*x_axis';
y_new =relative_coords*y_axis';
z_new =relative_coords*z_axis';
transformed_coords= [x_new,y_new,z_new];
end
