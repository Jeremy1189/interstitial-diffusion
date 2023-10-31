function [central1_ID_nn_set,central1_ID_count_set,relative_central1_coords] ...
    = update_interstitial_nn_kmc(central1,per_coords,lattice_constant,supersize,largest_nn)
% Description:
%       This function aims to obtain the atomIDs in nn shells
% input:
%      central1_ID:size=[1,1]
%      central2_ID: size=[1,1]
%      per_coords: size=[n,3]

% output:
%       central1_central2_sortID_nn_set: 
%       central1_central2_nn_count_set
%       central1_ID_nn_set:
%       central1_ID_count_set:
% time:
%          2020/11/15 (first version)
%%
central1_coord = central1;
L = length(per_coords);
relative_central1_coords = per_coords - repmat(central1_coord,[L,1]);
Len = lattice_constant*supersize;
boundary = [0,Len;0,Len;0,Len];
relative_central1_coords=dis_boundary_check(relative_central1_coords,boundary);
% for num_coords=1:length(relative_central1_coords)
%     relative_central1_coords(num_coords,:)=boundary_check(relative_central1_coords(num_coords,:),boundary);
% % end
fix_point_trans=1e10;
% relative_central1_coords= round(relative_central1_coords.*fix_point_trans)./fix_point_trans;
% % periodic boundary check
% min_per_coords= min(per_coords);
% % max_per_coords= max(per_coords);
% % lattice_constant=3.488;
% % supersize=10;
% max_per_coords=lattice_constant*supersize*ones(1,3);
% box_size =  max_per_coords-min_per_coords;
% upper_boundary = round(1/2.*box_size.*fix_point_trans)./fix_point_trans;
% for i= 1:length(box_size)
%     central1_index_modify= find(relative_central1_coords(:,i)>upper_boundary(i));
%     relative_central1_coords(central1_index_modify,i)= relative_central1_coords(central1_index_modify,i)-box_size(i);
%     % less than and equal to the box_size
%     central1_index_modify= find(relative_central1_coords(:,i)<=-upper_boundary(i));
%      relative_central1_coords(central1_index_modify,i)= relative_central1_coords(central1_index_modify,i)+box_size(i);
% end
% calculated the distance
central1_nn_distance = sqrt( sum(relative_central1_coords.^2,2) );
% fixed
% fix_point_trans=1e10;
central1_nn_distance = round(central1_nn_distance.*fix_point_trans)./fix_point_trans;
%unique distance
%central1
unique_central1_dis_set= unique(central1_nn_distance);
% largest_nn=10;
central1_ID_nn_set =zeros(largest_nn,1);% output1
central1_ID_count_set= zeros(largest_nn,1);% remove central1ancy itself distance 
count_central1=0;
for num_unique = 1:largest_nn
    %central1
    index_central1_set=find(central1_nn_distance ==unique_central1_dis_set(num_unique));
    central1_ID_count_set(num_unique) = length(index_central1_set);
    central1_ID_nn_set(count_central1+1 : count_central1+ length(index_central1_set) ) = index_central1_set;%output2
    count_central1 = count_central1 + length(index_central1_set);
end
end

