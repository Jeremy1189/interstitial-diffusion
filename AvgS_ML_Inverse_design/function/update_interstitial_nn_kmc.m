function [central1_ID_nn_set,central1_ID_count_set,relative_central1_coords] ...
    = update_interstitial_nn_kmc(central1,per_coords,boundary,largest_nn)
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
% Len = lattice_constant*supersize;
% boundary = [0,Len;0,Len;0,Len];
relative_central1_coords=dis_boundary_check(relative_central1_coords,boundary);
% for num_coords=1:length(relative_central1_coords)
%     relative_central1_coords(num_coords,:)=boundary_check(relative_central1_coords(num_coords,:),boundary);
% % end
fix_point_trans=1e6;
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

