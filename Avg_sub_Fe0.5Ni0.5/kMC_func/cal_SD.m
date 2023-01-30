function [SD, SD_element] = cal_SD(Rt,R0, cross_boundary_record, interstitial_per_types, boundary)
%CAL_MSD Summary of this function goes here
L_dim= boundary(:,2)-boundary(:,1);
Dt =(Rt-R0);
Dt=Dt+ cross_boundary_record.* repmat(L_dim',size(cross_boundary_record,1),1);
sum_Dt=sum(Dt.^2,2);%
% MSD =sum(sum_Dt)./N;
SD =sum(sum_Dt);
% index of each element
element_type_set = unique(interstitial_per_types);
type_index_cell=cell(1,length(element_type_set));
for num_ratio_type=1:length(element_type_set)
    curent_index= find(interstitial_per_types==element_type_set(num_ratio_type));
    type_index_cell{num_ratio_type}=curent_index;
end
%partitial MSD
SD_element=zeros(1,length(type_index_cell));
for num_ratio_types=1:length(type_index_cell)
    current_index=type_index_cell{num_ratio_types};
    if ~isempty(current_index)
        sum_cur=sum_Dt(current_index);
%         MSD_element(num_ratio_types)=sum(sum_cur)./length(current_index);
        SD_element(num_ratio_types)=sum(sum_cur);
    else
        SD_element(num_ratio_types)=0;
    end
end
end

