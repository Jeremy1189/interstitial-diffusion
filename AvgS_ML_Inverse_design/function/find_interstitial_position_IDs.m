function [ID_1,postions_1,NN1_ID_1,ID_2,postions_2,NN1_ID_2] = find_interstitial_position_IDs(interstitial_structure,lattice_constant,super_size_length,largest_nn)
%FIND_INTERSTITIAL_POSITION_IDS Summary of this function goes here
%   Detailed explanation goes here
%input 
% % interstitial_structure: [N,5], N atoms, 5: [IDs, types, x,y,z]

% Output
% ID_1: dumbbell interstitial_1  ID
% position_1: dumbbell interstitial_1  postion (x ,y ,z)
% ID_2: dumbbell interstitial_2 ID
% position_2: dumbbell interstitial_2  postion (x ,y ,z)

%%
IDs=interstitial_structure(:,1);
% types=interstitial_structure(:,2);
per_coords= interstitial_structure(:,3:5);
L=size(per_coords,1);
count=0;
ID_index=zeros(100,1);
positions=zeros(100,3);
NN1_ID_set=zeros(100,4);
for i=1:L
    cur_coords= per_coords(i,:);
    
    [central1_ID_nn_set,central1_ID_count_set,~] ...
        = update_interstitial_nn_kmc(cur_coords,per_coords,lattice_constant,super_size_length,largest_nn+1);% the nearest atom is itself
    NN1_IDs=central1_ID_nn_set(2:central1_ID_count_set(2)+1);% 2:central1_ID_count_set(2)+1 is the total 1st nn atoms
    if length(NN1_IDs)==4
        count=count+1;
        ID_index(count)=i;
        positions(count,:)=cur_coords;
        NN1_ID_set(count,:)=NN1_IDs';
    end
     
end
   if count>2
      
      disp('There are ',num2str(count/2),' dumbell interstitials!')
   end  
   ID_temp=IDs(ID_index(1:count));
   positions_temp=positions(1:count,:);
   ID_1= ID_temp(1);
   postions_1=positions_temp(1,:);
   ID_2=ID_temp(2);
   postions_2=positions_temp(2,:);
   NN1_ID_1=NN1_ID_set(1,:);
   NN1_ID_2=NN1_ID_set(2,:);
end

