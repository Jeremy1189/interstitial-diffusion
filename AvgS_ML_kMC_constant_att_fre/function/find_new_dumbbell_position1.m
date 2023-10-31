function [interstitial_new_target_coord,interstitial_new_mig_coord] = ...
    find_new_dumbbell_position1(initerstitial_per_coords,target_interstitial_ID,mig_interstitial_ID,set_dis,boundary)
%FIND_NEW_DUMBBELL_POSITION Summary of this function goes here
%   Detailed explanation goes here
   %% there are six possiple position for forming a new dumbell
        %interstitiall pair
        migration_interstitial_coord= initerstitial_per_coords(mig_interstitial_ID,:);
        new_interstitial_old_coord = initerstitial_per_coords(target_interstitial_ID,:);
        
%         possible_pos=[new_interstitial_old_coord+[set_dis,0,0];...
%             new_interstitial_old_coord+[-set_dis,0,0];...
%             new_interstitial_old_coord+[0,set_dis,0];...
%             new_interstitial_old_coord+[0,-set_dis,0];...
%             new_interstitial_old_coord+[0,0,set_dis];...
%             new_interstitial_old_coord+[0,0,-set_dis]];
        possible_pos=[
            new_interstitial_old_coord+[0,set_dis,0];...
            new_interstitial_old_coord+[0,-set_dis,0];...
            ];

        % boundary check
        for s=1:size( possible_pos,1)
            possible_pos(s,:)=boundary_check( possible_pos(s,:),boundary);
        end
        dif_dis_set= possible_pos-repmat(migration_interstitial_coord,size(possible_pos,1),1);%6 possible positions
        dif_dis_set= dis_boundary_check(dif_dis_set, boundary);
        sum_dif_dis= sum(dif_dis_set.^2,2);
        %         sum_dis_set= sum(sum_dif_dis);
        [~,min_index]=min(sum_dif_dis);
        new_interstitial_new_coord=possible_pos(min_index,:);
        interstitial_new_mig_coord= boundary_check(new_interstitial_new_coord,boundary);% final position of the migration atom
        % the new target position
        direction_distance=new_interstitial_old_coord-new_interstitial_new_coord;
        direction_distance= dis_boundary_check(direction_distance, boundary);
        interstitial_new_I2_coord= new_interstitial_new_coord + 2*direction_distance;
        interstitial_new_target_coord= boundary_check(interstitial_new_I2_coord,boundary);% final position of the target atom (to form the new interstitial)

end

