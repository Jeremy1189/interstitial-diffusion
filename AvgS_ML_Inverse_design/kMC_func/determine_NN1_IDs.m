function NN1_IDs = determine_NN1_IDs(dumbbell_coords,Int_ID,int_pos_set, NN1_ID_pos_cells)
%DETERMINE_NN1_IDS Summary of this function goes here
%   Detailed explanation goes here
fix_point=1e8;
interstitial_pos = dumbbell_coords(Int_ID,:);
interstitial_pos=round(interstitial_pos.*fix_point)./fix_point;
[~,~,index]= intersect(interstitial_pos,int_pos_set,'rows');
NN1_positions = NN1_ID_pos_cells{index};
NN1_positions=round(NN1_positions.*fix_point)./fix_point;
[~,~,NN1_IDs] =intersect(NN1_positions,dumbbell_coords,'rows');

end

