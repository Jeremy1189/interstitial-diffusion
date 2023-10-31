function interstitial_structure=creat_dumbell_interstitial(perfect_structure,...
 insert_pos_atom_ID, insert_interstitial_type,set_dis, boundary)
% perfect_structure size : N*5, N is the number of atom
Interstitial_coord_temp=perfect_structure(insert_pos_atom_ID,3:5);
I1_coord=Interstitial_coord_temp+[0,set_dis,0]; %[010] dumbbell interstitial
I1_coord = boundary_check(I1_coord,boundary);
perfect_structure(insert_pos_atom_ID,3:5)=I1_coord;
I2_coord= Interstitial_coord_temp+[0,-set_dis,0];% [010] 
I2_coord = boundary_check(I2_coord,boundary);
I2_index=perfect_structure(end,1)+1;
interstitial_structure=[perfect_structure;...
      I2_index,insert_interstitial_type,I2_coord];
end