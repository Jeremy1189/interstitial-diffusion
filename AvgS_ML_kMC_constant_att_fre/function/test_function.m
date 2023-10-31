% test interstitial_coords_trans
clc;
clear;

addpath(genpath(pwd));
%% test creat_dumbell_interstitial
load([pwd,'/load_Data/per55']);
perfect_structure=per55;
insert_pos_atom_ID=2221;
% insert_pos_atom_ID=1;% boundary test
insert_interstitial_type=perfect_structure(insert_pos_atom_ID,2);
set_dis = 1.25;
lattice_constant=perfect_structure(5,3);
supersize=10;
L=lattice_constant*supersize;
boundary = [0,L;0,L;0,L];
interstitial_structure=creat_dumbell_interstitial(perfect_structure,...
 insert_pos_atom_ID, insert_interstitial_type,set_dis, boundary);
%% test coord_transform
initia_I1_ID_back=2221;
initia_I2_ID_mig=4001;
final_mig_I3_ID=2222;
transformed_coords =coord_transform(interstitial_structure,initia_I1_ID_back,initia_I2_ID_mig,...
    final_mig_I3_ID,set_dis,lattice_constant,supersize);


