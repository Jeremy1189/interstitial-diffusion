% initialization
clc;
clear;
addpath(genpath(pwd));
close all;
rng('default');

%% initialization data
%% load initial perfect structure,% initial parameters setting
load([pwd,'/load_Data/per55']);
lattice_cos_ini= per55(5,3);
per_coords=per55(:,3:5);
fix_point=1e8;
kMC_init_para.frac_coord_supercell=round((per_coords./lattice_cos_ini).*fix_point)./fix_point;
kMC_init_para.v0=13.54e12;%HZ
kMC_init_para.T=1100;
kMC_init_para.lattice_constant=3.5;
kMC_init_para.L_size=5e4;
kMC_init_para.super_size_length=10;
kMC_init_para.dumbell_interval_dis=1.0704;
kMC_init_para.Kb= 0.025852/300 ;% ev
%boundary
L = kMC_init_para.lattice_constant*kMC_init_para.super_size_length;
kMC_init_para.boundary = [0,L;0,L;0,L];
%perfect coordinates
kMC_init_para.per_coords = kMC_init_para.frac_coord_supercell.*kMC_init_para.lattice_constant;
% stored interstitial positions and corresponding NN1 positions
[ ~,kMC_init_para.int_pos_set, kMC_init_para.NN1_ID_pos_cells] =...
    cal_distance_matrix_int(kMC_init_para.per_coords,kMC_init_para.boundary, kMC_init_para.dumbell_interval_dis);

avg_barrier_set = [0.286652866000000,0.484983508054523,0.231786671902269,0.445508313146233,0.0631189392405063,0.150463577981651,0.0650107817638266,0.153396858741259];
avg_barrier_min_set=[0.499854078066214	0.650000000000000	0.396614282656689	0.497629984229275	0.0100000000000000	0.400000000000000	0.0100000000000000	0.384581296967475];
avg_barrier_max_set=[0.320000000000000	0.350000000000000	0.100000000000000	0.209426702661019	0.150237896349729	0.165264334435154	0.0371069399379546	0.132588068840425];
for num_times=1
for ratio=0.5
kMC_init_para.ratio=[ratio, 1-ratio];
element_type_set =cumsum(ones(1,length(kMC_init_para.ratio)));
kMC_init_para.per_type= random_ratio_type(length(per_coords)+1,element_type_set,ratio);
%% intialization of differential evolution algorithm

% NP = 35;% number of atoms in population
% delete(gcp('nocreate'))
% parpool(36)
% There are 8 dimensions with the lu, which corresponds the dumbbell types as follows
% types_set=[ Iterstitial_type1,Iterstitial_type2,Iterstitial_type3];
% matrix_set=[1,1,1;            
%             1,1,2;
%             1,2,1;            
%             1,2,2;
%             2,1,1;            
%             2,1,2;
%             2,2,1;
%             2,2,2;];
% lu=[zeros(1,8)+0.01;ones(1,8).*1];% this means all the migration energy of each dumbbell types are in the range [0,1]
% lu = [0.2, 0.35, 0.1, 0.2, 0.01, 0.1, 0.01, 0.1;...
%       0.5, 0.65, 0.4, 0.5, 0.31, 0.4, 0.31, 0.4];
% initialization the population
pop=[avg_barrier_set;avg_barrier_min_set;avg_barrier_max_set];

% initial fitness value  
[fitpop,tracks_cell, type_set_cell]  = fitness_D(repmat(pop,4,1),kMC_init_para);
disp(fitpop)

end
save([pwd,'/save_data/A',num2str(kMC_init_para.ratio(1)),'B_',num2str(num_times),'1.mat'])

end
