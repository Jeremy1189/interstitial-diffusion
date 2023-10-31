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
kMC_init_para.L_size=1.5e4;
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
delete(gcp('nocreate'))
parpool(20)
for num_times=1:10
for ratio=0.05
kMC_init_para.ratio=[ratio, 1-ratio];
element_type_set =cumsum(ones(1,length(kMC_init_para.ratio)));
kMC_init_para.per_type= random_ratio_type(length(per_coords)+1,element_type_set,ratio);
%% intialization of differential evolution algorithm
NP = 20;% number of atoms in population
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
switch ratio
    case 0.5
        lu = [0.2, 0.35, 0.1, 0.2, 0.01, 0.1, 0.01, 0.1;...
            0.5, 0.65, 0.4, 0.5, 0.31, 0.4, 0.31, 0.4];
    case 0.2
        lu = [0.2, 0.35, 0.1, 0.2, 0.01, 0.1, 0.01, 0.1;...
            0.5, 0.65, 0.4, 0.5, 0.31, 0.4, 0.31, 0.4];
    case 0.8
        lu = [0.2, 0.45, 0.1, 0.4, 0.01, 0.1, 0.01, 0.1;...
            0.5, 0.75, 0.4, 0.7, 0.31, 0.4, 0.31, 0.4];
    case 0.05
        lu = [0.2, 0.35, 0.1, 0.2, 0.01, 0.1, 0.01, 0.1;...
            0.5, 0.65, 0.4, 0.5, 0.31, 0.4, 0.31, 0.4];
end

% initialization the population
n = size(lu, 2);
pop = ones(NP, 1) * lu(1, :) + rand(NP, n) .* (ones(NP, 1) * (lu(2, :) - lu(1, :)));
%% modify random pop
pop=DE_boundary_check(pop,lu);
constrain_violate = constrain_handle_check(pop);
% initial fitness value 
violate_scale=100;
fitpop = fitness_D(pop,kMC_init_para).*(1+violate_scale.*sum(constrain_violate,2));

% contral parameters
F = 0.9;% large benifit for exploring
CR = 0.9;% large for more cross
% main loop
count_updated = 0;
results = [];
best_pop_set = [];
best_fitness=min(fitpop);
stop_condition_times =50;
while(1)
    % Implement DE for each subpopulation
    offSubpop = DE(pop, NP, lu, n, F, CR);
%     offSubpop = constrain_check(offSubpop,lu);
    offSubpop=DE_boundary_check(offSubpop,lu);
     constrain_violate = constrain_handle_check(offSubpop);
    % Evaluate the offspring subpopulation based on the first
    % fitness function
    fitOffSubpop = fitness_D(offSubpop,kMC_init_para).*(1+violate_scale.*sum(constrain_violate,2));
    
    % Implement selection between the parent subpopulation and
    % the offspring subpopulation
    new_offspring_index = find(fitpop >= fitOffSubpop);
    pop(new_offspring_index, :) = offSubpop(new_offspring_index, :);
    fitpop(new_offspring_index, :) = fitOffSubpop(new_offspring_index, :);
    [cur_best_fitness,cur_best_index]=min(fitpop);
    if cur_best_fitness< best_fitness
       count_updated=0;
       best_fitness =  cur_best_fitness;
       results = [results, best_fitness];
       best_pop_set = [best_pop_set; pop(cur_best_index,:)];
    else 
        count_updated=count_updated+1;
    end
    if count_updated>stop_condition_times
            break;
    end
end
save([pwd,'/save_data/A',num2str(kMC_init_para.ratio(1)),'B_',num2str(num_times),'_50.mat'])
end
end
