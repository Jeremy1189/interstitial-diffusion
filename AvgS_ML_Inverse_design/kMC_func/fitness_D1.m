function fitness_value = fitness_D1(pop,kMC_init_para)

 tic;
 avg_barrier =pop;
 [R, ~] = size(pop);
 D_set=zeros(R,1);

%% initial parameters setting
v0=kMC_init_para.v0;%HZ
T=kMC_init_para.T;
L_size=kMC_init_para.L_size;
dumbell_interval_dis=kMC_init_para.dumbell_interval_dis;
boundary = kMC_init_para.boundary;
per_coords = kMC_init_para.per_coords; 
interstitial_per_types=kMC_init_para.per_type;
int_pos_set=kMC_init_para.int_pos_set;
NN1_ID_pos_cells=kMC_init_para.NN1_ID_pos_cells;

%% 
int_ID_mig = 2221;
int_ID_back=4001;
% parfor i=1:R
 for i=1:R

% ratio_cur = ratio;
mig_energy = avg_barrier(i,:);
% coordinates of create intestitial
initerstitial_coords=[per_coords;per_coords(int_ID_mig,:)+[-dumbell_interval_dis,0,0]];
initerstitial_coords(int_ID_mig,:) = initerstitial_coords(int_ID_mig,:)+[dumbell_interval_dis,0,0];


%% kMC initialization
t=0;
R0=initerstitial_coords;
Rt=R0;
cross_boundary_record=zeros(size(R0));
% tracks dimensions: [t,I1_ID,I1_x,I1_y,I1_z,I2_ID,I2_x,I2_y,I2_z,MSD,MSD_Fe,MSD_Ni],
%     I1 is the mig interstitial, I2 is the back interstitial
tracks= zeros(L_size,4);
interstitial_ID_back=int_ID_back;
interstitial_ID_mig=int_ID_mig;
    
%% main loop
for count = 1: L_size
% for count = 10000    
    % calculate the MSD
%     tic;
    [SD, SD_element] = cal_SD(Rt,R0, cross_boundary_record, interstitial_per_types, boundary);
%     disp('cal_SD')
%     toc;
    % record tracks
    tracks(count,:)=[t, SD, SD_element];
    
    % dumbbell positions
%     interstitial_back_coord = initerstitial_coords(interstitial_ID_back,:);
%     interstitial_mig_coord = initerstitial_coords(interstitial_ID_mig,:);
%    
    %% determining the 1st NN IDs of dumbell 
    % each atom of dumbbell has 4 possible migration paths
    % back atom and its 4 nearest neighbour atoms
%       tic;
%     [central1_ID_nn_set,central1_ID_count_set,~] ...
%         = update_interstitial_nn_kmc(interstitial_back_coord,initerstitial_coords,lattice_constant,super_size_length,largest_nn+1);% the nearest atom is itself
%     NN_I1_ID=central1_ID_nn_set(2:central1_ID_count_set(2)+1);% 2:central1_ID_count_set(2)+1 is the total 1st nn atoms
%    
%     % migration atom and its 4 nearest neighbour atoms
%     [central2_ID_nn_set,central2_ID_count_set,~] ...
%         = update_interstitial_nn_kmc(interstitial_mig_coord,initerstitial_coords,lattice_constant,super_size_length,largest_nn+1);
%     NN_I2_ID=central2_ID_nn_set(2:central2_ID_count_set(2)+1);
    NN1_back_IDs = determine_NN1_IDs(initerstitial_coords,interstitial_ID_back,int_pos_set, NN1_ID_pos_cells);
    NN1_mig_IDs = determine_NN1_IDs(initerstitial_coords,interstitial_ID_mig,int_pos_set, NN1_ID_pos_cells);
%total mig_IDs
    NN1_mig_ID_set=[NN1_back_IDs;NN1_mig_IDs]';
%         disp('update_interstitial_nn_kmc')
%     toc;
    %% migration energy for each possible paths

    final_int_types1=interstitial_per_types(NN1_back_IDs);
    final_int_types2=interstitial_per_types(NN1_mig_IDs);
    % 8 migration paths
    initial_int_types1=interstitial_per_types(interstitial_ID_back)';
    initial_int_types2=interstitial_per_types(interstitial_ID_mig)';
    NN1_mig_energy_set=zeros(1,length(NN1_mig_ID_set));
    for num_IDs=1:length(NN1_mig_ID_set)
        if num_IDs<=length(NN1_back_IDs)
            NN1_mig_energy_set(num_IDs)= mig_energy_determine (mig_energy,initial_int_types2,initial_int_types1,final_int_types1(num_IDs));
        elseif  length(NN1_back_IDs )<num_IDs && num_IDs<=length(NN1_mig_ID_set)
            NN1_mig_energy_set(num_IDs)= mig_energy_determine (mig_energy,initial_int_types1,initial_int_types2,final_int_types2(num_IDs-4));
        else
            disp('Error types in FCC crystal!!!');
        end
    end
    %% kMC
%       tic;
    [t,mig_index, ~] = KMC(NN1_mig_energy_set,v0,T,t);
%             disp('kMC')
%     toc;
    %% after determining the migration atom, the initial setting mig atom may became the back atom
    target_interstitial_ID = NN1_mig_ID_set(mig_index);
    if mig_index<=length(NN1_mig_ID_set)/2
        temp_ID = interstitial_ID_mig;
        interstitial_ID_mig = interstitial_ID_back;
        interstitial_ID_back = temp_ID;
    end
%       tic;
    %update the back atom
    temp_coord =initerstitial_coords(interstitial_ID_mig,:);
    dif_int_coord = temp_coord-initerstitial_coords(interstitial_ID_back,:);
    dif_int_coord=dis_boundary_check(dif_int_coord,boundary);
    initerstitial_coords(interstitial_ID_back,:)=initerstitial_coords(interstitial_ID_back,:)+dif_int_coord/2;
    initerstitial_coords(interstitial_ID_back,:)= boundary_check( initerstitial_coords(interstitial_ID_back,:),boundary);
    
    % update the target atom and mig atom
    [interstitial_new_target_coord,interstitial_new_mig_coord] = ...
        find_new_dumbbell_position(initerstitial_coords,target_interstitial_ID,interstitial_ID_mig,dumbell_interval_dis,boundary);
    initerstitial_coords(target_interstitial_ID,:)=interstitial_new_target_coord;
    initerstitial_coords(interstitial_ID_mig,:)=interstitial_new_mig_coord;% uptate the migration atom
    
    interstitial_ID_back = target_interstitial_ID;% update the interstitial ID
%             disp('update_interstitial')
%     toc;
    %% cross boundary statistic
    initial_pos=Rt;
    final_pos=initerstitial_coords;
    cross_boundary_record = cross_boundary_statistic(initial_pos,final_pos,boundary,cross_boundary_record);
    %% updated current structure 
    Rt =final_pos;
end
%% calculate the diffusion coefficeint
D_set(i) = cal_diffusion_coef(tracks);
end
fitness_value = D_set(:,1);
toc;


