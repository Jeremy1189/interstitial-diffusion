% initialization
clc;
clear;
addpath(genpath(pwd));
close all;
rng('default');
%% load initial perfect structure
load([pwd,'/load_Data/per55']);
lattice_cos_ini= per55(5,3);
per_index=per55(:,1);
per_type=per55(:,2);
per_coords=per55(:,3:5);

fix_point=1e10;
frac_coord_supercell=round((per_coords./lattice_cos_ini).*fix_point)./fix_point;
%% initial parameters setting
NN1=12;%fcc 1st nn
ratio_temp=(0:0.1:1)';
ratio=[ratio_temp,1-ratio_temp];
lattice_constant_set_Ni =[3.4989,3.4919, 3.4888,3.4872,...
    3.4874,3.4882,3.4923,0.4945,3.4996,...
    3.5077,3.5195];
lat_cons_fit = polyfit(ratio_temp',lattice_constant_set_Ni,4);
% plot fitting result
% figure;
% plot(ratio_temp,lattice_constant_set_Ni,'or');
% hold on;
% x=0:1e-2:1;
% plot(x,polyval(lat_cons_fit,x),'--')
% xlabel('Ni_xFe{1-x}')
% ylabel(['lattice constant (',char(197),')'])

super_size_length=10;
dumbell_interval_dis=1.0704;
Kb= 0.025852/300 ;% ev
% D0=11.83e12;%HZ
v0=13.54e12;%HZ
T=1100;
% T_set=[500,600,700,900,1000,1100];
L_size=1e3;
% num_ratio=11;
%% energy
path_avg=pwd;
load([path_avg,'/../../NiFe_data/data_handle_code/avg_EBF_set.mat'],'avg_barrier_set');
pure_Ni_mig_energy=0.3414;
pure_Fe_mig_energy=0.2414;
avg_barrier=[ones(1,8)*pure_Fe_mig_energy;avg_barrier_set;ones(1,8)*pure_Ni_mig_energy];


%% create interstitial strucutre
tic;
n=5;
ratio_cur = ratio(n,:);
mig_energy = avg_barrier(n,:);
lattice_constant=polyval(lat_cons_fit,ratio_cur(1));
L = lattice_constant*super_size_length;
boundary = [0,L;0,L;0,L];
[interstitial_structure,int_ID_back,int_ID_mig] = ...
    create_intestitial(ratio_cur,frac_coord_supercell,lattice_constant,boundary,dumbell_interval_dis);
interstitial_per_index=interstitial_structure(:,1);
interstitial_per_types=interstitial_structure(:,2);
initerstitial_per_coords=interstitial_structure(:,3:5);
%% kMC initialization
t=0;
R0=initerstitial_per_coords;
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
    interstitial_back_coord = initerstitial_per_coords(interstitial_ID_back,:);
    interstitial_mig_coord = initerstitial_per_coords(interstitial_ID_mig,:);
    largest_nn=1;
    %% determining the 1st NN IDs of dumbell 
    % each atom of dumbbell has 4 possible migration paths
    % back atom and its 4 nearest neighbour atoms
%       tic;
    [central1_ID_nn_set,central1_ID_count_set,~] ...
        = update_interstitial_nn_kmc(interstitial_back_coord,initerstitial_per_coords,lattice_constant,super_size_length,largest_nn+1);% the nearest atom is itself
    NN_I1_ID=central1_ID_nn_set(2:central1_ID_count_set(2)+1);% 2:central1_ID_count_set(2)+1 is the total 1st nn atoms
   
    % migration atom and its 4 nearest neighbour atoms
    [central2_ID_nn_set,central2_ID_count_set,~] ...
        = update_interstitial_nn_kmc(interstitial_mig_coord,initerstitial_per_coords,lattice_constant,super_size_length,largest_nn+1);
    NN_I2_ID=central2_ID_nn_set(2:central2_ID_count_set(2)+1);
    %total mig_IDs
    NN1_mig_ID_set=[NN_I1_ID;NN_I2_ID]';
%         disp('update_interstitial_nn_kmc')
%     toc;
    %% migration energy for each possible paths

    final_int_types1=interstitial_per_types(NN_I1_ID);
    final_int_types2=interstitial_per_types(NN_I2_ID);
    % 8 migration paths
    initial_int_types1=interstitial_per_types(interstitial_ID_back)';
    initial_int_types2=interstitial_per_types(interstitial_ID_mig)';
    NN1_mig_energy_set=zeros(1,length(NN1_mig_ID_set));
    for num_IDs=1:length(NN1_mig_ID_set)
        if num_IDs<=length(NN_I1_ID)
            NN1_mig_energy_set(num_IDs)= mig_energy_determine (mig_energy,initial_int_types2,initial_int_types1,final_int_types1(num_IDs));
        elseif  length(NN_I1_ID)<num_IDs && num_IDs<=length(NN1_mig_ID_set)
            NN1_mig_energy_set(num_IDs)= mig_energy_determine (mig_energy,initial_int_types1,initial_int_types2,final_int_types2(num_IDs-4));
        else
            disp('Error types in FCC crystal!!!');
        end
    end
    %% kMC
%       tic;
    [t,mig_index, K_tot] = KMC(NN1_mig_energy_set,v0,T,t);
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
    temp_coord =initerstitial_per_coords(interstitial_ID_mig,:);
    dif_int_coord = temp_coord-initerstitial_per_coords(interstitial_ID_back,:);
    dif_int_coord=dis_boundary_check(dif_int_coord,boundary);
    initerstitial_per_coords(interstitial_ID_back,:)=initerstitial_per_coords(interstitial_ID_back,:)+dif_int_coord/2;
    initerstitial_per_coords(interstitial_ID_back,:)= boundary_check( initerstitial_per_coords(interstitial_ID_back,:),boundary);
    
    % update the target atom and mig atom
    [interstitial_new_target_coord,interstitial_new_mig_coord] = ...
        find_new_dumbbell_position(initerstitial_per_coords,target_interstitial_ID,interstitial_ID_mig,dumbell_interval_dis,boundary);
    initerstitial_per_coords(target_interstitial_ID,:)=interstitial_new_target_coord;
    initerstitial_per_coords(interstitial_ID_mig,:)=interstitial_new_mig_coord;% uptate the migration atom
    
    interstitial_ID_back = target_interstitial_ID;% update the interstitial ID
%             disp('update_interstitial')
%     toc;
    %% cross boundary statistic
    initial_pos=Rt;
    final_pos=initerstitial_per_coords;
    cross_boundary_record = cross_boundary_statistic(initial_pos,final_pos,boundary,cross_boundary_record);
    %% updated current structure 
    Rt =final_pos;
end
%% calculate the diffusion coefficeint
D_set = cal_diffusion_coef(tracks);
disp(D_set);
toc;


