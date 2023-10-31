function  [NN1_structure_cell,int_pos_set, NN1_ID_pos_cells] = cal_distance_matrix_int(per_coords,boundary, dumbell_interval_dis)

num_atoms = size(per_coords,1);
%% inset position of intestitial
fix_point=1e8;
largest_nn=2;
NN1_structure_cell = cell(num_atoms,1);
int_pos_set = zeros(num_atoms*6,3); % each position has 2 atoms in dumbbell with [100],[010],[001]
NN1_ID_pos_cells = cell(num_atoms*6,1);
count=0;
for i =1: num_atoms
    
    interstitial_centre = per_coords(i,:);
    perfect_structure = per_coords;
    % [100] dumbell
    interstitial_I1_100 = interstitial_centre + [dumbell_interval_dis,0,0];
    interstitial_I1_100 = boundary_check(interstitial_I1_100,boundary);
    perfect_structure(i,:) = interstitial_I1_100;
    interstitial_I2_100 = interstitial_centre + [-dumbell_interval_dis,0,0];
    interstitial_I2_100 = boundary_check(interstitial_I2_100,boundary);
    interstitial_structure_100 =[perfect_structure(:,:);interstitial_I2_100];
    [I1_IDs,~,~] = update_interstitial_nn_kmc(interstitial_I1_100,interstitial_structure_100,boundary,largest_nn);
    I1_IDs=I1_IDs(2:end);% the first value is itself
    interstitial.NN1_pos_I1_100= interstitial_structure_100(I1_IDs,:);
    [I2_IDs,~,~] = update_interstitial_nn_kmc(interstitial_I2_100,interstitial_structure_100,boundary,largest_nn);
    I2_IDs=I2_IDs(2:end);
    interstitial.NN1_pos_I2_100= interstitial_structure_100(I2_IDs,:);
    interstitial.interstitial_I1_100 = interstitial_I1_100;
    interstitial.interstitial_I2_100 = interstitial_I2_100;
    count=count+1;
    int_pos_set(count,:) =round(interstitial_I1_100.*fix_point)./fix_point; 
    NN1_ID_pos_cells{count} = round(interstitial.NN1_pos_I1_100.*fix_point)./fix_point;
    count=count+1;
    int_pos_set(count,:) = round(interstitial_I2_100.*fix_point)./fix_point;
    NN1_ID_pos_cells{count} = round(interstitial.NN1_pos_I2_100.*fix_point)./fix_point;
    
    % [010] dumbell
    interstitial_I1_010 = interstitial_centre + [0,dumbell_interval_dis,0];
    interstitial_I1_010 = boundary_check(interstitial_I1_010,boundary);
    perfect_structure(i,:) = interstitial_I1_010;
    interstitial_I2_010 = interstitial_centre + [0,-dumbell_interval_dis,0];
    interstitial_I2_010= boundary_check(interstitial_I2_010,boundary);
    interstitial_structure_010 =[perfect_structure(:,:);interstitial_I2_010];
    
    [I1_IDs,~,~] = update_interstitial_nn_kmc(interstitial_I1_010,interstitial_structure_010,boundary,largest_nn);
    I1_IDs=I1_IDs(2:end);
    interstitial.NN1_pos_I1_010= interstitial_structure_100(I1_IDs,:);
    [I2_IDs,~,~] = update_interstitial_nn_kmc(interstitial_I2_010,interstitial_structure_010,boundary,largest_nn);
    I2_IDs=I2_IDs(2:end);
    interstitial.NN1_pos_I2_010= interstitial_structure_010(I2_IDs,:);
    interstitial.interstitial_I1_010 = interstitial_I1_010;
    interstitial.interstitial_I2_010 = interstitial_I2_010;
        count=count+1;
    count=count+1;
    int_pos_set(count,:) =round(interstitial_I1_010.*fix_point)./fix_point; 
    NN1_ID_pos_cells{count} = round(interstitial.NN1_pos_I1_010.*fix_point)./fix_point;
    count=count+1;
    int_pos_set(count,:) = round(interstitial_I2_010.*fix_point)./fix_point;
    NN1_ID_pos_cells{count} = round(interstitial.NN1_pos_I2_010.*fix_point)./fix_point;
    
    
    % [001] dumbell
    interstitial_I1_001 = interstitial_centre + [0,0,dumbell_interval_dis];
    interstitial_I1_001 = boundary_check(interstitial_I1_001,boundary);
    perfect_structure(i,:) = interstitial_I1_001;
    interstitial_I2_001 = interstitial_centre + [0,0,-dumbell_interval_dis];
    interstitial_I2_001 = boundary_check(interstitial_I2_001,boundary);
    interstitial_structure_001 =[perfect_structure(:,:);interstitial_I2_001];
    [I1_IDs,~,~] = update_interstitial_nn_kmc(interstitial_I1_001,interstitial_structure_001,boundary,largest_nn);
    I1_IDs=I1_IDs(2:end);
    interstitial.NN1_pos_I1_001= interstitial_structure_100(I1_IDs,:);
    [I2_IDs,~,~] = update_interstitial_nn_kmc(interstitial_I2_001,interstitial_structure_001,boundary,largest_nn);
    I2_IDs=I2_IDs(2:end);
    interstitial.NN1_pos_I2_001= interstitial_structure_001(I2_IDs,:);
    interstitial.interstitial_I1_001 = interstitial_I1_001;
    interstitial.interstitial_I2_001 = interstitial_I2_001;
    count=count+1;
    int_pos_set(count,:) =round(interstitial_I1_001.*fix_point)./fix_point; 
    NN1_ID_pos_cells{count} = round(interstitial.NN1_pos_I1_001.*fix_point)./fix_point;
    count=count+1;
    int_pos_set(count,:) = round(interstitial_I2_001.*fix_point)./fix_point;
    NN1_ID_pos_cells{count} = round(interstitial.NN1_pos_I2_001.*fix_point)./fix_point;
    
    %%
    NN1_structure_cell{i} = interstitial;

end
