%%  initalization
% clc;
clear;
addpath(genpath(pwd));
rng('default')
%% load the ML model
% load ('model_train_7nn.mat','model_set');
% load NiFe_final_10nn1nodes3.5RDT0.35d10.5d2.mat model_set_EBF
load NiFe_EBF_final_10nn11nodes3.5RDT0.35d10.5d2.mat model_set_EBF

ML_model= model_set_EBF;
% load ('model_train_dif_E_7nn.mat','model_set');
% load ('EBF_model.mat','model_set_EBF');
% load ('delta_E_model.mat','model_set');
%  load NiFe_dif_E_final_10nn3nodes3.5RDT0.35d10.5d2.mat model_set_FE
load NiFe_dif_E_final_10nn5nodes3.5RDT0.35d10.5d2.mat model_set_FE
Delta_E_model=model_set_FE;
%% load initial perfect structure
load([pwd,'/load_Data/per55']);
lattice_cos_ini= per55(5,3);
per_index=per55(:,1);
per_type=per55(:,2);
per_coords=per55(:,3:5);

fix_point=1e10;
frac_per=round((per_coords./lattice_cos_ini).*fix_point)./fix_point;
%% initial parameters setting
NN1=12;%fcc 1st nn
ratio_temp=(0:0.1:1)';
ratio=[ratio_temp,1-ratio_temp];
lattice_constant_set_Ni =[3.4989,3.4919, 3.4888,3.4872,...
    3.4874,3.4882,3.4923,3.4945,3.4996,...
    3.5077,3.5195];
super_size_length=10;
set_dis=1.0704;
Kb= 0.025852/300 ;% ev
% D0=2.5e12;%HZ
% D0=11e12;%HZ
%% the setting of attemp frequency
% if the attemp frequency is a set for varies components, for exampe
% [13, 13.1, 13.2].*1e12 corresponds to FexNi, x=[0.1 0.2 0.3], Then, you can uncommon
% the line of D0_set and commom the line of D0
%
D0=13.54e12;%500 K MD, constant
% D0_set = [your setteing attemp frequency];
%%
% for T = 500:200:1100
T=900;
L_size=60e4;
% for num_ratio=4
parfor num_ratio=1:length(ratio)
    %% D0_set
    % D0 = D0_set(num_ratio);
    %% update the ratio and lattice constant
    current_ratio=ratio(num_ratio,:);
    lattice_constant=lattice_constant_set_Ni(num_ratio);

    %%   creating interstitial structure
    perfect_structure_coords=frac_per.*lattice_constant;
    element_type_set=1:length(current_ratio);
    % the random type according to the ratio
    per_type= random_ratio_type(length(frac_per),element_type_set,current_ratio);
    %%test display
    index_type_Ni=find(per_type==1);
    index_type_Fe=find(per_type==2);
    disp('****************Ratio test*************************')
    disp(['Ni ratios:', num2str(length(index_type_Ni)/length(frac_per))]);
    disp(['Fe ratios:', num2str(length(index_type_Fe)/length(frac_per))]);
    %% generate the interstitial configuration
    perfect_structure=[per_index,per_type',perfect_structure_coords];
    % insert_pos_atom_ID=1;% boundary test
    insert_pos_atom_ID=2221;
    if rand<=current_ratio(1)
        insert_interstitial_type=1;%
    else
        insert_interstitial_type=2;% pure Fe
    end
    %     if num_ratio==1
    %         insert_interstitial_type=2;% pure Fe
    %     else
    %         insert_interstitial_type=1;% insert a Ni atom as the initial dumbbell,since Ni-Ni or Fe-Ni are more stable than Fe-Fe
    %     end
    L=lattice_constant*super_size_length;
    boundary = [0,L;0,L;0,L];

    %     disp(['The position before inserting [010] dumbbell:  ',num2str( perfect_structure(insert_pos_atom_ID,3:5))])
    % the final atom inform is added at the end of the initial perfect structure
    interstitial_structure=creat_dumbell_interstitial(perfect_structure,...
        insert_pos_atom_ID, insert_interstitial_type,set_dis, boundary);

    %     disp(['The positio of interstitial 1 :  ',num2str( interstitial_structure(insert_pos_atom_ID,3:5))])
    %
    %     disp(['The positio of interstitial 2 :  ',num2str( interstitial_structure(length(frac_per)+1,3:5))])

    % end
    % for num_ratio=1
    %% the standard coordinates
    interstitial_per_index=interstitial_structure(:,1);
    interstitial_per_types=interstitial_structure(:,2);
    initerstitial_per_coords=interstitial_structure(:,3:5);
    initia_I1_ID_back=interstitial_structure(end,1);
    initia_I2_ID_mig=insert_pos_atom_ID;% the insert atom set as the last index in interstitial_structure
    final_mig_I3_ID=2222;% this ID must the 1st NN of the insert_pos_atom
    standard_transformed_coords =coord_transform(initerstitial_per_coords,initia_I1_ID_back,initia_I2_ID_mig,...
        final_mig_I3_ID,lattice_constant,super_size_length);

    %% the type_index_cell is for distinct the different type atom
    type_index_cell=cell(1,length(element_type_set));
    for num_type=1:length(element_type_set)
        curent_index= find(interstitial_per_types==element_type_set(num_type));
        type_index_cell{num_type}=curent_index;
    end

    %% caluculation the nearest neibours, obtain the sort ID
    % this ID sort must satisfy the same rule of the  coord_transform
    % function, that has a same sequence as the ML train
    dif_vector =initerstitial_per_coords(initia_I1_ID_back,:)-initerstitial_per_coords(initia_I2_ID_mig,:);
    dif_vector =dis_boundary_check(dif_vector,boundary);% consider difference of the boundary condition
    central_position_mig_former=initerstitial_per_coords(initia_I2_ID_mig,:)+0.5.*dif_vector;% the dumbbell central before migration
    central_position_mig_former=boundary_check(central_position_mig_former,boundary);
    central_position_mig_latter =initerstitial_per_coords(final_mig_I3_ID,:); % the new I1 and I3 will generate a new central position
    % the range of the local atomic enviroment
    NN_num1=10;%
    %     NN_num2=NN_num1+2;
    [central1_ID_nn_set,central1_ID_count_set,relative_central1_coords]...
        = update_interstitial_nn_kmc(central_position_mig_former,initerstitial_per_coords,lattice_constant,super_size_length,NN_num1);
    %     fraction1_coords=relative_central1_coords./lattice_constant;
    %     [central2_ID_nn_set,central2_ID_count_set,relative_central2_coords]...
    %         = update_interstitial_nn_kmc(central_position_mig_latter,initerstitial_per_coords,lattice_constant,super_size_length,NN_num2);
    %     sort_ID=unique([central1_ID_nn_set;central2_ID_nn_set]);
    sort_ID=central1_ID_nn_set;
    %% input relative coordinates
    input_standard_coords=standard_transformed_coords(sort_ID,:);

    %% kMC prameters initialization
    t=0;
    R0=initerstitial_per_coords;
    Rt=R0;
    cross_boundary_record=zeros(size(R0));
    tracks= zeros(L_size,12);%[t,I1_ID,I1_x,I1_y,I1_z,I2_ID,I2_x,I2_y,I2_z,MSD,MSD_Fe,MSD_Ni], I1 is the mig interstitial, I2 is the back interstitial
    mig_index_set = zeros(L_size,1);
    mig_type_set = zeros(L_size,1);
    k_tot_avg=zeros(L_size,1);
    NN1_mig_ID_count_set=zeros(L_size,8);
    NN1_mig_energy_count_set=zeros(L_size,8);
    NN1_EBF_count_set=zeros(L_size,8);
    NN1_delta_E_count_set=zeros(L_size,8);
    interstitial_ID_back=initia_I1_ID_back;
    interstitial_ID_mig=initia_I2_ID_mig;
    N=length(interstitial_per_types);

    %% kMC loop

    tic;
    for count=1:L_size

        if mod(count,1e4)==0
            disp(count)
        end
        %% calculating MSD
        Rt= initerstitial_per_coords;
        Dt =(Rt-R0);

        Dt=Dt+L.* cross_boundary_record;
        sum_Dt=sum(Dt.^2,2);%

        MSD =sum(sum_Dt)./N;
        %partitial MSD
        MSD_element=zeros(1,length(type_index_cell));
        for num_type=1:length(type_index_cell)
            current_index=type_index_cell{num_type};
            if ~isempty(current_index)
                sum_cur=sum_Dt(current_index);
                %                      r_index=find(Co_index==vac_ID);
                MSD_element(num_type)=sum(sum_cur)./length(current_index);
            else
                MSD_element(num_type)=0;
            end
        end

        interstitial_back_coord = initerstitial_per_coords(interstitial_ID_back,:);
        interstitial_mig_coord = initerstitial_per_coords(interstitial_ID_mig,:);
        tracks(count,:)=[t,interstitial_ID_back,interstitial_back_coord,...
            interstitial_ID_mig,interstitial_mig_coord,MSD,MSD_element];

        %% obtaining the migtration energy of the 1st NN of the two interstitial
        largest_nn=1;
        % I1 has 4 possible migration paths
        [central1_ID_nn_set,central1_ID_count_set,~] ...
            = update_interstitial_nn_kmc(interstitial_back_coord,initerstitial_per_coords,lattice_constant,super_size_length,largest_nn+1);% the nearest atom is itself
        NN_I1_ID=central1_ID_nn_set(2:central1_ID_count_set(2)+1);% 2:central1_ID_count_set(2)+1 is the total 1st nn atoms
        interstitial_back_energy=zeros(1,length(NN_I1_ID));
        energy_temp_back_set = zeros(length(NN_I1_ID),3);
        for num_mig =1:length(NN_I1_ID)
            initia_I1_ID_back= interstitial_ID_mig;
            initia_I2_ID_mig= interstitial_ID_back;
            final_mig_I3_ID = NN_I1_ID(num_mig);

            transformed_coords = coord_transform(initerstitial_per_coords,initia_I1_ID_back,initia_I2_ID_mig,...
                final_mig_I3_ID,lattice_constant,super_size_length);
            %determining the same central ID as the standard coordinates
            dif_vector =initerstitial_per_coords(initia_I1_ID_back,:)-initerstitial_per_coords(initia_I2_ID_mig,:);
            dif_vector =dis_boundary_check(dif_vector,boundary);
            central_position_mig_former=initerstitial_per_coords(initia_I2_ID_mig,:)+0.5.*dif_vector;% the dumbbell central before migration
            central_position_mig_former=boundary_check(central_position_mig_former,boundary);
            central_position_mig_latter =initerstitial_per_coords(final_mig_I3_ID,:); % the new I1 and I3 will generate a new central position
            % the NN input of the migration before and after
            [central1_ID_nn_set,~,~]...
                = update_interstitial_nn_kmc(central_position_mig_former,initerstitial_per_coords,lattice_constant,super_size_length,NN_num1);

            %             [central2_ID_nn_set,~,~]...
            %                 = update_interstitial_nn_kmc(central_position_mig_latter,initerstitial_per_coords,lattice_constant,super_size_length,NN_num2);
            %             current_ID_set=unique([central1_ID_nn_set;central2_ID_nn_set]);
            current_ID_set=central1_ID_nn_set;
            current_coords_set = transformed_coords(current_ID_set,:);
            % sort the ID
            [C,index_current,index_input] = intersect(current_coords_set,input_standard_coords,'rows');%C= current_coords_set(index_current,:) and C = input_standard_coords(index_input,:)
            [sort_value, sort_standart_input_index]= sort(index_input);
            % adjust the current_ID has a same order with the standard one
            sort_standard_index_current= index_current(sort_standart_input_index);
            %             current_sort_ID=current_ID_set(sort_standart_input_index);% current_coords corresponds to the current_ID_set
            current_sort_ID=current_ID_set( sort_standard_index_current);
            % predict the energy
            current_input = interstitial_per_types(current_sort_ID);
            % encoding and predition
            %             input= one_hot_encode2(current_input');
            input=current_input;
            cur_mig_energy = ML_model(input);
            cur_delta_E_energy=Delta_E_model(input);

            %% modify 20220823
            %              initia_I1_ID_back= interstitial_ID_mig;
            %             initia_I2_ID_mig= interstitial_ID_back;
            %             final_mig_I3_ID = NN_I1_ID(num_mig);
            Iterstitial_type1=interstitial_per_types(initia_I1_ID_back);
            Iterstitial_type2=interstitial_per_types(initia_I2_ID_mig);
            Iterstitial_type3=interstitial_per_types(final_mig_I3_ID);
            %             mig_energy = mig_energy_determine (cur_mig_energy,cur_delta_E_energy, Iterstitial_type1,Iterstitial_type2,Iterstitial_type3);
            %             if Iterstitial_type1==Iterstitial_type3
            %                 mig_energy=cur_mig_energy;
            %             else
            %                 mig_energy=cur_mig_energy+cur_delta_E_energy;
            %             end
            mig_energy=cur_mig_energy;
            interstitial_back_energy(num_mig)= mig_energy;

            %            interstitial_back_energy(num_mig)= ML_model(current_input)+ Delta_E_model(current_input);
            energy_temp_back_set(num_mig,:) = [mig_energy,cur_mig_energy,cur_delta_E_energy];
        end

        % I2 has another 4 possible migration paths

        [central2_ID_nn_set,central2_ID_count_set,relative_central2_coords] ...
            = update_interstitial_nn_kmc(interstitial_mig_coord,initerstitial_per_coords,lattice_constant,super_size_length,largest_nn+1);
        NN_I2_ID=central2_ID_nn_set(2:end);
        interstitial_mig_energy=zeros(1,length(NN_I2_ID));
        energy_temp_mig_set = zeros(length(NN_I2_ID),3);
        for num_mig =1:length(NN_I2_ID)
            initia_I1_ID_back= interstitial_ID_back;
            initia_I2_ID_mig= interstitial_ID_mig;
            final_mig_I3_ID = NN_I2_ID(num_mig);
            transformed_coords = coord_transform(initerstitial_per_coords,initia_I1_ID_back,initia_I2_ID_mig,...
                final_mig_I3_ID,lattice_constant,super_size_length);
            %sort ID
            dif_vector =initerstitial_per_coords(initia_I1_ID_back,:)-initerstitial_per_coords(initia_I2_ID_mig,:);
            dif_vector =dis_boundary_check(dif_vector,boundary);
            central_position_mig_former=initerstitial_per_coords(initia_I2_ID_mig,:)+0.5.*dif_vector;% the dumbbell central before migration
            central_position_mig_former=boundary_check(central_position_mig_former,boundary);

            central_position_mig_latter =initerstitial_per_coords(final_mig_I3_ID,:); % the new I1 and I3 will generate a new central position
            % the NN input of the migration before and after
            [central1_ID_nn_set,~,~]...
                = update_interstitial_nn_kmc(central_position_mig_former,initerstitial_per_coords,lattice_constant,super_size_length,NN_num1);

            %             [central2_ID_nn_set,~,~]...
            %                 = update_interstitial_nn_kmc(central_position_mig_latter,initerstitial_per_coords,lattice_constant,super_size_length,NN_num2);
            %             current_ID_set=unique([central1_ID_nn_set;central2_ID_nn_set]);
            current_ID_set=central1_ID_nn_set;
            current_coords_set = transformed_coords(current_ID_set,:);
            % sort the ID
            [C,index_current,index_input] = intersect(current_coords_set,input_standard_coords,'rows');%C= current_coords_set(index_current,:) and C = input_standard_coords(index_input,:)
            [sort_value, sort_standart_input_index]= sort(index_input);
            % adjust the current_ID has a same order with the standard one
            sort_standard_index_current= index_current(sort_standart_input_index);
            %             current_sort_ID=current_ID_set(sort_standart_input_index);% current_coords corresponds to the current_ID_set
            current_sort_ID=current_ID_set( sort_standard_index_current);
            % predict the energy
            current_input = interstitial_per_types(current_sort_ID);
            % encoding and predition
            %             input= one_hot_encode2(current_input');
            input=current_input;
            cur_mig_energy = ML_model(input);
            cur_delta_E_energy=Delta_E_model(input);

            %% modify 20220823
            Iterstitial_type1=interstitial_per_types(initia_I1_ID_back);
            Iterstitial_type2=interstitial_per_types(initia_I2_ID_mig);
            Iterstitial_type3=interstitial_per_types(final_mig_I3_ID);
            %             mig_energy = mig_energy_determine (cur_mig_energy,cur_delta_E_energy, Iterstitial_type1,Iterstitial_type2,Iterstitial_type3);
            %             if Iterstitial_type1==Iterstitial_type3
            %                 mig_energy=cur_mig_energy;
            %             else
            %                 mig_energy=cur_mig_energy+cur_delta_E_energy;
            %             end
            mig_energy=cur_mig_energy;
            interstitial_mig_energy(num_mig)= mig_energy;
            %            interstitial_back_energy(num_mig)= ML_model(current_input)+ Delta_E_model(current_input);
            energy_temp_mig_set(num_mig,:) = [mig_energy,cur_mig_energy, cur_delta_E_energy];
        end
        NN1_mig_ID_set=[NN_I1_ID;NN_I2_ID]';
        NN1_mig_ID_count_set(count,:)=NN1_mig_ID_set;
        %%% modify the predict energy
        NN1_mig_energy_set=[interstitial_back_energy,interstitial_mig_energy];
        %index modify 20220614
        index_mo=find(NN1_mig_energy_set<=1e-2);
        %         NN1_mig_energy_set(index_mo)=0.01.*rand(1,length(index_mo));
        NN1_mig_energy_set(index_mo)=0.01;
        NN1_mig_energy_count_set(count,:)= NN1_mig_energy_set;
        NN1_EBF_count_set(count,:)=[energy_temp_back_set(:,2);energy_temp_mig_set(:,2)]';
        NN1_delta_E_count_set(count,:)=[energy_temp_back_set(:,3);energy_temp_mig_set(:,3)]';
        %             %             NN1_energy_set(count,:)= NN1_mig_energy_set;
        %% kMC
        K_set= D0.*exp(-(NN1_mig_energy_set)./(Kb*T));
        K_tot =sum(K_set);
        k_tot_avg(count)=K_tot;
        cum_k_set= cumsum(K_set);
        roulette_k_set = cum_k_set./K_tot;
        r1 =rand(1);
        mig_index = find(r1-roulette_k_set <0,1);
        r2 =rand(1);
        t = t + -1/K_tot* log(r2);

        %% after determining the migration atom, the initialsetting mig atom may became the back atom
        target_interstitial_ID = NN1_mig_ID_set(mig_index);
        mig_index_set(count)=mig_index;
        mig_type_set(count)=interstitial_per_types(NN1_mig_ID_set(mig_index));
        if mig_index<=length(NN1_mig_ID_set)/2
            temp_ID = interstitial_ID_mig;
            interstitial_ID_mig = interstitial_ID_back;
            interstitial_ID_back = temp_ID;
        end

        %update the back atom
        temp_coord =initerstitial_per_coords(interstitial_ID_mig,:);
        dif_int_coord = temp_coord-initerstitial_per_coords(interstitial_ID_back,:);
        dif_int_coord=dis_boundary_check(dif_int_coord,boundary);
        initerstitial_per_coords(interstitial_ID_back,:)=initerstitial_per_coords(interstitial_ID_back,:)+dif_int_coord/2;
        initerstitial_per_coords(interstitial_ID_back,:)= boundary_check( initerstitial_per_coords(interstitial_ID_back,:),boundary);

        % update the target atom and mig atom
        [interstitial_new_target_coord,interstitial_new_mig_coord] = ...
            find_new_dumbbell_position(initerstitial_per_coords,target_interstitial_ID,interstitial_ID_mig,set_dis,boundary);
        initerstitial_per_coords(target_interstitial_ID,:)=interstitial_new_target_coord;
        initerstitial_per_coords(interstitial_ID_mig,:)=interstitial_new_mig_coord;% uptate the migration atom

        interstitial_ID_back = target_interstitial_ID;% update the interstitial ID
        %% cross boundary statistic
        initial_pos=Rt;
        final_pos=initerstitial_per_coords;
        cross_boundary_record = cross_boundary_statistic(initial_pos,final_pos,boundary,cross_boundary_record);

    end
    %% plot
    t_set=tracks(:,1);
    MSD_t=tracks(:,10);
    MSD_Fe=tracks(:,11);
    MSD_Ni=tracks(:,12);
    p0= polyfit(t_set, MSD_t,1);
    p1= polyfit(t_set, MSD_Fe,1);
    p2= polyfit(t_set, MSD_Ni,1);
    x0=0:max(t_set)/1e3:max(t_set);
    y0= polyval(p0,x0);
    y1= polyval(p1,x0);
    y2= polyval(p2,x0);
    D_all=p0(1)./6;
    D_Fe=p1(1)./6;
    D_Ni=p2(1)./6;
    Ni_num=length(find(interstitial_per_types==1));
    Fe_num=length(find(interstitial_per_types==2));
    Ni_mig_num= length(find(mig_type_set==1));
    Fe_mig_num= length(find(mig_type_set==2));
    unit_a=lattice_constant^2/2;
    f_tracer_all= MSD_t(end)/((length(MSD_t)/N) * unit_a);
    if  ~isempty(Fe_mig_num)

        f_tracer_Fe= MSD_Fe(end)/( Fe_mig_num/Fe_num * unit_a);
    else
        f_tracer_Fe=0;
    end
    if  ~isempty(Ni_mig_num)

        f_tracer_Ni= MSD_Ni(end)/(Ni_mig_num/Ni_num * unit_a);
    else
        f_tracer_Ni=0;
    end
    %     jump frequency
    mean_JF = mean(1./k_tot_avg);
    JF_all = 1/mean_JF;
    Mig_prob_Ni=Ni_mig_num/length(mig_type_set);
    Mig_prob_Fe=Fe_mig_num/length(mig_type_set);
    JF_Ni = Mig_prob_Ni./current_ratio(:,1)*JF_all;
    JF_Fe = Mig_prob_Fe./current_ratio(:,2)*JF_all;

    %% store
    interstial_struct_cell{num_ratio}=[interstitial_per_index,interstitial_per_types,initerstitial_per_coords];
    R0_cell{num_ratio}= R0;
    Rt_cell{num_ratio}= Rt;
    cross_boundary_cell{num_ratio}= cross_boundary_record;
    tracks_cell{num_ratio}=tracks;
    k_tot_avg_cell{num_ratio}= k_tot_avg;
    mig_index_cell{num_ratio}=mig_index_set;
    NN1_mig_ID_count_cell{num_ratio}=NN1_mig_ID_count_set;
    NN1_mig_energy_cell{num_ratio}=NN1_mig_energy_count_set;
    NN1_EBF_count_cell{num_ratio}= NN1_EBF_count_set;
    NN1_delta_E_count_cell{num_ratio}=NN1_delta_E_count_set;
    D_cell{num_ratio}=[D_all,D_Fe,D_Ni];
    f_cell{num_ratio}=[f_tracer_all,f_tracer_Ni,f_tracer_Fe];
    JF_cell{num_ratio}=[JF_all,JF_Ni,JF_Fe];
    fit_parameter_P0{num_ratio}=p0;
    fit_parameter_P1{num_ratio}=p1;
    fit_parameter_P2{num_ratio}=p2;

end
% datestr('now','HHMMSS')
save (['Interstitial_all_ratio_',num2str(T),'K.mat']);
% end

