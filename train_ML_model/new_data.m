load EBF_set.mat result_set type_total_set
result_set0=[];
type_data=[];
[R,C]=size(result_set);
set=[1:8];
for i=1:R
    if i==8
        continue;
    end
    result_set0=[result_set0;cell2mat(result_set(i,set)')];
    type_data=[type_data,cell2mat(type_total_set(i,set))];
end
result_data=result_set0';
dim=5;
tmp_per=read_data_tmp([pwd,'/load_Data/tmp2'],dim);
% tmp_per=read_data_tmp([pwd,'/load_Data/per55'],dim);
tmp_ID=tmp_per(1,:);
[sort_v,sort_index]= sort(tmp_ID);
per=tmp_per(:,sort_index)';
%% Create Interstitial
% load per55; % perfect structure
I1_ID =2221;
I2_ID=4001;
I3_ID=2222;
set_dis=1.05;
lattice_constant=per(5,3);
supersize=10;
I1_coord=per(I1_ID,3:5);
I2_coord=per(I2_ID,3:5);
I3_coord=per(I3_ID,3:5);
%% create the initial interstitial structure
% per =[per55;...
%     [I2_ID,I2_type,I2_coord]];
per_coords=per(:,3:5);
per_ID_set = per(:,1);
%% caluculation the nearest neibours, obtain the sort ID
% central_position_mig_former= (I1_coord+I2_coord)./2;% the dumbbell central before migration
dif_vector =I1_coord-I2_coord;
L=lattice_constant*supersize;
boundary = [0,L;0,L;0,L];
dif_vector =dis_boundary_check(dif_vector,boundary);
central_position=I2_coord+0.5.*dif_vector;% the dumbbell central before migration
central_position=boundary_check(central_position,boundary);
NN_num1=10;
% NN_num2=NN_num1+2;
[central1_ID_nn_set,central1_ID_count_set,relative_central1_coords]...
    = update_interstitial_nn_kmc(central_position,per_coords,lattice_constant,supersize,NN_num1);

%% load the input and output data
perf_set=[];
I2_types = result_data(3,:);
type_data(I1_ID,:) = result_data(1,:);
type_data(I3_ID,:) = result_data(10,:);
types = [type_data;I2_types];
EBF = result_data(16,:);
% dis_set=result_data1([5:7,12:14],:)';%% 2221 4001 2222 2221 4001 2222
dis_set=result_data([5,6,12,14,7,13],:)';
RDT=result_data(18,:);
flag=0;
% close all;
plot_EBF_RDT_d1_d2(flag,EBF,RDT,dis_set);

input_type=types(central1_ID_nn_set,:);
%% filter conditions
RDT_tol=3.5;
cond_1=0.35;
cond_2=0.5;
valid_index1=find(RDT<RDT_tol & RDT>1.1);
valid_index2=find(dis_set(:,1)<cond_1 & dis_set(:,2)<cond_1 & dis_set(:,3)<cond_1...
    & dis_set(:,4)<cond_1 & dis_set(:,5)<cond_2 & dis_set(:,6)<cond_2);
%             valid_index2=find(dis_set(:,1)<cond_1 & dis_set(:,2)<cond_2 & dis_set(:,3)<cond_1...
%                 & dis_set(:,4)<cond_2 );
valid_index3=find(EBF>0);
i_s=intersect(valid_index1,valid_index2);
valid_index=intersect(i_s,valid_index3);

%%

output1=EBF(valid_index);
plot_EBF_RDT_d1_d2(flag,output1,RDT(valid_index),dis_set(valid_index,:));


input1=input_type(:,valid_index);
