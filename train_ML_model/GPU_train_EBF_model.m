
%% initial
clc;clear;close;
addpath(genpath(pwd));
%% load data_set1
load load_data/model_train_dif_ratio.mat result_data type_data
% % save model_train_dif_ratio.mat result_data type_data
% % result_data=cell2mat(result_data_cell(:,2));
% % type_data=cell2mat(result_data_cell(:,3));
% % result_data1=cell2mat(result_data_cell(:,4)');
dim=5;
tmp_per=read_data_tmp([pwd,'/load_data/tmp2'],dim);
% tmp_per=read_data_tmp([pwd,'/load_Data/per55'],dim);
tmp_ID=tmp_per(1,:);
[sort_v,sort_index]= sort(tmp_ID);
per=tmp_per(:,sort_index)';
% Create Interstitial
load per55; % perfect structure
I1_ID =2221;
I2_ID=4001;
I3_ID=2222;
set_dis=1.05;
lattice_constant=per(5,3);
supersize=10;
I1_coord=per(I1_ID,3:5);
I2_coord=per(I2_ID,3:5);
I3_coord=per(I3_ID,3:5);
disp([I1_coord;I2_coord;I3_coord])
I2_type=1;% it just for satisfy the structure of the per, and won't be used in the later calculation
% create the initial interstitial structure
per_coords=per(:,3:5);
per_ID_set = per(:,1);
% caluculation the nearest neibours, obtain the sort ID
central_position_mig_former= (I1_coord+I2_coord)./2;% the dumbbell central before migration
dif_vector =I1_coord-I2_coord;
L=lattice_constant*supersize;
boundary = [0,L;0,L;0,L];
dif_vector =dis_boundary_check(dif_vector,boundary);
central_position=I2_coord+0.5.*dif_vector;% the dumbbell central before migration
central_position=boundary_check(central_position,boundary);
RDT_tol=3.5;
cond_1=0.35;
cond_2=0.5;
% load the input and output data
perf_set=[];
I2_types = result_data(2,:);
types = [type_data;I2_types];
EBF = result_data(11,:);
dis_set=result_data([3,4,8,9],:)';
RDT=result_data(13,:);
flag=0;
for NN_num1=10
    
    
        [central1_ID_nn_set,central1_ID_count_set,relative_central1_coords]...
            = update_interstitial_nn_kmc(central_position,per_coords,lattice_constant,supersize,NN_num1);
%         close all;
        plot_EBF_RDT_d1_d2(flag,EBF,RDT,dis_set);

        input_type=types(central1_ID_nn_set,:);
        % filter conditions

        valid_index1=find(RDT<RDT_tol & RDT>1.1);
        valid_index2=find(dis_set(:,1)<cond_1 & dis_set(:,2)<cond_2 & dis_set(:,3)<cond_1...
            & dis_set(:,4)<cond_2 );
        valid_index3=find(EBF>0);
        i_s=intersect(valid_index1,valid_index2);
        valid_index=intersect(i_s,valid_index3);

        %
        output0=EBF(valid_index);
        plot_EBF_RDT_d1_d2(flag,output0,RDT(valid_index),dis_set(valid_index,:));
        input0=input_type(:,valid_index);

end
% input=input0;
% output=output0;

%% add new generated data: data_set2

new_data;

%%
input=[input0,input1];
output=[output0,output1];
% input=input1;
% output=output1;
% save train_data_10nn_RDT4.5_d1_035_d2_05_nodes3.mat input output nodes_num
% load train_data_10nn_RDT4.5_d1_035_d2_05_nodes3.mat input output nodes_num 
%% GPU train: finding the suitable trainning parameters 
N=length(output);
input=input(:,1:N);
output=output(1:N);
% % input= dec2bin_alloy3bit(input);
% %             input= dec2bin_alloy3(input);
for nodes_num=-3:2:20

nodes=[10,8,4]+nodes_num;
% nodes=[10,10,7,4];
% nodes=[12,8,4];
% nodes=10;
% nodes=[50,30,20];
% nodes=[12,6];
%             nodes=[12];
% nodel
net = feedforwardnet(nodes,'trainscg');
% net = feedforwardnet(nodes,'trainbr');
% net = cascadeforwardnet(nodes,'trainbr');
% net = feedforwardnet(nodes,'traingdx');
net.trainParam.showWindow=1;
net.trainParam.epochs=500;
net.trainParam.max_fail=6;
net.divideParam.trainRatio=0.9;
net.divideParam.valRatio=0.1;
net.divideParam.testRatio=0;
% model_set_EBF= train(net,input, output,'UseParallel','yes');
model_set_EBF= train(net,input, output,'UseGPU','yes');

predict_data=model_set_EBF(input);
R=corrcoef(predict_data,output);
error_mae=abs(predict_data-output);
mae=mean(error_mae);

%             figure;
%             plot(1:length(predict_data),abs(predict_data-output),'.');
%
%             index= find(error_mae>0.15);
%             input_test= input_type(index,1:14);
%             index_Fe_Cr=find(input_test(:,1)==2&input_test(:,2)==3);
%             R_Fe_Cr=length(index_Fe_Cr)./length(input_test(:,1));
perf_set=[perf_set;RDT_tol,cond_1,cond_2,R(2),mae];
file_name=[pwd,'/train_model_GPU/NiFe_final_',num2str(NN_num1),'nn',num2str(nodes_num),'nodes',num2str(RDT_tol),'RDT',...
    num2str(cond_1),'d1',num2str(cond_2),'d2','.mat'];
save(file_name, 'input', 'output', 'model_set_EBF','R','mae','perf_set')
end

