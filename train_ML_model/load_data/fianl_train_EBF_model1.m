
%% initial
clc;clear;close;
% data_path='/home/biaoxu4/interstitial/data_new_tol_5';
% addpath(genpath(data_path));
addpath(genpath(pwd));
%% load the input and output data
% load NiFe_interstitial_55.mat result type_set
% load data_interstitial1.mat result_data type_data
% load NiFeCr_interstitial_data.mat result_data_cell
load NiFe_interstitial_data8RDT0.6d10.8d2.mat result_data_cell 
result_data=cell2mat(result_data_cell(:,2));
type_data=cell2mat(result_data_cell(:,3));
result_data1=cell2mat(result_data_cell(:,4)');
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
NN_num1=7;
% NN_num2=NN_num1+2;
[central1_ID_nn_set,central1_ID_count_set,relative_central1_coords]...
    = update_interstitial_nn_kmc(central_position,per_coords,lattice_constant,supersize,NN_num1);

%% load the input and output data
perf_set=[];
I2_types = result_data(:,3);
types = [type_data,I2_types];
EBF = result_data(:,11);
dif_E=result_data(:,15)-result_data(:,14);
dis_set=result_data1([5:7,12:14],:)';%% 2221 4001 2222 2221 4001 2222
RDT=result_data(:,13);

input_type=types(:,central1_ID_nn_set);
%% filter conditions
for RDT_tol=3.5
    for cond_1=0.45
        for cond_2=0.5
            valid_index1=find(RDT<RDT_tol);
            valid_index2=find(dis_set(:,1)<cond_1 & dis_set(:,2)<cond_1 & dis_set(:,3)<cond_2...
                & dis_set(:,4)<cond_1 & dis_set(:,5)<cond_2 & dis_set(:,6)<cond_1);
            valid_index3=find(EBF>0);
            i_s=intersect(valid_index1,valid_index2);
            valid_index=intersect(i_s,valid_index3);
            
            %%            
            output=EBF(valid_index);
            output1=dif_E(valid_index);
            input=input_type(valid_index,:);
            N=length(output);
            input=input(1:N,:);
            output=output(1:N);
            % input= dec2bin_alloy3bit(input);
%             input= dec2bin_alloy3(input);
            input=input-1;
%             nodes=[10,8,4];
            % nodes=[10,10,7,4];
            % nodes=[12,8,4];
            nodes=10;
            % nodes=[50,30,20];
            % nodes=[12,6];
%             nodes=[12];
            % nodel
%             net = feedforwardnet(nodes,'trainscg');
            net = feedforwardnet(nodes,'trainbr');
            % net = cascadeforwardnet(nodes,'trainbr');
            % net = feedforwardnet(nodes,'traingdx');
            net.trainParam.showWindow=1;
            net.trainParam.epochs=500;
            net.trainParam.max_fail=8;
            net.divideParam.trainRatio=0.9;
            net.divideParam.valRatio=0.1;
            net.divideParam.testRatio=0;
            model_set_EBF= train(net,input', output','UseParallel','yes');
            model_set_FE= train(net,input', output1','UseParallel','yes');
%             model_set_EBF= train(net,input', output','UseGPU','yes');
            
            predict_data=model_set_EBF(input');
            R=corrcoef(predict_data,output);
            predict_data1=model_set_FE(input');
            R1=corrcoef(predict_data1,output1);            
            error_mae1=abs(predict_data1-output1');
            mae1=mean(error_mae1);
            
            figure;
            plot(1:length(predict_data),abs(predict_data-output'),'.');
            error_mae=abs(predict_data-output');
            mae=mean(error_mae);  
            index= find(error_mae>0.15);
            input_test= input_type(index,1:14);
            index_Fe_Cr=find(input_test(:,1)==2&input_test(:,2)==3);
            R_Fe_Cr=length(index_Fe_Cr)./length(input_test(:,1));
            perf_set=[perf_set;RDT_tol,cond_1,cond_2,R(2),mae,R1(2),mae1];
            file_name=[pwd,'/NiFeCr_data_',num2str(NN_num1),'nn',num2str(RDT_tol),'RDT',...
                num2str(cond_1),'d1',num2str(cond_2),'d2','.mat'];
            save(file_name, 'input', 'output', 'model_set_EBF','model_set_FE','R','mae','perf_set','R1','mae1')
        end
    end
end

