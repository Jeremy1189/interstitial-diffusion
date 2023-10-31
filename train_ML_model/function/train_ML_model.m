function [optimal_model,model_set_cell,min_mse,R2]=train_ML_model(input, output)
% input train

% output train
num_s=length(output);
num_kfold=10;
nodes=10;
indices= crossvalind('kfold',num_s,num_kfold);
R2_train=zeros(num_kfold,1);
R2_test=zeros(num_kfold,1);
mse=zeros(num_kfold,1);
bad_energy=[];
MAE_set=[];
test_input_set=[];
test_output_set=[];
bad_input_type=[];
test_predict_set=[];
for i = 1:num_kfold
    test_ind = (indices == i);
    train_ind = ~test_ind;
    test_input_data= input(:,test_ind);
    test_output_data=output(:,test_ind);
    train_input_data=input(:,train_ind);
    train_output_data=output(:,train_ind);
%      net = feedforwardnet(nodes,'trainlm');
    net = feedforwardnet(nodes,'trainbr');
%     net.layers{2}.transferFcn = 'tansig';
%     net = cascadeforwardnet(nodes,'trainbfg');
%     net = cascadeforwardnet(nodes,'trainbr');
%     net = patternnet(nodes,'trainbr');
    net.trainParam.showWindow=1;
    net.trainParam.epochs=500;
    net.trainParam.max_fail=13;
    net.divideParam.trainRatio=0.9;
    net.divideParam.valRatio=0.1;
    net.divideParam.testRatio=0;
    model_set= train(net,train_input_data, train_output_data,'useParallel','yes');
    train_predict_data= model_set( train_input_data);
    R2_train(i) = pearson_R2(train_predict_data,train_output_data);
    
    test_predict_data = model_set(test_input_data);
    model_set_cell{i}=model_set;
    R2_test(i) =pearson_R2(test_predict_data,test_output_data);
    MAE_current=test_predict_data-test_output_data;
    MAE_set=[MAE_set,MAE_current];
    test_input_set=[test_input_set,test_input_data];
    test_output_set=[test_output_set,test_output_data];
    test_predict_set=[test_predict_set,test_predict_data];
    index_bad= find(abs(MAE_current)>0.18);    
    index_bad_set{i}=index_bad;
    bad_energy=[bad_energy, test_output_data(index_bad)];
    bad_E_set{i}=[test_output_data(index_bad);test_predict_set(index_bad)];
    
    bad_input_type=[bad_input_type,test_input_data(:,index_bad)];
   
    mse(i)=1/length(test_output_data)*sum(abs(test_predict_data-test_output_data));
    
%     %
%     figure;
%     plot([0,1.5],[0,1.5],'--',"LineWidth",2)
%      hold on
%    scatter(test_predict_data(index_bad),test_output_data(index_bad));
%    xlabel('True migration energy (eV)')
%    ylabel('Prediction migration energy')
%        title(['CV number:',num2str(i)]);
% hold off
% set(gca,'FontSize',12,'FontName','Arial','FontWeight','bold','LineWidth',2)
end
[min_mse,min_index]= min(mse);
% optimal_model = model_set_cell{min_index};
 R2= R2_test(min_index);
    model_set= train(net,input, output,'useParallel','yes');
   predict_data= model_set(input);
    R2_model = pearson_R2(predict_data,output);
    optimal_model =model_set;