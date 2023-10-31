function read_data = read_data_tmp(file_path,data_dim)
%current_path = C:\Users\86135\Desktop\NiFeCr_based pure Ni\FeNi\test_nn_effect;
count =0;
% data_result=zeros(12,num_samples);% the initial size of data
fid = fopen(file_path);
while ~feof(fid)
    tline = fgetl(fid);
    if ~isempty(tline)
        temp_data=sscanf(tline,'%f');
        
        if length(temp_data)==data_dim
            count= count +1;
            data_result(:,count)= temp_data';
%         elseif ~isempty(temp_data)
%             count= count +1;
%             data_result(:,count)= zeros(1,data_dim);
        else
            disp("Not a data line: ")
            disp(tline)
        end
    end
end
fclose(fid);
read_data=data_result;
end
