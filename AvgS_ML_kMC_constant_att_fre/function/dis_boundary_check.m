function output_dif_matrix= dis_boundary_check(input_dif_matrix, boundary)
% description

% input_matrix: 
   % size: m*3 ,m rows and 3(dimentions)
%boundary size is [3*2], in which the '2' means: [lowest_boundary,upper_boundary],3 means the dimensions;

   % output_matrix: same with the input_matrix
   
%%
n= length(boundary);
fix_point=1e10;
input_dif_matrix=round(input_dif_matrix.*fix_point)./fix_point;
for num_dim = 1:n
    Len= boundary(num_dim,2)-boundary(num_dim,1);
    index1=find(input_dif_matrix(:,num_dim)> Len/2);
    if ~isempty(index1)        
        input_dif_matrix(index1,num_dim)=input_dif_matrix(index1,num_dim)-Len;
    end
    index2=find(input_dif_matrix(:,num_dim)<-Len/2);
    if ~isempty(index2)
        input_dif_matrix(index2,num_dim)=input_dif_matrix(index2,num_dim)+Len;
    end
end
output_dif_matrix=input_dif_matrix;