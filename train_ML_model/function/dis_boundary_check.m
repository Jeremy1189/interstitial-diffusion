function output_dif_matrix= dis_boundary_check(input_dif_matrix, boundary)
% description

% input_matrix: 
   % size: m*3 ,m rows and 3(dimentions)
%boundary size is [3*2], in which the '2' means: [lowest_boundary,upper_boundary],3 means the dimensions;

   % output_matrix: same with the input_matrix
   
%%
[~,n]= size(input_dif_matrix);
fix_point=1e10;
input_dif_matrix=round(input_dif_matrix.*fix_point)./fix_point;
for num_dim = 1:n
    Len= boundary(num_dim,2)-boundary(num_dim,1);
    index = find(abs(input_dif_matrix(:,num_dim))>Len/2);
    input_dif_matrix(index,num_dim)=Len-input_dif_matrix(index,num_dim);    
end
output_dif_matrix=input_dif_matrix;