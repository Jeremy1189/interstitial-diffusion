function modified_vector= boundary_check(input_vector,boundary)
%input
%input_vector：3 row 1 colomn vector,the check vector
% boundary: 3 row 2 colomn vector, [lowest_boundary_x, upper_boundary_x;...
%                                      lowest_boundary_y, upper_boundary_y;
%                                      lowest_boundary_z, upper_boundary_z];
% output
% the periodic boundary modified vactor
%  modified_vector
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
modified_vector=zeros(size(input_vector));
% lowest_boundary_x = boundary(1,1)-lattice_constant/4;
% lowest_boundary_y = boundary(2,1)-lattice_constant/4;
% lowest_boundary_z = boundary(3,1)-lattice_constant/4;
% upper_boundary_x=boundary(1,2)+lattice_constant/4;
% upper_boundary_y=boundary(2,2)+lattice_constant/4;
% upper_boundary_z=boundary(3,2)+lattice_constant/4;
lowest_boundary_x = boundary(1,1);
lowest_boundary_y = boundary(2,1);
lowest_boundary_z = boundary(3,1);
upper_boundary_x=boundary(1,2);
upper_boundary_y=boundary(2,2);
upper_boundary_z=boundary(3,2);
Len_x=boundary(1,2)-boundary(1,1);
Len_y=boundary(2,2)-boundary(2,1);
Len_z=boundary(3,2)-boundary(3,1);
input_vector=round(input_vector.*1e10)./1e10;
x=input_vector(1);
if x >=lowest_boundary_x && x<=upper_boundary_x
   modified_vector(1)=input_vector(1);
elseif x <lowest_boundary_x
   modified_vector(1)=input_vector(1)+Len_x;
elseif x >upper_boundary_x
    modified_vector(1)=input_vector(1)-Len_x;
else
    disp('x of the input_vector is an illegal input')
end

%%
y=input_vector(2);
if y >=lowest_boundary_y && y<=upper_boundary_x
   modified_vector(2)=input_vector(2);
elseif y <lowest_boundary_y
   modified_vector(2)=input_vector(2)+Len_y;
elseif y >upper_boundary_y
    modified_vector(2)=input_vector(2)-Len_y;
else
    disp('y of the input_vector is an illegal input')
end
%%
z=input_vector(3);
if z >=lowest_boundary_z && z<=upper_boundary_x
   modified_vector(3)=input_vector(3);
elseif z <lowest_boundary_z
   modified_vector(3)=input_vector(3)+Len_z;
elseif z >upper_boundary_z
    modified_vector(3)=input_vector(3)-Len_z;
else
    disp('z of the input_vector is an illegal input')
end
end

