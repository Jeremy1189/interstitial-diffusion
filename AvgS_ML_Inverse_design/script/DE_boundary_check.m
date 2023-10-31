function pop=DE_boundary_check(pop,lu)
% if the pop exceed the lower boundary, it equals to the lower boundary; a
% same rule for the one out of the upper boundary, it equals to the upper
% boundary.
[R,~]= size(pop);
% lower
lower_boundary = lu(1,:);
M1=repmat(lower_boundary,R,1);

dif_lower = pop-M1;
index_lower= find(dif_lower<0);
pop(index_lower)=M1(index_lower);
% upper
upper_boundary = lu(2,:);
M2=repmat(upper_boundary,R,1);
dif_upper = pop-M2;
index_upper= find(dif_upper>0);
pop(index_upper)=M2(index_upper);

end