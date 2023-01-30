function cross_boundary_set = cross_boundary_statistic(initial_pos,final_pos,boundary,cross_boundary_old_set)

dif_pos =final_pos-initial_pos;
fix_point=1e10;
dif_pos=round(dif_pos.*fix_point)./fix_point;
n= length(boundary);
cross_boundary_set=cross_boundary_old_set;
for num_dim = 1:n
    Len= boundary(num_dim,2)-boundary(num_dim,1);
    index1=find(dif_pos(:,num_dim)> Len/2);% left cross
    if ~isempty(index1)        
        cross_boundary_set(index1,num_dim)=cross_boundary_set(index1,num_dim)-1;
    end
    index2=find(dif_pos(:,num_dim)< -Len/2);
    if ~isempty(index2)
        cross_boundary_set(index2,num_dim)=cross_boundary_set(index2,num_dim)+1;
    end
end
end
