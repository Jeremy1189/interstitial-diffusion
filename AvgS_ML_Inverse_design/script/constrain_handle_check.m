function constrain_violate = constrain_handle_check(pop)
[R,~]= size(pop);

%% we assume that '2' is the lower migration energy barrier, '1' is the higer migration energy barrier
% then [1-1-1] interstitial types have a higher energy barrier than '2'

% %% bounday check
% %lower
% lower_boundary = lu(1,:);
% M1=repmat(lower_boundary,R,1);
% 
% dif_lower = pop-M1;
% index_lower= find(dif_lower<0);
% pop(index_lower)=M1(index_lower);
% % upper
% upper_boundary = lu(2,:);
% M2=repmat(upper_boundary,R,1);
% dif_upper = pop-M2;
% index_upper= find(dif_upper>0);
% pop(index_upper)=M2(index_upper);
% other constraints
dif(:,1)= pop(:,2)-pop(:,1);
dif(:,2)= pop(:,4)-pop(:,3);
dif(:,3)= pop(:,6)-pop(:,5);
dif(:,4)= pop(:,8)-pop(:,7);
dif(:,5)= pop(:,1)-pop(:,8);
constrain_violate=zeros(R,5);
for i =1:R
    dif_cur = dif(i,:);
    dif_index = find(dif_cur<0,1);  
    if ~isempty(dif_index)
        constrain_violate(i,dif_index)=abs(dif_cur(dif_index));
%         pop(i,:) =  lu(1, :) + rand(1, L) .*  (lu(2, :) - lu(1, :));
          
    end
end

