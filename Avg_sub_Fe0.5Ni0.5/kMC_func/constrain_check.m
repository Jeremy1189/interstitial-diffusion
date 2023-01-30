function pop= constrain_check(pop,lu)
[R,L]= size(pop);
%% we assume that '2' is the lower migration energy barrier, '1' is the higer migration energy barrier
% then [1-1-1] interstitial types have a higher energy barrier than '2'

dif(:,1)= pop(:,2)-pop(:,1);
dif(:,2)= pop(:,4)-pop(:,3);
dif(:,3)= pop(:,6)-pop(:,5);
dif(:,4)= pop(:,8)-pop(:,7);
dif(:,5)= pop(:,1)-pop(:,8);
for i =1:R
    dif_cur = dif(i,:);
    dif_index = find(dif_cur<0,1);  
    if ~isempty(dif_index)
        pop(i,:) =  lu(1, :) + rand(1, L) .*  (lu(2, :) - lu(1, :));
    end
end

