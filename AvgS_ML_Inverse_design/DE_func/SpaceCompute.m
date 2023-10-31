function [Space_Sec, Space_D] = SpaceCompute(X)

L_samp = size(X, 1);

Space_D = zeros(L_samp, 1);

Space_Sec = zeros(L_samp, 1);

for i = 1 : L_samp
    
    Temp_Samp = X;
    
    X1 = X(i ,:);
    
    Temp_Samp(i, :) = [];
    
    X2 = Temp_Samp;
    
    Euler_D = Euler_Distance(X1, X2);
    
    Space_sort = sort(Euler_D);
    
    Space_D(i) = Space_sort(1);
    
    Space_Sec(i) = Space_sort(2);
    
end