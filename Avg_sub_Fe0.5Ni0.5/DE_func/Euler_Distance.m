%% compute the Euler distance
% X1 : a row vector of D dimension
% X2 : a matrix of D dimension

function Euler_D = Euler_Distance(X1, X2)

L1 = size(X2, 1);

Euler_D = sqrt(sum((repmat(X1, L1, 1) - X2) .^2, 2));
