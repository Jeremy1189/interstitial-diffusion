function Norm_X = normalization(X, var)

L = size(X, 1);

Lower = var(1, :);

Upper = var(2, :);

Norm_X = (X - repmat(Lower, L, 1)) ./ repmat(Upper - Lower, L, 1);