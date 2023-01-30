function offpop = DE(pop, NS, lu, n, F, CR)

offpop = zeros(NS, n);

for i = 1 : NS
    
    % Choose the indices for mutation
    indexSet = 1 : NS;
    indexSet(i) = [];
    
    % Choose the first Index
    temp = floor(rand * (NS - 1)) + 1;
    nouse(1) = indexSet(temp);
    indexSet(temp) = [];
    
    % Choose the second index
    temp = floor(rand * (NS - 2)) + 1;
    nouse(2) = indexSet(temp);
    indexSet(temp) = [];
    
    % Choose the third index
    temp = floor(rand * (NS - 3)) + 1;
    nouse(3) = indexSet(temp);
    
    % Generate the mutant vector
    %     V = pop(nouse(1), :) + F .* (pop(nouse(2), :) - pop(nouse(3), : ));
    V = pop(i, :) + rand .* (pop(nouse(1), :) - pop(i, : )) +  F .* (pop(nouse(2), :) - pop(nouse(3), :));
    
    % Handle the elements of the vector which violate the boundary
    vioLow = find(V < lu(1, : ));
    V(1, vioLow) = lu(1, vioLow);
    
    
    vioUpper = find(V > lu(2, : ));
    V(1, vioUpper) =  lu(2, vioUpper);
    
    % Implement the binomial crossover
    jRand = floor(rand * n) + 1;
    t = rand(1, n) < CR;
    t(1, jRand) = 1;
    t_ = 1 - t;
    U = t .* V + t_ .* pop(i, :);
    
    offpop(i,  : ) = U;
    
end
