function random_ratio_type_set= random_ratio_type(num_atoms,element_type_set,ratio_set)
%% input:
% %     num_atoms: total atoms
% %     element_type_set: [1,2,3,4,5] means five different elements 
% %    ratio_set: [0.1,0.1,0.3,0.3,0.2] the ratio of five elements 
%% output
% %   random_ratio_type_set
%%
        L=length(element_type_set);
        num_atoms_all=zeros(L,1);
        for i= 1:L-1
            num_atoms_all(i)= round(num_atoms*ratio_set(i));
        end
        num_atoms_all(L)= num_atoms-sum(num_atoms_all); % the number of the remain element
        cumsum_set = cumsum(num_atoms_all);
        per_type(1:cumsum_set(1))=element_type_set(1);
        for j = 2: L
            per_type(  (cumsum_set(j-1)+1) :cumsum_set(j))=element_type_set(j);
        end
        rand_index = randperm(num_atoms);
        random_ratio_type_set = per_type(rand_index);