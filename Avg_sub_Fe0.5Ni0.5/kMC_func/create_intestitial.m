function [interstitial_structure,int_ID_back,int_ID_mig] = create_intestitial(ratio,frac_coord_supercell,lattice_constant,boundary,dumbell_interval_dis)
%%
% ratio: 1X N, N is the number of elements, eg. Ni0.5Fe0.5 has a ratio [0.5,0.5]
% frac_coord_supercell: nX3, the fractional coordinates of the perfect supercell
%                       n is the number of atoms in perfect supercell 
% lattice constant: 1 dimensions, eg: 3.5195(Ni)
%
% output: 
%         interstitial_structure: (n+1)X5, the structure of created interstitial, added
%                                  an atom into the supercell at the centre
                                
%% perfect structure
num_atoms = size(frac_coord_supercell,1);
element_type_set =cumsum(ones(1,length(ratio)));
per_coords=frac_coord_supercell.*lattice_constant;
per_index=cumsum(ones(num_atoms,1));
per_type= random_ratio_type(length(per_index),element_type_set,ratio);
perfect_structure = [per_index,per_type',per_coords];

%% inset position of intestitial
supercell_centre = ones(1,3)./2;

dif_temp = frac_coord_supercell-repmat(supercell_centre,num_atoms,1);
[~,ID_centre_atoms_pos] = min(sum(abs(dif_temp),2));
% centre_atoms_pos = frac_coord_supercell(ID_centre_atoms_pos,:).*lattice_constant;

%% inset types of intestital
index_type= find(rand(1)<cumsum(ratio),1);
insert_interstitial_type = element_type_set(index_type);
%% interstitial structure
interstitial_structure=creat_dumbell_interstitial(perfect_structure,...
        ID_centre_atoms_pos, insert_interstitial_type,dumbell_interval_dis, boundary);

int_ID_back=interstitial_structure(end,1);
int_ID_mig=ID_centre_atoms_pos;

end

