function mig_energy = mig_energy_determine (cur_mig_energy,Iterstitial_type1,Iterstitial_type2,Iterstitial_type3)
%% therer are 8 paths for the interstitial jumps, the different type..
%  interrrstitial corresponding to different energy , the detainl can be
%  found in the Table II of the "Ferasat, Keyvan, et al. "Accelerated...
%kinetic Monte Carlo: A case study; vacancy and dumbbell interstitial...
%diffusion traps in concentrated solid solution alloys." ...
%The Journal of Chemical Physics 153.7 (2020): 074109."
types_set=[ Iterstitial_type1,Iterstitial_type2,Iterstitial_type3];
matrix_set=[1,1,1;            
            1,1,2;
            1,2,1;            
            1,2,2;
            2,1,1;            
            2,1,2;
            2,2,1;
            2,2,2;];
 dif = repmat(types_set,length(cur_mig_energy),1)-matrix_set;
 sum_dif = sum(abs(dif),2);
 index=find(sum_dif==0);
 mig_energy=cur_mig_energy(index); %#ok<FNDSB>


end