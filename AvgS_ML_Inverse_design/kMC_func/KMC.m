function [t,mig_index, K_tot] = KMC(NN1_mig_energy_set,v0,T,t)
%KMC Summary of this function goes here
%   Detailed explanation goes here
        Kb = 8.6173e-5;
        K_set= v0.*exp(-(NN1_mig_energy_set)./(Kb*T));
        K_tot =sum(K_set);
        cum_k_set= cumsum(K_set);
        roulette_k_set = cum_k_set./K_tot;
        r1 =rand(1);
        mig_index = find(r1-roulette_k_set <0,1);
        r2 =rand(1);
        t = t + -1/K_tot* log(r2);
end

